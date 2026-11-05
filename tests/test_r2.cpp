// R2 gate (DESIGN.md section 3), at the pre-committed tolerance:
//   current-current Chebyshev moment table mu_nm (Weisse et al., RMP 78,
//   275 (2006), arXiv:cond-mat/0504627, Eq. 155: mu_nm = (1/N) Tr[v_x
//   T_n(H') v_x T_m(H')]) at L = 6 with the deterministic construction —
//   n_random >= N selects the exact unit-basis trace — against the
//   ED-derived table
//     mu_nm = (1/N) sum_ab |<a|v_x|b>|^2 T_n(e'_a) T_m(e'_b)
//   at matched rescaling, linf < 1e-12.
// Plus an exact mu_00 anchor, mu_00 = (1/N) ||W||_F^2 = 2(L-1)/L at
// t = 1. The stochastic path (n_random = 32 < N) is measured against the
// same ED table and printed next to the kpm-tightbinding noise envelope
// 4/sqrt(2*N*R) — information for R3's regression bound, not a gate
// (the 4-sigma discipline applies only if the stochastic trace is used
// for the gate; it is not). Every check prints its measured value; any
// failure exits nonzero.

#include "kubo/kubo.hpp"

#include "kpm/lattice.hpp"

#include <Eigen/Dense>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

int failures = 0;

void check(bool ok, const char* name, const std::string& detail) {
  std::printf("[%s] %s: %s\n", ok ? "pass" : "FAIL", name, detail.c_str());
  if (!ok)
    ++failures;
}

Eigen::MatrixXd dense(const kubo::CsrMatrix& A) {
  Eigen::MatrixXd M = Eigen::MatrixXd::Zero(A.n, A.n);
  for (int row = 0; row < A.n; ++row)
    for (int p = A.row_ptr[row]; p < A.row_ptr[row + 1]; ++p)
      M(row, A.col_idx[p]) = A.values[p];
  return M;
}

struct Eigensystem {
  Eigen::VectorXd evals;
  Eigen::MatrixXd evecs;
};

Eigensystem eigensolve(const kubo::CsrMatrix& H) {
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(dense(H));
  if (es.info() != Eigen::Success)
    throw std::runtime_error("test eigensolve failed");
  return {es.eigenvalues(), es.eigenvectors()};
}

} // namespace

int main() {
  char buf[256];

  const int L = 6;
  const int N = L * L;
  const kubo::CsrMatrix H = kpm::square_lattice(L);
  const kubo::CsrMatrix W = kubo::velocity_x(L, 1.0);

  kubo::KuboParams p; // n_cheb = 512, fixed seed; pin the rescaling
  p.a = 4.0 * 1.01;   // kpm-tightbinding eta = 0.01 margin over bound 4
  p.b = 0.0;
  p.n_random = N; // deterministic construction: exact unit-basis trace

  const kubo::MomentTable mt = kubo::current_moments(H, W, p);
  const int nc1 = mt.n_cheb + 1;
  std::snprintf(buf, sizeof buf, "n_cheb=%d mu.size=%zu want %d", mt.n_cheb,
                mt.mu.size(), nc1 * nc1);
  check(mt.n_cheb == p.n_cheb &&
            mt.mu.size() == static_cast<size_t>(nc1) * nc1,
        "MomentTable shape", buf);

  // ED-derived deterministic table at matched rescaling. T_n by the
  // recurrence, never cos(n arccos e) (house rule, kpm-tightbinding
  // docs/KPM.md).
  const Eigensystem es = eigensolve(H);
  const Eigen::MatrixXd WP = dense(W) * es.evecs;
  const Eigen::MatrixXd w = es.evecs.transpose() * WP;
  std::vector<double> tn(static_cast<size_t>(N) * nc1);
  for (int a = 0; a < N; ++a) {
    const double e = (es.evals[a] - p.b) / p.a;
    tn[static_cast<size_t>(a) * nc1] = 1.0;
    tn[static_cast<size_t>(a) * nc1 + 1] = e;
    for (int k = 1; k < mt.n_cheb; ++k)
      tn[static_cast<size_t>(a) * nc1 + k + 1] =
          2.0 * e * tn[static_cast<size_t>(a) * nc1 + k] -
          tn[static_cast<size_t>(a) * nc1 + k - 1];
  }
  std::vector<double> mu_ed(static_cast<size_t>(nc1) * nc1, 0.0);
  for (int a = 0; a < N; ++a)
    for (int b = 0; b < N; ++b) {
      const double c = w(a, b) * w(a, b) / N;
      for (int nn = 0; nn < nc1; ++nn) {
        const double cn = c * tn[static_cast<size_t>(a) * nc1 + nn];
        for (int mm = 0; mm < nc1; ++mm)
          mu_ed[static_cast<size_t>(nn) * nc1 + mm] +=
              cn * tn[static_cast<size_t>(b) * nc1 + mm];
      }
    }

  double linf = 0.0;
  int arg_n = 0, arg_m = 0;
  for (int nn = 0; nn < nc1; ++nn)
    for (int mm = 0; mm < nc1; ++mm) {
      const double d =
          std::abs(mt.mu[static_cast<size_t>(nn) * nc1 + mm] -
                   mu_ed[static_cast<size_t>(nn) * nc1 + mm]);
      if (d > linf) {
        linf = d;
        arg_n = nn;
        arg_m = mm;
      }
    }
  std::snprintf(buf, sizeof buf,
                "L=%d n_cheb=%d linf = %.3e (bound 1e-12) at (n,m)=(%d,%d)",
                L, mt.n_cheb, linf, arg_n, arg_m);
  check(linf < 1e-12, "R2 moment table vs ED-derived table", buf);

  // Exact anchor: mu_00 = (1/N) sum_ij W_ij^2 = 2L(L-1)/N at t = 1.
  const double mu00_want = 2.0 * (L - 1) / L;
  std::snprintf(buf, sizeof buf, "mu_00 = %.15f want 2(L-1)/L = %.15f err %.3e",
                mt.mu[0], mu00_want, std::abs(mt.mu[0] - mu00_want));
  check(std::abs(mt.mu[0] - mu00_want) < 1e-12,
        "mu_00 = (1/N) ||W||_F^2 anchor", buf);

  double asym = 0.0;
  for (int nn = 0; nn < nc1; ++nn)
    for (int mm = 0; mm < nc1; ++mm)
      asym = std::max(asym, std::abs(mt.mu[static_cast<size_t>(nn) * nc1 + mm] -
                                     mt.mu[static_cast<size_t>(mm) * nc1 + nn]));
  std::printf("[info] table symmetry max|mu_nm - mu_mn| = %.3e\n", asym);

  // Stochastic path at the default R = 32, same seed: noise level vs the
  // kpm-tightbinding envelope, measured for R3's regression bound.
  kubo::KuboParams ps = p;
  ps.n_random = 32;
  const kubo::MomentTable mts = kubo::current_moments(H, W, ps);
  double linf_s = 0.0;
  for (size_t i = 0; i < mu_ed.size(); ++i)
    linf_s = std::max(linf_s, std::abs(mts.mu[i] - mu_ed[i]));
  const double env = 4.0 / std::sqrt(2.0 * N * ps.n_random);
  std::printf("[info] stochastic R=%d: linf vs ED = %.3e (4/sqrt(2NR) = %.4f)\n",
              ps.n_random, linf_s, env);

  std::printf("%s (%d failure%s)\n",
              failures == 0 ? "ALL CHECKS PASSED" : "FAILURES", failures,
              failures == 1 ? "" : "s");
  return failures == 0 ? 0 : 1;
}
