// R3 gates (DESIGN.md section 3), at the pre-committed tolerances:
//   f-sum re-gated on the reconstruction: int sigma_1 d-omega of the KPM
//     reconstruction vs (pi/2N)(-<K_x>) from the same eigensystem,
//     relative error < 1% (pre-committed; Jackson damping redistributes
//     weight, the integral is the conserved thing).
//   reconstruction vs the R1 Lehmann oracle at L = 8: binned L1 within
//     the regression bound pinned below.
//
// The L1 pin (the kpm-tightbinding L1 discipline: a regression bound,
// not an accuracy estimate). Trace policy: deterministic — R2 measured
// the stochastic R=32 current-current table running hot (linf 0.192 vs
// the DOS-calibrated envelope 4/sqrt(2NR) = 0.0833 at L=6) while the
// deterministic table matched ED at linf 1.725e-13, so this gate runs
// with n_random = N = 64 (exact unit-basis trace) and the pipeline is
// seed-independent: the L1 is a fixed number per toolchain, dominated
// by Jackson smoothing against the finite-L discreteness of the oracle,
// not by noise. Bring-up measurement at these settings (L=8,
// n_cheb=512, a=4.04, n_grid=64 bins over [0,8], the KPM reconstruction
// bin-averaged over 32 sub-samples per bin — the matched-bin discipline
// of kpm-tightbinding's test_kpm.cpp): L1 = 0.0775. Pinned bound:
// 0.16, ~2x measured, recorded before the gate was treated as
// authoritative; never repinned around an observed outcome.
// Every check prints its measured value; any failure exits nonzero.

#include "kubo/kubo.hpp"

#include "kpm/lattice.hpp"

#include <Eigen/Dense>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <numbers>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

constexpr double kL1Bound = 0.16; // pinned, see the header comment

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

// H restricted to its x bonds (horizontal neighbors, same lattice row):
// the kinetic-energy operator K_x entering the f-sum rule.
Eigen::MatrixXd x_bond_part(const kubo::CsrMatrix& H, int L) {
  Eigen::MatrixXd Hx = Eigen::MatrixXd::Zero(H.n, H.n);
  for (int row = 0; row < H.n; ++row)
    for (int p = H.row_ptr[row]; p < H.row_ptr[row + 1]; ++p) {
      const int col = H.col_idx[p];
      if (std::abs(row - col) == 1 && row / L == col / L)
        Hx(row, col) = H.values[p];
    }
  return Hx;
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

// Occupied-state count near quarter filling with a clear gap at e_f.
int pick_n_occ(const Eigen::VectorXd& evals) {
  const int N = static_cast<int>(evals.size());
  int n_occ = -1;
  double gap = 0.0;
  const int lo = std::max(1, N / 4 - 3);
  const int hi = std::min(N - 1, N / 4 + 3);
  for (int k = lo; k <= hi; ++k) {
    const double g = evals[k] - evals[k - 1];
    if (g > gap) {
      gap = g;
      n_occ = k;
    }
  }
  if (n_occ < 1 || gap < 1e-6)
    throw std::runtime_error("no gapped filling near N/4");
  return n_occ;
}

} // namespace

int main() {
  char buf[256];
  const double pi = std::numbers::pi;

  const int L = 8;
  const int N = L * L;
  const kubo::CsrMatrix H = kpm::square_lattice(L);
  const kubo::CsrMatrix W = kubo::velocity_x(L, 1.0);
  const Eigensystem es = eigensolve(H);
  const int n_occ = pick_n_occ(es.evals);
  const double e_f = 0.5 * (es.evals[n_occ - 1] + es.evals[n_occ]);
  std::printf("[info] L=%d n_occ=%d e_f=%+.6f (gap [%.6f, %.6f])\n", L, n_occ,
              e_f, es.evals[n_occ - 1], es.evals[n_occ]);

  // f-sum reference: (pi/2N)(-<K_x>) over the occupied sea, from the
  // same eigensystem (the R1 identity).
  const Eigen::MatrixXd Hx = x_bond_part(H, L);
  double K_x = 0.0;
  for (int n = 0; n < n_occ; ++n) {
    const Eigen::VectorXd hxv = Hx * es.evecs.col(n);
    K_x += es.evecs.col(n).dot(hxv);
  }
  const double fsum = (pi / (2.0 * N)) * (-K_x);

  kubo::KuboParams p;
  p.a = 4.0 * 1.01; // kpm-tightbinding eta = 0.01 margin over bound 4
  p.b = 0.0;
  p.n_cheb = 512;
  p.n_random = N; // deterministic construction: exact unit-basis trace
  const kubo::MomentTable mt = kubo::current_moments(H, W, p);

  const double omega_max = 8.0;
  const int n_grid = 64;
  const int sub = 32;
  const int n_fine = n_grid * sub;
  const auto sigma_kpm = kubo::optical_conductivity(mt, p, e_f, omega_max, n_fine);
  const auto sigma_ed = kubo::kubo_lehmann_ed(H, W, e_f, omega_max, n_grid);

  const double dw = omega_max / n_grid;
  const double dwf = omega_max / n_fine;

  double l1 = 0.0;
  for (int k = 0; k < n_grid; ++k) {
    double avg = 0.0;
    for (int s = 0; s < sub; ++s)
      avg += sigma_kpm[k * sub + s].second;
    avg /= sub;
    l1 += std::abs(avg - sigma_ed[k].second) * dw;
  }
  std::snprintf(buf, sizeof buf,
                "L=%d n_cheb=%d n_grid=%d (x%d sub-samples) L1 = %.6f (pinned bound %.2f)",
                L, p.n_cheb, n_grid, sub, l1, kL1Bound);
  check(l1 < kL1Bound, "R3 reconstruction vs Lehmann oracle, binned L1", buf);

  double total = 0.0, smin = 0.0;
  for (const auto& [om, s] : sigma_kpm) {
    (void)om;
    total += s * dwf;
    smin = std::min(smin, s);
  }
  std::snprintf(buf, sizeof buf,
                "int sigma_1 d-omega = %.15f vs (pi/2N)(-K_x) = %.15f rel = %.3e (bound 1e-2)",
                total, fsum, std::abs(total - fsum) / fsum);
  check(std::abs(total - fsum) / fsum < 1e-2,
        "R3 f-sum re-gated on the reconstruction", buf);

  std::printf("[info] min sigma_kpm = %.3e over %d fine points\n", smin, n_fine);

  // Positivity, gated (correctness pass 2026-11-08): the Jackson-damped
  // reconstruction is nonnegative by construction — the product kernel is
  // positive and the table's density is positive — and STATE.md records
  // sigma >= 0 exactly on the fine grid; measured min exactly 0.000e+00
  // here. Teeth verified by perturbation: stripping the Jackson damping
  // (g_n := 1) drives min sigma to -2.078 (Gibbs) while this gate's
  // sibling f-sum goes red at 2.0e-2.
  std::snprintf(buf, sizeof buf, "min sigma_kpm = %.3e over %d fine points",
                smin, n_fine);
  check(smin >= 0.0, "R3 sigma(omega) >= 0 on the fine grid", buf);

  std::printf("%s (%d failure%s)\n",
              failures == 0 ? "ALL CHECKS PASSED" : "FAILURES", failures,
              failures == 1 ? "" : "s");
  return failures == 0 ? 0 : 1;
}
