// R1 gates (DESIGN.md section 3), at the pre-committed tolerances:
//   f-sum rule: total Lehmann weight vs the kinetic-energy expectation
//     from the same eigensystem, relative error < 1e-12. The identity
//     behind it: sum over n occ, m unocc of |v_nm|^2/omega_mn equals
//     -(1/2) <K_x>, since [x,[H,x]] = -H_x and occupied-occupied pairs
//     cancel; with sigma_1 as declared in kubo.hpp that reads
//     integral sigma_1 d-omega = (pi/2N) (-<K_x>).
//   sigma(omega) >= 0 at every grid point.
//   closure: the sum over final states reconstructs the <v_x^2>
//     diagonal to 1e-12.
// Plus the section 2 cross-gate (ED exact DOS vs the committed
// kpm-tightbinding reference CSVs, L = 8..14) and an exact L=2 anchor:
// one occupied state, a single transition at omega = 2 with weight pi/8.
// Every check prints its measured value; any failure exits nonzero.

#include "kubo/kubo.hpp"

#include "kpm/exact_oracle.hpp"
#include "kpm/histogram.hpp"
#include "kpm/lattice.hpp"

#include <Eigen/Dense>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <numbers>
#include <stdexcept>
#include <string>
#include <utility>
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

std::vector<std::pair<double, double>> read_dos_csv(const std::string& path) {
  std::ifstream f(path);
  if (!f)
    throw std::runtime_error("cannot open " + path);
  std::vector<std::pair<double, double>> out;
  std::string line;
  std::getline(f, line); // header
  while (std::getline(f, line)) {
    if (line.empty())
      continue;
    const std::size_t comma = line.find(',');
    out.emplace_back(std::stod(line.substr(0, comma)),
                     std::stod(line.substr(comma + 1)));
  }
  return out;
}

} // namespace

int main() {
  char buf[256];
  const double pi = std::numbers::pi;

  // velocity_x is W = [H, x]: dense commutator check against the
  // definition, plus nnz and antisymmetry, exact in integer arithmetic.
  for (int L : {2, 4, 6, 8}) {
    const kubo::CsrMatrix H = kpm::square_lattice(L);
    const kubo::CsrMatrix W = kubo::velocity_x(L, 1.0);
    const int N = L * L;
    Eigen::VectorXd x(N);
    for (int s = 0; s < N; ++s)
      x[s] = s % L;
    const Eigen::MatrixXd comm =
        dense(H) * x.asDiagonal() - x.asDiagonal() * dense(H);
    const Eigen::MatrixXd Wd = dense(W);
    const double def_err = (Wd - comm).cwiseAbs().maxCoeff();
    const double asym_err = (Wd + Wd.transpose()).cwiseAbs().maxCoeff();
    std::snprintf(buf, sizeof buf, "L=%d nnz=%d (want %d) max|W-[H,x]|=%.3e max|W+W^T|=%.3e",
                  L, W.nnz(), 2 * L * (L - 1), def_err, asym_err);
    check(W.nnz() == 2 * L * (L - 1) && def_err == 0.0 && asym_err == 0.0,
          "velocity_x is exactly [H,x], antisymmetric, x-bond sparse", buf);
  }

  // L=2 anchor: spectrum {-2,0,0,2}, one occupied state at e_f = -1.
  // W psi_0 lands entirely in the zero eigenspace: a single transition
  // at omega = 2 with ||W psi_0||^2 = 1, hence weight (pi/N)(1/2) = pi/8.
  {
    const int L = 2;
    const kubo::CsrMatrix H = kpm::square_lattice(L);
    const kubo::CsrMatrix W = kubo::velocity_x(L, 1.0);
    const double omega_max = 8.0;
    const int n_grid = 1024;
    const double dw = omega_max / n_grid;
    const auto sigma = kubo::kubo_lehmann_ed(H, W, -1.0, omega_max, n_grid);
    double total = 0.0, smin = 0.0, peak_w = 0.0, peak_wc = 0.0;
    for (const auto& [wc, s] : sigma) {
      total += s * dw;
      smin = std::min(smin, s);
      if (s * dw > peak_w) {
        peak_w = s * dw;
        peak_wc = wc;
      }
    }
    // Numerical dust: the degenerate zero eigenspace is rotated
    // arbitrarily, leaving ~1e-33 of weight at omega = 4. Require the
    // remainder outside the peak bin to be <= 1e-24.
    const double dust = total - peak_w;
    const double expect = pi / 8.0;
    std::snprintf(buf, sizeof buf,
                  "total weight=%.15f want pi/8=%.15f rel=%.3e dust=%.3e",
                  total, expect, std::abs(total - expect) / expect, dust);
    check(std::abs(total - expect) / expect < 1e-12 && dust <= 1e-24,
          "L=2 anchor: single transition, weight pi/8", buf);
    std::snprintf(buf, sizeof buf, "peak bin center omega=%.6f (bin width %.6f), weight=%.15f",
                  peak_wc, dw, peak_w);
    check(std::abs(peak_wc - 2.0) < dw && std::abs(peak_w - expect) / expect < 1e-12,
          "L=2 anchor: transition at omega=2", buf);
    std::snprintf(buf, sizeof buf, "min sigma = %.3e", smin);
    check(smin >= 0.0, "L=2 sigma >= 0", buf);
  }

  // f-sum, positivity, closure at L = 4, 6, 8.
  for (int L : {4, 6, 8}) {
    const int N = L * L;
    const kubo::CsrMatrix H = kpm::square_lattice(L);
    const kubo::CsrMatrix W = kubo::velocity_x(L, 1.0);
    const Eigensystem es = eigensolve(H);
    const int n_occ = pick_n_occ(es.evals);
    const double e_f = 0.5 * (es.evals[n_occ - 1] + es.evals[n_occ]);

    const double omega_max = 8.0;
    const int n_grid = 1024;
    const double dw = omega_max / n_grid;
    const auto sigma = kubo::kubo_lehmann_ed(H, W, e_f, omega_max, n_grid);

    double total = 0.0, smin = 0.0;
    for (const auto& [wc, s] : sigma) {
      (void)wc;
      total += s * dw;
      smin = std::min(smin, s);
    }

    // Kinetic energy of the x bonds over the occupied sea, from the
    // same eigensystem.
    const Eigen::MatrixXd Hx = x_bond_part(H, L);
    double K_x = 0.0;
    for (int n = 0; n < n_occ; ++n)
      K_x += es.evecs.col(n).dot(Hx * es.evecs.col(n));
    const double fsum = (pi / (2.0 * N)) * (-K_x);
    const double rel = std::abs(total - fsum) / std::abs(fsum);
    std::snprintf(buf, sizeof buf,
                  "L=%d n_occ=%d e_f=%+.6f weight=%.15f (pi/2N)(-K_x)=%.15f rel=%.3e",
                  L, n_occ, e_f, total, fsum, rel);
    check(rel < 1e-12, "f-sum: Lehmann weight vs kinetic energy", buf);

    std::snprintf(buf, sizeof buf, "L=%d min sigma = %.3e over %d bins", L, smin, n_grid);
    check(smin >= 0.0, "sigma(omega) >= 0 on the grid", buf);

    // Closure: for every eigenstate n, sum_m |<n|W|m>|^2 reconstructs
    // <n|v_x^2|n> = ||W phi_n||^2.
    const Eigen::MatrixXd WP = dense(W) * es.evecs;
    const Eigen::MatrixXd w = es.evecs.transpose() * WP;
    double closure = 0.0;
    for (int n = 0; n < N; ++n)
      closure = std::max(closure,
                         std::abs(w.row(n).squaredNorm() - WP.col(n).squaredNorm()));
    std::snprintf(buf, sizeof buf, "L=%d max_n |sum_m w_nm^2 - ||W phi_n||^2| = %.3e",
                  L, closure);
    check(closure < 1e-12, "closure: <v_x^2> diagonal from the final-state sum", buf);
  }

  // Cross-gate (DESIGN.md section 2): the H entering the module is the
  // verified one — ED exact DOS vs the committed kpm-tightbinding CSVs.
  for (int L = 8; L <= 14; ++L) {
    const auto ref = read_dos_csv(std::string(KPM_TB_DATA_DIR) + "/exact_dos_L" +
                                  std::to_string(L) + ".csv");
    const auto mine =
        kpm::dos_histogram(kpm::exact_eigenvalues(kpm::square_lattice(L)), 0.1);
    if (ref.size() != mine.size())
      throw std::runtime_error("reference CSV size mismatch at L=" + std::to_string(L));
    double linf = 0.0;
    for (std::size_t k = 0; k < ref.size(); ++k)
      linf = std::max({linf, std::abs(ref[k].first - mine[k].first),
                       std::abs(ref[k].second - mine[k].second)});
    std::snprintf(buf, sizeof buf, "L=%2d linf vs exact_dos_L%d.csv = %.3e", L, L, linf);
    check(linf < 1e-12, "exact-DOS cross-gate vs kpm-tightbinding reference", buf);
  }

  std::printf("%s (%d failure%s)\n",
              failures == 0 ? "ALL CHECKS PASSED" : "FAILURES", failures,
              failures == 1 ? "" : "s");
  return failures == 0 ? 0 : 1;
}
