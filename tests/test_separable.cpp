// Separable-spectrum accounting gate (correctness pass 2026-11-08).
//
// The R4 leakage decomposition and the clean-metal figure's exact
// overlay rest on the separable spectrum of the open square lattice:
// E[p,q] = e[p] + e[q] with the 1D open-chain modes
// e[p] = -2 cos(pi (p+1)/(L+1)), phi_p[j] = sqrt(2/(L+1))
// sin(pi (p+1)(j+1)/(L+1)), and W (x bonds only) diagonal in the
// transverse index: <p q|W|p' q'> = w[p,p'] delta_qq'. These are the
// standing-wave modes and current matrix elements of Rigol & Shastry,
// arXiv:0801.0592 — the reference for the recorded finite-L leakage
// physics (the Drude content of a finite open system lives in
// finite-frequency peaks at omega ~ 2 pi t l/L, l odd, collapsing to
// delta(0) only in the thermodynamic limit; D_L of an open system is
// identically zero at finite L).
//
// Gates at L = 8 against the R1 ED oracle (bounds 1e-12; values measured
// 2026-11-08):
//   spectrum: sorted chain E[p,q] vs ED eigenvalues, linf 8.9e-15;
//   total Lehmann pair weight: chain accounting vs kubo_lehmann_ed on
//     the CSR operators, |diff| 1.0e-15 (the recorded R4 validation
//     claim, re-verified);
//   Drude weight: chain exact D = (pi/2N)(-K_x) vs the kinetic-energy
//     expectation of the same ED eigensystem, |diff| 5.6e-17;
//   parity selection rule: w[p,p'] = 0 for even p-p' exactly in the sine
//     algebra (Rigol-Shastry), max |w| over even pairs 7.4e-16.
// Teeth verified by perturbation (scratch tree): periodic-mode chain
// energies turn the spectrum, total-weight, and Drude gates red.
//
// Info lines only, never gated: the exact regular weight above
// omega = 0.4 at fresh L (48, 96, 160 join the recorded 32/64/128/256)
// — the ~1/L decay data any future finite-L leakage bound (Option B,
// STATE.md 2026-11-08) must predict at fresh L before any re-gating.
// This suite asserts nothing about it; a fitted bound would be a
// description, not a gate.
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

// Exact 1D open-chain modes; w[a][b] = <a|W_1D|b> from the eigenvector
// sum (the same construction as tests/test_r4.cpp).
struct Chain {
  std::vector<double> e;
  std::vector<std::vector<double>> w;
};

Chain make_chain(int L) {
  Chain c;
  c.e.resize(L);
  std::vector<std::vector<double>> phi(L, std::vector<double>(L));
  for (int p = 0; p < L; ++p) {
    c.e[p] = -2.0 * std::cos(std::numbers::pi * (p + 1) / (L + 1));
    for (int j = 0; j < L; ++j)
      phi[p][j] = std::sqrt(2.0 / (L + 1)) *
                  std::sin(std::numbers::pi * (p + 1) * (j + 1) / (L + 1));
  }
  c.w.assign(L, std::vector<double>(L, 0.0));
  for (int a = 0; a < L; ++a)
    for (int b = 0; b < L; ++b) {
      double s = 0.0;
      for (int j = 0; j < L - 1; ++j)
        s += phi[a][j + 1] * phi[b][j] - phi[a][j] * phi[b][j + 1];
      c.w[a][b] = s;
    }
  return c;
}

// Quarter filling by the pick_n_occ rule (biggest gap near N/4) on the
// exact spectrum; returns e_f mid-gap.
double pick_ef(const Chain& c, int L) {
  const int N = L * L;
  std::vector<double> E(N);
  for (int p = 0; p < L; ++p)
    for (int q = 0; q < L; ++q)
      E[p * L + q] = c.e[p] + c.e[q];
  std::sort(E.begin(), E.end());
  const int lo = std::max(1, N / 4 - 3), hi = std::min(N - 1, N / 4 + 3);
  int n_occ = -1;
  double gap = 0.0;
  for (int k = lo; k <= hi; ++k)
    if (E[k] - E[k - 1] > gap) {
      gap = E[k] - E[k - 1];
      n_occ = k;
    }
  if (n_occ < 1)
    throw std::runtime_error("no gapped filling near N/4");
  return 0.5 * (E[n_occ - 1] + E[n_occ]);
}

double exact_drude(const Chain& c, int L, double e_f) {
  const int N = L * L;
  double K_x = 0.0;
  for (int p = 0; p < L; ++p)
    for (int q = 0; q < L; ++q)
      if (c.e[p] + c.e[q] < e_f)
        K_x += c.e[p];
  return (std::numbers::pi / (2.0 * N)) * (-K_x);
}

double exact_above(const Chain& c, int L, double e_f, double ostar,
                   double& total_out) {
  const int N = L * L;
  double above = 0.0, total = 0.0;
  for (int q = 0; q < L; ++q)
    for (int p = 0; p < L; ++p) {
      if (!(c.e[p] + c.e[q] < e_f))
        continue;
      for (int pp = 0; pp < L; ++pp) {
        if (!(c.e[pp] + c.e[q] > e_f))
          continue;
        const double om = c.e[pp] - c.e[p];
        const double wt =
            (std::numbers::pi / N) * c.w[p][pp] * c.w[p][pp] / om;
        total += wt;
        if (om > ostar)
          above += wt;
      }
    }
  total_out = total;
  return above;
}

} // namespace

int main() {
  char buf[256];
  const double pi = std::numbers::pi;

  const int L = 8;
  const int N = L * L;
  const Chain c = make_chain(L);
  const double e_f = pick_ef(c, L);

  const kubo::CsrMatrix H = kpm::square_lattice(L);
  const kubo::CsrMatrix W = kubo::velocity_x(L, 1.0);
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(dense(H));
  if (es.info() != Eigen::Success)
    throw std::runtime_error("test eigensolve failed");
  const Eigen::VectorXd& eval = es.eigenvalues();

  // Spectrum: sorted chain sums vs ED eigenvalues.
  {
    std::vector<double> E(N);
    for (int p = 0; p < L; ++p)
      for (int q = 0; q < L; ++q)
        E[p * L + q] = c.e[p] + c.e[q];
    std::sort(E.begin(), E.end());
    double linf = 0.0;
    for (int k = 0; k < N; ++k)
      linf = std::max(linf, std::abs(E[k] - eval[k]));
    std::snprintf(buf, sizeof buf, "L=%d linf(chain E[p+q], ED evals) = %.3e (bound 1e-12)",
                  L, linf);
    check(linf < 1e-12, "separable spectrum vs ED eigenvalues", buf);
  }

  // Parity selection rule: w[p,p'] = 0 for even p-p' (Rigol-Shastry).
  {
    double even_max = 0.0;
    for (int a = 0; a < L; ++a)
      for (int b = 0; b < L; ++b)
        if ((a - b) % 2 == 0)
          even_max = std::max(even_max, std::abs(c.w[a][b]));
    std::snprintf(buf, sizeof buf, "L=%d max|w[p,p']| over even p-p' = %.3e (bound 1e-12)",
                  L, even_max);
    check(even_max < 1e-12, "chain current matrix-element parity rule", buf);
  }

  // Total Lehmann pair weight: chain accounting vs the ED oracle on the
  // CSR operators. Binning-independent (every transition lies below 8).
  {
    double total_chain = 0.0;
    exact_above(c, L, e_f, 0.0, total_chain);
    const double omega_max = 8.0;
    const int n_grid = 4096;
    const auto sigma = kubo::kubo_lehmann_ed(H, W, e_f, omega_max, n_grid);
    const double dw = omega_max / n_grid;
    double total_ed = 0.0;
    for (const auto& [om, s] : sigma)
      total_ed += s * dw;
    std::snprintf(buf, sizeof buf,
                  "L=%d chain total weight %.15f vs ED %.15f, |diff| = %.3e (bound 1e-12)",
                  L, total_chain, total_ed, std::abs(total_chain - total_ed));
    check(std::abs(total_chain - total_ed) < 1e-12,
          "chain Lehmann total weight vs kubo_lehmann_ed", buf);
  }

  // Drude weight: chain exact vs the kinetic-energy expectation of the
  // same ED eigensystem (the R1 f-sum identity, both sides at once).
  {
    const Eigen::MatrixXd Hx = x_bond_part(H, L);
    double K_x = 0.0;
    for (int n = 0; n < N; ++n) {
      if (!(eval[n] < e_f))
        continue;
      K_x += es.eigenvectors().col(n).dot(Hx * es.eigenvectors().col(n));
    }
    const double D_ed = (pi / (2.0 * N)) * (-K_x);
    const double D_chain = exact_drude(c, L, e_f);
    std::snprintf(buf, sizeof buf,
                  "L=%d chain D %.15f vs ED kinetic %.15f, |diff| = %.3e (bound 1e-12)",
                  L, D_chain, D_ed, std::abs(D_chain - D_ed));
    check(std::abs(D_chain - D_ed) < 1e-12,
          "chain exact Drude weight vs ED kinetic energy", buf);
  }

  // Info only: the ~1/L decay of the exact regular weight above
  // omega = 0.4 — fresh L (48, 96, 160) alongside the recorded values.
  // Validation data for the owed Option B bound; this suite gates none
  // of it (a bound fitted to the run it gates is a description).
  std::printf("[info] exact regular weight above omega=0.4 (x L in parens):\n");
  for (int Lx : {32, 48, 64, 96, 128, 160, 256}) {
    const Chain cx = make_chain(Lx);
    const double efx = pick_ef(cx, Lx);
    double tot = 0.0;
    const double above = exact_above(cx, Lx, efx, 0.4, tot);
    std::printf("[info]   L=%3d: %.6f (%.2f%% of total %.6f), x L = %.4f\n", Lx,
                above, 100.0 * above / tot, tot, above * Lx);
  }

  std::printf("%s (%d failure%s)\n",
              failures == 0 ? "ALL CHECKS PASSED" : "FAILURES", failures,
              failures == 1 ? "" : "s");
  return failures == 0 ? 0 : 1;
}
