// R4 gates (DESIGN.md section 3), the clean-metal figure's numbers:
//   Drude: the KPM spectral weight within the stated resolution of
//     omega = 0 matches the kinetic-energy Drude weight D = (pi/2N)(-K_x)
//     (drude_weight), |W[0, Delta] - D| < 0.045 (pin derivation below).
//   leakage: the regular-part weight above the stated resolution matches
//     the clean-metal zero within the leakage bound computed from the
//     Jackson kernel at the chosen n_cheb (the ridge device below).
//
// Config: L = 128, quarter filling (e_f mid-gap of the exact separable
// spectrum), n_cheb = 128, a = 4.04, stochastic trace R = 32 at the
// default seed — the deterministic trace costs N x n_cheb^2 SpMVs, so the
// stochastic path is mandatory at this L (STATE.md). Stated resolution
// Delta = 4 kernel half-widths = 4 pi a / n_cheb = 0.3968: the Jackson
// tail past 4 half-widths carries 0.5% of a Drude delta's weight
// (measured from the kernel by the ridge device), so "within resolution"
// covers the kernel's body. Fine grid 4096 over [0, 8] (dw = 0.002,
// 50x finer than a half-width — the R3 grid lesson).
//
// The Drude pin (the R3 L1 discipline: a regression bound pinned from
// measured levels, not an accuracy estimate; and per STATE.md nothing
// here pins from the DOS noise envelope — the stochastic current-current
// path runs hot, R2 measured linf 0.192 vs the envelope 0.0833 at L=6).
// Bring-up at the gate seed: W[0, Delta] = 0.493352, D = 0.514969
// (drude_weight at n_cheb = 1024, R = 32; vs the exact separable-
// spectrum D = 0.513196, rel err 3.5e-3 — the Jackson-smoothed Fermi
// edge against the finite-L level spacing). Offset 0.0216; pinned bound
// 0.045, ~2.1x measured. The offset is systematic, not noise: the true
// finite-L regular weight above Delta (exact Lehmann accounting on the
// separable spectrum, printed below) is 0.0179 and the kernel moves it
// out of the window. Seed spread at bring-up, 8 seeds: W[0, Delta] in
// [0.476, 0.496] (sigma ~ 0.008), worst |W - D| = 0.0389 < 0.045.
//
// The leakage gate fails against the pre-committed bound. Recorded
// plainly per the house rule; the bound is not adjusted. Measured
// regular weight at the gate seed 0.023225 vs the kernel-leakage bound
// 0.0026 — 9x over.
// Decomposition: true finite-L regular weight 0.0179 (exact) + kernel
// leakage 0.0026 + stochastic ~0.0007. The clean-metal zero is the
// thermodynamic-limit statement; the open-boundary lattice carries
// genuine finite-frequency weight — boundary-assisted absorption plus
// the near-Fermi transition clusters at omega_d ~ 4.4 d / L — decaying
// only as ~1/L (exact weight above omega = 0.4: 14.1% of the total at
// L = 32, 7.1% at 64, 3.5% at 128). Kill criterion 1 (DESIGN.md
// section 4): the honest DC output is the gated Drude weight; a finite
// DC number exists only as sigma(omega) at omega >= Delta, stated.
// Every check prints its measured value; any failure exits nonzero.

#include "kubo/kubo.hpp"

#include "kpm/lattice.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <numbers>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

constexpr double kDrudeBound = 0.045; // pinned, see the header comment

int failures = 0;

void check(bool ok, const char* name, const std::string& detail) {
  std::printf("[%s] %s: %s\n", ok ? "pass" : "FAIL", name, detail.c_str());
  if (!ok)
    ++failures;
}

// Exact 1D open-chain modes. The open square lattice is separable:
// E[p,q] = e[p] + e[q], and W (x bonds only) is diagonal in q:
// <p q|W|p' q'> = w[p,p'] delta_qq'. Validated against kubo_lehmann_ed
// at L=8 (total weight to 1e-15, bins to a bin-edge rounding flip).
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
double pick_ef(const Chain& c, int L, int& n_occ_out) {
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
  n_occ_out = n_occ;
  return 0.5 * (E[n_occ - 1] + E[n_occ]);
}

// Exact D = (pi/2N)(-K_x): -K_x = -sum_occ e[p], the x-part of the
// Kronecker sum.
double exact_drude(const Chain& c, int L, double e_f) {
  const int N = L * L;
  double K_x = 0.0;
  for (int p = 0; p < L; ++p)
    for (int q = 0; q < L; ++q)
      if (c.e[p] + c.e[q] < e_f)
        K_x += c.e[p];
  return (std::numbers::pi / (2.0 * N)) * (-K_x);
}

// Exact Lehmann pair weight above ostar: the finite-L "regular" content
// of the open lattice.
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

// Leakage measuring device: the clean metal's moment table is the
// diagonal ridge over all k of a periodic grid, mu_nm = (1/Np) sum_k
// v_k^2 T_n(e'_k) T_m(e'_k) with v_k = 2 sin k_x (t=1), occupation via
// the Fermi window at reconstruction. The reconstruction is the
// Jackson-broadened Drude delta; its weight above Delta is what the
// kernel alone leaks at this n_cheb. Device only — the measured side of
// the gate is the open-lattice table.
void ridge_leak(int nc, double a, double delta, double& above_out,
                double& total_out) {
  const int Lp = 128;
  const int Np = Lp * Lp;
  std::vector<double> ek(Np), v2(Np);
  for (int ix = 0; ix < Lp; ++ix)
    for (int iy = 0; iy < Lp; ++iy) {
      const double kx = 2.0 * std::numbers::pi * ix / Lp;
      const double ky = 2.0 * std::numbers::pi * iy / Lp;
      ek[ix * Lp + iy] = -2.0 * (std::cos(kx) + std::cos(ky));
      v2[ix * Lp + iy] = 4.0 * std::sin(kx) * std::sin(kx);
    }
  std::vector<double> es = ek;
  std::sort(es.begin(), es.end());
  const double e_f = 0.5 * (es[Np / 4 - 1] + es[Np / 4]);

  const int nc1 = nc + 1;
  kubo::MomentTable ridge;
  ridge.n_cheb = nc;
  ridge.mu.assign(static_cast<size_t>(nc1) * nc1, 0.0);
  for (int s = 0; s < Np; ++s) {
    const double ep = ek[s] / a;
    std::vector<double> T(nc1);
    T[0] = 1.0;
    T[1] = ep;
    for (int n = 1; n < nc; ++n)
      T[n + 1] = 2.0 * ep * T[n] - T[n - 1];
    const double c = v2[s] / Np;
    for (int n = 0; n < nc1; ++n) {
      const double cn = c * T[n];
      for (int m = 0; m < nc1; ++m)
        ridge.mu[static_cast<size_t>(n) * nc1 + m] += cn * T[m];
    }
  }
  kubo::KuboParams p;
  p.a = a;
  p.b = 0.0;
  p.n_cheb = nc;
  const double omega_max = 8.0;
  const int n_fine = 4096;
  const auto sig = kubo::optical_conductivity(ridge, p, e_f, omega_max, n_fine);
  const double dw = omega_max / n_fine;
  double above = 0.0, total = 0.0;
  for (const auto& [om, s] : sig) {
    total += s * dw;
    if (om > delta)
      above += s * dw;
  }
  above_out = above;
  total_out = total;
}

} // namespace

int main() {
  char buf[256];

  const int L = 128;
  const int N = L * L;
  const double a = 4.0 * 1.01; // kpm-tightbinding eta = 0.01 margin over 4
  const int n_cheb = 128;
  const double hw = std::numbers::pi * a / n_cheb; // kernel half-width
  const double Delta = 4.0 * hw;                   // stated resolution

  const Chain c = make_chain(L);
  int n_occ = 0;
  const double e_f = pick_ef(c, L, n_occ);
  const double D_exact = exact_drude(c, L, e_f);
  double W_true_total = 0.0;
  const double W_true_above = exact_above(c, L, e_f, Delta, W_true_total);
  std::printf("[info] config: L=%d N=%d quarter filling n_occ=%d e_f=%+.6f "
              "n_cheb=%d R=32 a=%.2f hw=%.4f Delta=4hw=%.4f\n",
              L, N, n_occ, e_f, n_cheb, a, hw, Delta);
  std::printf("[info] exact separable-spectrum: D=%.6f, true regular weight "
              "above Delta=%.6f (%.2f%% of total %.6f)\n",
              D_exact, W_true_above, 100.0 * W_true_above / W_true_total,
              W_true_total);

  const kubo::CsrMatrix H = kpm::square_lattice(L);
  const kubo::CsrMatrix W = kubo::velocity_x(L, 1.0);

  kubo::KuboParams pd;
  pd.a = a;
  pd.b = 0.0;
  pd.n_cheb = 1024;
  pd.n_random = 32;
  const double D = kubo::drude_weight(H, W, e_f, pd);
  std::printf("[info] drude_weight (n_cheb=%d R=%d) = %.6f vs exact %.6f, "
              "rel %.3e\n",
              pd.n_cheb, pd.n_random, D, D_exact, std::abs(D - D_exact) / D_exact);

  kubo::KuboParams p;
  p.a = a;
  p.b = 0.0;
  p.n_cheb = n_cheb; // fixed seed: the stochastic trace at large L
  const kubo::MomentTable mt = kubo::current_moments(H, W, p);
  const double mu00_want = 2.0 * (L - 1) / L;
  std::printf("[info] stochastic mu_00 = %.6f vs exact 2(L-1)/L = %.6f "
              "(hot-path indicator, not gated)\n",
              mt.mu[0], mu00_want);

  const double omega_max = 8.0;
  const int n_fine = 4096;
  const double dw = omega_max / n_fine;
  const auto sigma = kubo::optical_conductivity(mt, p, e_f, omega_max, n_fine);

  double W_below = 0.0, W_above = 0.0;
  for (const auto& [om, s] : sigma) {
    if (om <= Delta)
      W_below += s * dw;
    else
      W_above += s * dw;
  }

  const double drude_err = std::abs(W_below - D);
  std::snprintf(buf, sizeof buf,
                "W[0,Delta] = %.6f vs D = %.6f, |diff| = %.6f (pinned bound %.3f)",
                W_below, D, drude_err, kDrudeBound);
  check(drude_err < kDrudeBound,
        "R4 Drude: weight within resolution vs kinetic-energy Drude weight",
        buf);

  double leak_abs = 0.0, ridge_total = 0.0;
  ridge_leak(n_cheb, a, Delta, leak_abs, ridge_total);
  const double leak_bound = leak_abs * (D / ridge_total);
  std::snprintf(buf, sizeof buf,
                "W(Delta,8] = %.6f vs kernel-leakage bound %.6f "
                "(ridge device: %.6f of %.6f total past 4hw at n_cheb=%d)",
                W_above, leak_bound, leak_abs, ridge_total, n_cheb);
  check(W_above < leak_bound,
        "R4 leakage: regular-part weight vs the clean-metal zero", buf);

  std::printf("[info] decomposition of the regular weight: true finite-L "
              "%.6f + kernel leakage %.6f; stochastic seed spread at "
              "bring-up ~0.0007\n",
              W_true_above, leak_abs);
  std::printf("[info] f-sum consistency: W[0,Delta]+W(Delta,8] = %.6f vs "
              "D = %.6f\n",
              W_below + W_above, D);

  std::printf("%s (%d failure%s)\n",
              failures == 0 ? "ALL CHECKS PASSED" : "FAILURES", failures,
              failures == 1 ? "" : "s");
  return failures == 0 ? 0 : 1;
}
