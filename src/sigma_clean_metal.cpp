// R4 figure data: the clean-metal sigma_1(omega) at large L, quarter
// filling — the KPM reconstruction (fine grid) and the exact Lehmann
// spectrum of the separable open lattice (binned) at the same e_f, plus
// the gated Drude weight. Defaults are the test_r4 gate config, so the
// figure is the gated object; deterministic for a fixed seed.

#include "kubo/kubo.hpp"

#include "kpm/lattice.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <numbers>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

// Exact 1D open-chain modes: E[p,q] = e[p] + e[q], <p q|W|p' q'> =
// w[p,p'] delta_qq' (derivation in tests/test_r4.cpp).
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

void write_csv(const std::string& path,
               const std::vector<std::pair<double, double>>& xy) {
  std::ofstream f(path);
  if (!f)
    throw std::runtime_error("cannot open " + path);
  f << "omega,sigma\n";
  f.precision(17);
  for (const auto& [om, s] : xy)
    f << om << ',' << s << '\n';
  if (!f)
    throw std::runtime_error("write failed: " + path);
}

} // namespace

int main(int argc, char** argv) {
  int L = 128;
  kubo::KuboParams p;
  p.a = 4.0 * 1.01;
  p.b = 0.0;
  p.n_cheb = 128;
  p.n_random = 32;
  int n_grid = 4096;
  for (int i = 1; i < argc; ++i) {
    try {
      if (std::strncmp(argv[i], "--L=", 4) == 0)
        L = std::stoi(argv[i] + 4);
      else if (std::strncmp(argv[i], "--n-cheb=", 9) == 0)
        p.n_cheb = std::stoi(argv[i] + 9);
      else if (std::strncmp(argv[i], "--R=", 4) == 0)
        p.n_random = std::stoi(argv[i] + 4);
      else if (std::strncmp(argv[i], "--seed=", 7) == 0)
        p.seed = std::stoull(argv[i] + 7);
      else if (std::strncmp(argv[i], "--grid=", 7) == 0)
        n_grid = std::stoi(argv[i] + 7);
      else
        throw std::invalid_argument(std::string("unknown argument: ") + argv[i]);
    } catch (const std::exception& e) {
      std::fprintf(stderr, "sigma_clean_metal: %s\n", e.what());
      std::fprintf(stderr, "usage: sigma_clean_metal [--L=128] [--n-cheb=128] "
                           "[--R=32] [--seed=N] [--grid=4096]\n");
      return 1;
    }
  }

  try {
    const double omega_max = 8.0;
    const Chain c = make_chain(L);
    int n_occ = 0;
    const double e_f = pick_ef(c, L, n_occ);
    const kubo::CsrMatrix H = kpm::square_lattice(L);
    const kubo::CsrMatrix W = kubo::velocity_x(L, 1.0);

    const auto t0 = std::chrono::steady_clock::now();
    const kubo::MomentTable mt = kubo::current_moments(H, W, p);
    const auto t1 = std::chrono::steady_clock::now();
    const auto sigma = kubo::optical_conductivity(mt, p, e_f, omega_max, n_grid);
    const auto t2 = std::chrono::steady_clock::now();

    kubo::KuboParams pd;
    pd.a = p.a;
    pd.b = p.b;
    pd.n_cheb = 1024;
    pd.n_random = 32;
    const double D = kubo::drude_weight(H, W, e_f, pd);

    // Exact Lehmann spectrum of the same open lattice, binned coarse
    // enough that the near-Fermi transition clusters show.
    const int n_bin = 256;
    const double dwb = omega_max / n_bin;
    std::vector<std::pair<double, double>> exact(n_bin);
    for (int k = 0; k < n_bin; ++k)
      exact[k] = {(k + 0.5) * dwb, 0.0};
    const int N = L * L;
    for (int q = 0; q < L; ++q)
      for (int pi_ = 0; pi_ < L; ++pi_) {
        if (!(c.e[pi_] + c.e[q] < e_f))
          continue;
        for (int pp = 0; pp < L; ++pp) {
          if (!(c.e[pp] + c.e[q] > e_f))
            continue;
          const double om = c.e[pp] - c.e[pi_];
          const int k = static_cast<int>(std::floor(om / dwb));
          if (k < n_bin)
            exact[k].second +=
                (std::numbers::pi / N) * c.w[pi_][pp] * c.w[pi_][pp] / om / dwb;
        }
      }

    const std::string base = std::string(KUBO_DATA_DIR) +
                             "/sigma_clean_metal_L" + std::to_string(L);
    write_csv(base + ".csv", sigma);
    write_csv(base + "_exact.csv", exact);

    const double dwf = omega_max / n_grid;
    const double Delta = 4.0 * std::numbers::pi * p.a / p.n_cheb;
    double below = 0.0, above = 0.0;
    for (const auto& [om, s] : sigma) {
      if (om <= Delta)
        below += s * dwf;
      else
        above += s * dwf;
    }
    std::printf("clean metal L=%d N=%d n_occ=%d e_f=%+.6f n_cheb=%d R=%d "
                "seed=%llu\n",
                L, N, n_occ, e_f, p.n_cheb, p.n_random,
                (unsigned long long)p.seed);
    std::printf("Drude weight D = %.6f, stated resolution Delta = 4hw = %.4f: "
                "W[0,Delta] = %.6f, W(Delta,%.1f] = %.6f\n",
                D, Delta, below, omega_max, above);
    std::printf("moments %.1f s, reconstruction %.1f s -> %s.csv (+ _exact)\n",
                std::chrono::duration<double>(t1 - t0).count(),
                std::chrono::duration<double>(t2 - t1).count(), base.c_str());
  } catch (const std::exception& e) {
    std::fprintf(stderr, "sigma_clean_metal failed: %s\n", e.what());
    return 1;
  }
  return 0;
}
