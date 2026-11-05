#include "kubo/kubo.hpp"

#include "kpm/kpm.hpp"

#include <Eigen/Dense>

#include <algorithm>
#include <cmath>
#include <numbers>
#include <random>

namespace kubo {

namespace {

Eigen::MatrixXd dense(const CsrMatrix& A) {
  Eigen::MatrixXd M = Eigen::MatrixXd::Zero(A.n, A.n);
  for (int row = 0; row < A.n; ++row)
    for (int p = A.row_ptr[row]; p < A.row_ptr[row + 1]; ++p)
      M(row, A.col_idx[p]) = A.values[p];
  return M;
}

// y = (A x - b x) / a, the rescaled Hamiltonian applied to x (same idiom
// as kpm-tightbinding's src/kpm.cpp).
void spmv_scaled(const CsrMatrix& A, double a, double b,
                 const std::vector<double>& x, std::vector<double>& y) {
  for (int row = 0; row < A.n; ++row) {
    double acc = 0.0;
    for (int p = A.row_ptr[row]; p < A.row_ptr[row + 1]; ++p)
      acc += A.values[p] * x[A.col_idx[p]];
    y[row] = (acc - b * x[row]) / a;
  }
}

double dot(const std::vector<double>& u, const std::vector<double>& v) {
  double acc = 0.0;
  for (size_t i = 0; i < u.size(); ++i)
    acc += u[i] * v[i];
  return acc;
}

} // namespace

CsrMatrix velocity_x(int L, double t) {
  if (L < 2)
    throw std::invalid_argument("velocity_x: L must be >= 2");
  CsrMatrix W;
  W.n = L * L;
  W.row_ptr.reserve(W.n + 1);
  W.row_ptr.push_back(0);
  W.col_idx.reserve(2 * L * (L - 1));
  W.values.reserve(2 * L * (L - 1));
  for (int i = 0; i < L; ++i) {
    for (int j = 0; j < L; ++j) {
      const int s = i * L + j;
      if (j > 0) {
        W.col_idx.push_back(s - 1);
        W.values.push_back(t);
      }
      if (j < L - 1) {
        W.col_idx.push_back(s + 1);
        W.values.push_back(-t);
      }
      W.row_ptr.push_back(static_cast<int>(W.col_idx.size()));
    }
  }
  return W;
}

std::vector<std::pair<double, double>>
kubo_lehmann_ed(const CsrMatrix& H, const CsrMatrix& V, double e_f,
                double omega_max, int n_grid) {
  if (H.n != V.n)
    throw std::invalid_argument("kubo_lehmann_ed: H and V dimension mismatch");
  if (H.n < 1)
    throw std::invalid_argument("kubo_lehmann_ed: empty H");
  if (n_grid < 1 || omega_max <= 0.0)
    throw std::invalid_argument("kubo_lehmann_ed: bad omega grid");

  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(dense(H));
  if (es.info() != Eigen::Success)
    throw std::runtime_error("kubo_lehmann_ed: eigensolve failed");
  const Eigen::VectorXd& eval = es.eigenvalues();
  // w_nm = <n|W|m> is real since W is real antisymmetric and the
  // eigenvectors of real symmetric H are real; |<n|v_x|m>|^2 = w_nm^2.
  // Split into two products: the nested triple product trips a GCC 16
  // -Wmaybe-uninitialized false positive inside Eigen's evaluators.
  const Eigen::MatrixXd WP = dense(V) * es.eigenvectors();
  const Eigen::MatrixXd w = es.eigenvectors().transpose() * WP;

  const double dw = omega_max / n_grid;
  std::vector<std::pair<double, double>> sigma(n_grid);
  for (int k = 0; k < n_grid; ++k)
    sigma[k] = {(k + 0.5) * dw, 0.0};

  const double pref = std::numbers::pi / H.n;
  for (int n = 0; n < H.n; ++n) {
    if (!(eval[n] < e_f))
      continue;
    for (int m = 0; m < H.n; ++m) {
      if (!(eval[m] > e_f))
        continue;
      const double omega = eval[m] - eval[n];
      const double weight = pref * w(n, m) * w(n, m) / omega;
      const int k = static_cast<int>(std::floor(omega / dw));
      if (k < n_grid)
        sigma[k].second += weight / dw;
    }
  }
  return sigma;
}

MomentTable current_moments(const CsrMatrix& H, const CsrMatrix& V,
                            const KuboParams& p) {
  if (H.n != V.n)
    throw std::invalid_argument("current_moments: H and V dimension mismatch");
  if (H.n < 1)
    throw std::invalid_argument("current_moments: empty H");
  if (p.n_cheb < 1)
    throw std::invalid_argument("current_moments: n_cheb must be >= 1");
  if (p.n_random < 1)
    throw std::invalid_argument("current_moments: n_random must be >= 1");
  if (p.a <= 0.0)
    throw std::invalid_argument(
        "current_moments: a unset, caller supplies the spectral bound");

  const int n = H.n;
  const int nc = p.n_cheb;
  const int nc1 = nc + 1;
  // Exact trace over the unit basis once n_random reaches the dimension
  // (cheap at the R2 sizes and the sharper gate); below that, the
  // stochastic trace over random +/-1 vectors.
  const bool exact = p.n_random >= n;
  const int n_vec = exact ? n : p.n_random;

  MomentTable table;
  table.n_cheb = nc;
  table.mu.assign(static_cast<size_t>(nc1) * nc1, 0.0);

  std::mt19937_64 rng(p.seed);
  std::vector<double> r(n), vr(n), vm(n);
  std::vector<double> t_prev(n), t_cur(n), t_next(n);
  std::vector<double> d_prev(n), d_cur(n), d_next(n);
  for (int v = 0; v < n_vec; ++v) {
    if (exact) {
      std::fill(r.begin(), r.end(), 0.0);
      r[v] = 1.0;
    } else {
      for (int i = 0; i < n; ++i)
        r[i] = (rng() & 1ULL) ? 1.0 : -1.0;
    }
    kpm::spmv(V, r, vr);
    // Column m of the table gets <V r, T_n(H') V T_m(H') r> for all n.
    auto accumulate = [&](int m, const std::vector<double>& tm) {
      kpm::spmv(V, tm, vm);
      table.mu[m] += dot(vr, vm); // n = 0
      spmv_scaled(H, p.a, p.b, vm, d_cur);
      table.mu[nc1 + m] += dot(vr, d_cur); // n = 1
      d_prev = vm;
      for (int k = 1; k < nc; ++k) {
        spmv_scaled(H, p.a, p.b, d_cur, d_next);
        for (int i = 0; i < n; ++i)
          d_next[i] = 2.0 * d_next[i] - d_prev[i];
        table.mu[static_cast<size_t>(k + 1) * nc1 + m] += dot(vr, d_next);
        std::swap(d_prev, d_cur);
        std::swap(d_cur, d_next);
      }
    };
    // Outer recurrence T_m(H') r, same shape as kpm::chebyshev_moments.
    t_prev = r;
    spmv_scaled(H, p.a, p.b, t_prev, t_cur);
    accumulate(0, t_prev);
    accumulate(1, t_cur);
    for (int m = 1; m < nc; ++m) {
      spmv_scaled(H, p.a, p.b, t_cur, t_next);
      for (int i = 0; i < n; ++i)
        t_next[i] = 2.0 * t_next[i] - t_prev[i];
      accumulate(m + 1, t_next);
      std::swap(t_prev, t_cur);
      std::swap(t_cur, t_next);
    }
  }
  const double scale =
      exact ? 1.0 / n : 1.0 / (static_cast<double>(n) * n_vec);
  for (double& x : table.mu)
    x *= scale;
  return table;
}

std::vector<std::pair<double, double>>
optical_conductivity(const MomentTable& table, const KuboParams& p,
                     double e_f, double omega_max, int n_grid) {
  const int nc = table.n_cheb;
  const int nc1 = nc + 1;
  if (nc < 1 || table.mu.size() != static_cast<size_t>(nc1) * nc1)
    throw std::invalid_argument("optical_conductivity: bad moment table");
  if (p.a <= 0.0)
    throw std::invalid_argument(
        "optical_conductivity: a unset, must match the table's rescaling");
  if (n_grid < 1 || omega_max <= 0.0)
    throw std::invalid_argument("optical_conductivity: bad omega grid");

  // Jackson damping and the h_nm normalization fold into weights
  // w_nm = h_n h_m g_n g_m mu_nm; the table itself stays raw.
  const std::vector<double> g = kpm::jackson_kernel(nc);
  std::vector<double> w(static_cast<size_t>(nc1) * nc1);
  for (int n = 0; n < nc1; ++n) {
    const double hn = n == 0 ? 1.0 : 2.0;
    for (int m = 0; m < nc1; ++m) {
      const double hm = m == 0 ? 1.0 : 2.0;
      w[static_cast<size_t>(n) * nc1 + m] =
          table.mu[static_cast<size_t>(n) * nc1 + m] * g[n] * g[m] * hn * hm;
    }
  }

  // j~ = F / (pi^2 sqrt(1-eps^2) sqrt(1-eta^2)) with
  // F(eps, eta) = sum_nm w_nm T_n(eps) T_m(eta), evaluated on a uniform
  // midpoint grid in rescaled energy; two chained contractions keep it
  // O(nc1 * G^2) instead of O(nc1^2 * G^2). T_m by the recurrence, never
  // cos(m arccos eps) (house rule).
  const int G = 4 * nc1;
  const double de = 2.0 / G;
  std::vector<double> tg(static_cast<size_t>(nc1) * G);
  for (int j = 0; j < G; ++j) {
    const double eps = -1.0 + (j + 0.5) * de;
    tg[j] = 1.0;
    tg[G + j] = eps;
    for (int m = 1; m < nc; ++m)
      tg[static_cast<size_t>(m + 1) * G + j] =
          2.0 * eps * tg[static_cast<size_t>(m) * G + j] -
          tg[static_cast<size_t>(m - 1) * G + j];
  }
  std::vector<double> u(static_cast<size_t>(nc1) * G, 0.0);
  for (int n = 0; n < nc1; ++n)
    for (int m = 0; m < nc1; ++m) {
      const double wnm = w[static_cast<size_t>(n) * nc1 + m];
      for (int j = 0; j < G; ++j)
        u[static_cast<size_t>(n) * G + j] +=
            wnm * tg[static_cast<size_t>(m) * G + j];
    }
  std::vector<double> dens(static_cast<size_t>(G) * G, 0.0);
  for (int n = 0; n < nc1; ++n)
    for (int i = 0; i < G; ++i) {
      const double tni = tg[static_cast<size_t>(n) * G + i];
      for (int j = 0; j < G; ++j)
        dens[static_cast<size_t>(i) * G + j] +=
            tni * u[static_cast<size_t>(n) * G + j];
    }

  // sigma_1(omega) = 1/(pi a^2 omega') int_{e_f'-omega'}^{e_f'} d-eps
  // F(eps, eps+omega') / (sqrt(1-eps^2) sqrt(1-(eps+omega')^2));
  // midpoint quadrature over the Fermi window with fractional edge
  // cells, F linearly interpolated in eta (smooth on the grid scale).
  const double ef = (e_f - p.b) / p.a;
  const double dw = omega_max / n_grid;
  std::vector<std::pair<double, double>> sigma(n_grid);
  for (int k = 0; k < n_grid; ++k) {
    const double om = (k + 0.5) * dw;
    const double omp = om / p.a;
    const double lo = ef - omp, hi = ef;
    const int i0 = std::max(0, static_cast<int>(std::floor((lo + 1.0) / de)));
    const int i1 = std::min(G - 1, static_cast<int>(std::floor((hi + 1.0) / de)));
    double acc = 0.0;
    for (int i = i0; i <= i1; ++i) {
      const double cell_lo = -1.0 + i * de;
      const double cov =
          std::min(cell_lo + de, hi) - std::max(cell_lo, lo);
      if (cov <= 0.0)
        continue;
      const double eps = cell_lo + 0.5 * de;
      const double eta = eps + omp;
      if (eta >= 1.0 - 0.5 * de)
        continue; // the support of j~ ends at the rescaled band edge 1
      const double pos = (eta + 1.0) / de - 0.5;
      const int j0 = static_cast<int>(std::floor(pos));
      const double frac = pos - j0;
      const double fv = dens[static_cast<size_t>(i) * G + j0] * (1.0 - frac) +
                        dens[static_cast<size_t>(i) * G + j0 + 1] * frac;
      acc += cov * fv / std::sqrt((1.0 - eps * eps) * (1.0 - eta * eta));
    }
    sigma[k] = {om, acc / (std::numbers::pi * p.a * om)};
  }
  return sigma;
}

double drude_weight(const CsrMatrix& H, const CsrMatrix& V, double e_f,
                    const KuboParams& p) {
  if (H.n != V.n)
    throw std::invalid_argument("drude_weight: H and V dimension mismatch");
  if (H.n < 1)
    throw std::invalid_argument("drude_weight: empty H");
  if (p.n_cheb < 1)
    throw std::invalid_argument("drude_weight: n_cheb must be >= 1");
  if (p.n_random < 1)
    throw std::invalid_argument("drude_weight: n_random must be >= 1");
  if (p.a <= 0.0)
    throw std::invalid_argument(
        "drude_weight: a unset, caller supplies the spectral bound");

  const int n = H.n;
  // K_x = H restricted to the support of V (the x bonds). H's rows are
  // not column-sorted (kpm::square_lattice pushes the y bonds first), so
  // intersect by scan, not merge.
  CsrMatrix Kx;
  Kx.n = n;
  Kx.row_ptr.assign(n + 1, 0);
  Kx.col_idx.reserve(V.nnz());
  Kx.values.reserve(V.nnz());
  for (int row = 0; row < n; ++row) {
    for (int ph = H.row_ptr[row]; ph < H.row_ptr[row + 1]; ++ph)
      for (int pv = V.row_ptr[row]; pv < V.row_ptr[row + 1]; ++pv)
        if (V.col_idx[pv] == H.col_idx[ph]) {
          Kx.col_idx.push_back(H.col_idx[ph]);
          Kx.values.push_back(H.values[ph]);
          break;
        }
    Kx.row_ptr[row + 1] = static_cast<int>(Kx.col_idx.size());
  }

  // Moments of the K_x-weighted DOS: mu_n = (1/N)(1/R) sum_r <K_x r,
  // T_n(H') r>. K_x is symmetric, so no adjoint is needed.
  const int nc = p.n_cheb;
  std::vector<double> mu(nc + 1, 0.0);
  std::mt19937_64 rng(p.seed);
  std::vector<double> r(n), kr(n), t_prev(n), t_cur(n), t_next(n);
  for (int v = 0; v < p.n_random; ++v) {
    for (int i = 0; i < n; ++i)
      r[i] = (rng() & 1ULL) ? 1.0 : -1.0;
    kpm::spmv(Kx, r, kr);
    mu[0] += dot(kr, r);
    spmv_scaled(H, p.a, p.b, r, t_cur);
    mu[1] += dot(kr, t_cur);
    t_prev = r;
    for (int k = 1; k < nc; ++k) {
      spmv_scaled(H, p.a, p.b, t_cur, t_next);
      for (int i = 0; i < n; ++i)
        t_next[i] = 2.0 * t_next[i] - t_prev[i];
      mu[k + 1] += dot(kr, t_next);
      std::swap(t_prev, t_cur);
      std::swap(t_cur, t_next);
    }
  }
  const double scale = 1.0 / (static_cast<double>(n) * p.n_random);
  for (double& x : mu)
    x *= scale;

  // Jackson-damped reconstruction of the K_x-weighted DOS, integrated to
  // e_f' on a uniform midpoint grid with a fractional edge cell (the R3
  // grid lesson). D = (pi/2N)(-K_x) with K_x the per-site integral, so
  // N cancels.
  const std::vector<double> g = kpm::jackson_kernel(nc);
  const int G = 4 * (nc + 1);
  const double de = 2.0 / G;
  const double ef = (e_f - p.b) / p.a;
  double kx_per_site = 0.0;
  for (int j = 0; j < G; ++j) {
    const double cell_lo = -1.0 + j * de;
    const double cov = std::min(cell_lo + de, ef) - cell_lo;
    if (cov <= 0.0)
      continue;
    const double eps = cell_lo + 0.5 * de;
    double t0 = 1.0, t1 = eps;
    double acc = g[0] * mu[0] + 2.0 * g[1] * mu[1] * t1;
    for (int k = 1; k < nc; ++k) {
      const double t2 = 2.0 * eps * t1 - t0;
      acc += 2.0 * g[k + 1] * mu[k + 1] * t2;
      t0 = t1;
      t1 = t2;
    }
    kx_per_site += cov * acc / (std::numbers::pi * std::sqrt(1.0 - eps * eps));
  }
  return -0.5 * std::numbers::pi * kx_per_site;
}

} // namespace kubo
