#include "kubo/kubo.hpp"

#include <Eigen/Dense>

#include <cmath>
#include <numbers>

namespace kubo {

namespace {

Eigen::MatrixXd dense(const CsrMatrix& A) {
  Eigen::MatrixXd M = Eigen::MatrixXd::Zero(A.n, A.n);
  for (int row = 0; row < A.n; ++row)
    for (int p = A.row_ptr[row]; p < A.row_ptr[row + 1]; ++p)
      M(row, A.col_idx[p]) = A.values[p];
  return M;
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

} // namespace kubo
