#pragma once

#include "kpm/lattice.hpp"

#include <cstdint>
#include <stdexcept>
#include <utility>
#include <vector>

namespace kubo {

/// CSR sparse matrix, FP64. DESIGN.md section 6 decision, taken at R1:
/// kpm-tightbinding is a build dependency pinned at a commit (see
/// CMakeLists.txt), so this is kpm::CsrMatrix itself, not a copy.
using kpm::CsrMatrix;

/// Moment-run parameters. a = 0 means unset: the caller must supply a
/// spectral bound; the R2 gate pins a = 4*1.01, the kpm-tightbinding
/// safety margin (eta = 0.01) over the square-lattice bound 4t.
struct KuboParams {
  int n_cheb = 512;
  int n_random = 32;
  double a = 0.0;
  double b = 0.0;
  std::uint64_t seed = 0x9E3779B9ULL;
};

/// 2D Chebyshev moment table mu_nm, row-major, (n_cheb+1)^2 entries.
struct MomentTable {
  int n_cheb = 0;
  std::vector<double> mu;
};

/// R1: velocity operator v_x = i[H, x] on the LxL lattice with hopping t
/// (H as in kpm::square_lattice: -t on nearest-neighbor bonds, open
/// boundaries, site index s = i*L + j, x_s = j). v_x is purely imaginary
/// in the site basis, so the real CSR stores W = [H, x] = -i v_x, real
/// antisymmetric: W_{ss'} = H_{ss'} (x_{s'} - x_s), nonzero only between
/// x-neighbors, W = -t along +x. Every observable in this module is
/// quadratic in v_x, and |<n|v_x|m>|^2 = <n|W|m>^2. Lattice constant 1.
CsrMatrix velocity_x(int L, double t);

/// R1: Kubo-Greenwood sigma_1(omega) in the Lehmann representation from
/// dense ED of H with velocity V (as stored by velocity_x), at Fermi
/// energy e_f. Spinless, e = hbar = 1, N = H.n sites:
///   sigma_1(omega) = (pi/N) sum_{n occ, m unocc}
///                    |<n|v_x|m>|^2 / (e_m - e_n) delta(omega - (e_m-e_n)).
/// Returned binned on n_grid uniform bins over [0, omega_max], values at
/// bin centers: sum_k sigma_k (omega_max/n_grid) is exactly the total
/// pair weight inside the window (peaks above omega_max are dropped, so
/// weight conservation wants omega_max >= the spectral width). States
/// exactly at e_f count as unoccupied. Oracle for all later rungs.
std::vector<std::pair<double, double>>
kubo_lehmann_ed(const CsrMatrix& H, const CsrMatrix& V, double e_f,
                double omega_max, int n_grid);

/// R2: current-current Chebyshev moment table of H and the velocity V
/// (as stored by velocity_x), Weisse et al., RMP 78, 275 (2006),
/// arXiv:cond-mat/0504627, Eq. 155 (fetched and verified against the
/// arXiv text at R2):
///   mu_nm = (1/N) Tr[v_x T_n(H') v_x T_m(H')],   H' = (H - b)/a,
/// the moments of the current matrix-element density j(x, y) (their
/// Eq. 153) that R3 reconstructs sigma(omega) from. With v_x = iW and
/// real antisymmetric W the estimator stays real and positive on the
/// diagonal:
///   mu_nm = (1/N) sum_r <W r | T_n(H') | W T_m(H') r>
///         = (1/N) sum_ab |<a|v_x|b>|^2 T_n(e'_a) T_m(e'_b),
/// symmetric in n, m. Trace evaluation: exact over the N unit vectors
/// once n_random >= N (the R2 gate uses this deterministic
/// construction); below that, the stochastic average over n_random
/// random +/-1 vectors (their Eq. 40, seeded). Row n, column m of the
/// table sits at mu[n*(n_cheb+1) + m], n, m = 0..n_cheb.
MomentTable current_moments(const CsrMatrix& H, const CsrMatrix& V,
                            const KuboParams& p);

/// R3: optical conductivity sigma_1(omega) at Fermi energy e_f from the
/// R2 moment table, Weisse et al. Eqs. 154 and 156-157 evaluated at
/// T = 0 (verified against the arXiv text at R3):
///   sigma_1(omega) = (pi/omega) int_{e_f-omega}^{e_f} j(x, x+omega) dx,
/// in the house pi-convention of kubo_lehmann_ed (Weisse carry their own
/// prefactor; the R1 oracle is the authority this is gated against). j is
/// the Jackson-damped 2D Chebyshev reconstruction of the table (Eq. 154,
/// h_nm = h_n h_m with h_0 = 1, h_n = 2 per Eq. 95), evaluated on a
/// uniform midpoint grid in rescaled energy and integrated over the
/// Fermi window by midpoint quadrature with fractional edge cells. The
/// table stays raw; the Jackson kernel (kpm::jackson_kernel) is applied
/// here at reconstruction, as in kpm-tightbinding's kpm_dos. p supplies
/// the rescaling and must be the one the table was built with (only a, b
/// are read). Returned on n_grid uniform bins over [0, omega_max],
/// values at bin centers — the same convention as kubo_lehmann_ed, so
/// the two compare bin by bin.
std::vector<std::pair<double, double>>
optical_conductivity(const MomentTable& table, const KuboParams& p,
                     double e_f, double omega_max, int n_grid);

/// R4: the kinetic-energy Drude weight D = (pi/2N)(-<K_x>) of the clean
/// metal at Fermi energy e_f — the only honest "DC" output (DESIGN.md
/// section 2). K_x is H restricted to the x bonds, recovered as the
/// support of V as stored by velocity_x; the occupied-sea trace
/// Tr[K_x theta(e_f - H)] is evaluated from stochastic Chebyshev moments
/// of the (signed) K_x-weighted DOS, Jackson-damped and integrated to
/// e_f on a fine midpoint grid. Cost is O(n_cheb) SpMVs per random
/// vector — the DOS price, not the 2D-table price. Accuracy is set by
/// the Jackson-smoothed Fermi edge against the finite-L level spacing:
/// measured rel err 3.5e-3 at L=128, n_cheb=1024, R=32 vs the exact
/// separable-spectrum value (test_r4). Deterministic for a fixed seed.
double drude_weight(const CsrMatrix& H, const CsrMatrix& V, double e_f,
                    const KuboParams& p);

} // namespace kubo
