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
/// spectral bound (4t covers the square lattice at mu = 0); the first
/// moment rung (R2) pins the value actually used.
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

/// R2: current-current Chebyshev moment table of H and V (construction
/// per Weisse et al., verified there before use).
inline MomentTable current_moments(const CsrMatrix&, const CsrMatrix&,
                                   const KuboParams&) {
  throw std::logic_error("current_moments: not implemented (DESIGN.md R2)");
}

/// R3: optical conductivity on an omega grid from the moment table at
/// Fermi energy e_f.
inline std::vector<std::pair<double, double>>
optical_conductivity(const MomentTable&, double) {
  throw std::logic_error("optical_conductivity: not implemented (DESIGN.md R3)");
}

/// R4: Drude weight of the clean metal from H and V at Fermi energy
/// e_f — the only honest "DC" output (DESIGN.md section 2).
inline double drude_weight(const CsrMatrix&, const CsrMatrix&, double) {
  throw std::logic_error("drude_weight: not implemented (DESIGN.md R4)");
}

} // namespace kubo
