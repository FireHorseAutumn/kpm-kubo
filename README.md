# kpm-kubo

Kubo–Greenwood optical and DC conductivity on the verified KPM
tight-binding core: current–current Chebyshev moments, the velocity
operator on the square lattice, DC conductivity as the ω→0 limit.

Boundary, stated once: this is a module, not a model. The Hamiltonian
stays the square-lattice tight-binding core; the observable is what is
new. And a clean metal has no finite DC conductivity — it has a Drude
weight. Any DC number this code ever prints comes with a stated
broadening or disorder scale, or it does not print.

## Relationship to kpm-tightbinding

The planned end state is a module of
[kpm-tightbinding](https://github.com/FireHorseAutumn/kpm-tightbinding),
turning it from a DOS tool into a spectral and transport tool. That repo
is public and untouched by this work: nothing is copied, vendored, or
modified there, and the merge — if the gates earn it — is rung R5,
landing in that repo on its own gates. This repo is the local proving
ground.

## Status

R3 done (2026-11-05). R1 implemented `velocity_x` — the velocity operator
stored as the real antisymmetric W = [H, x], v_x = iW — and
`kubo_lehmann_ed`, the dense-ED Lehmann oracle for σ(ω) at L ≤ 8, with
gates f-sum relative error 2.5e-15 (bound 1e-12), closure of the <v_x²>
diagonal 3.6e-15 (bound 1e-12), σ(ω) ≥ 0 exactly, an exact L=2 anchor
(single transition at ω = 2, weight π/8), and the exact-DOS cross-gate
against the kpm-tightbinding reference CSVs (linf 0, L = 8..14). R2 added
`current_moments`, the current–current Chebyshev moment table
μ_nm = (1/N) Tr[v_x T_n(H′) v_x T_m(H′)] (Weisse et al., RMP 78, 275,
Eq. 155), deterministic exact trace at these sizes; gated at L = 6,
n_cheb = 512 against the ED-derived table at matched rescaling, linf
1.725e-13 (bound 1e-12). R3 adds `optical_conductivity`, σ₁(ω) from the
table with Jackson damping (Weisse Eqs. 154 + 157 at T = 0); gated at
L = 8, n_cheb = 512: f-sum on the reconstruction 2.25e-3 (bound 1e-2),
binned L1 vs the Lehmann oracle 0.0775 (regression bound 0.16, pinned
from the bring-up measurement at ~2× — derivation in the test_r3
header). R4 adds `drude_weight` — the kinetic-energy Drude weight
D = (π/2N)(−⟨K_x⟩) from a stochastic K_x-weighted DOS trace (O(n_cheb)
SpMVs per vector, not the 2D-table price) — and the clean-metal figure
at L = 128, quarter filling, n_cheb = 128, stochastic trace R = 32
(mandatory at this L). The Drude gate is green: weight within the stated
resolution Δ = 4πa/n_cheb = 0.397 of ω = 0 matches D at 0.0216 against
the pinned 0.045. The leakage gate fails against its pre-committed
bound: regular-part weight above Δ is 0.0232 against the kernel-leakage
bound 0.0026 — the
excess is genuine finite-L physics of the open lattice (boundary-assisted
absorption plus the near-Fermi transition clusters; exact value 0.0179,
decaying ~1/L), not noise, not a reconstruction bug. Per kill criterion
1 the DC output is the gated Drude weight; a finite DC number exists
only as σ(ω) at ω ≥ Δ, stated. Figure: `./build/sigma_clean_metal &&
python3 tools/plot_sigma_clean_metal.py` regenerates
data/sigma_clean_metal_L128.{csv,png} deterministically.
kpm-tightbinding is a CMake build dependency pinned at commit
c596c86; see DESIGN.md §6 for the decision. Build and test: `cmake -S . -B build && cmake --build
build && ctest --test-dir build` (host-only, needs Eigen3; clones the
pinned kpm-tightbinding commit from GitHub into the build directory —
point CMake's `KPM_TIGHTBINDING_REPO` at a local checkout to use one
instead).
