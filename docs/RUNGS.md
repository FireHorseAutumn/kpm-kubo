# Rung record

The full record for kpm-kubo: every rung, every gate, every pinned
bound. The README carries the headline; this file carries the evidence.
Host-only suite (Eigen3); kpm-tightbinding is a CMake build dependency
pinned at `c596c86` (populate-only, that repo byte-untouched; the
decision is DESIGN.md §6).

## R1: velocity operator and the Lehmann oracle

- `velocity_x`: the velocity operator stored as the real antisymmetric
  W = [H, x], v_x = iW.
- `kubo_lehmann_ed`: the dense-ED Lehmann oracle for σ(ω) at L ≤ 8.
- Gates: f-sum relative error 2.5e-15 (bound 1e-12); closure of the
  <v_x²> diagonal 3.6e-15 (bound 1e-12); σ(ω) ≥ 0 exactly; the exact
  L=2 anchor (a single transition at ω = 2, weight π/8); the
  exact-DOS cross-gate against the kpm-tightbinding reference CSVs,
  linf 0 at L = 8..14.

## R2: current–current moments

- `current_moments`: the current–current Chebyshev moment table
  μ_nm = (1/N) Tr[v_x T_n(H′) v_x T_m(H′)] (Weisse et al., RMP 78, 275,
  Eq. 155); deterministic exact trace at these sizes.
- Gate: L = 6, n_cheb = 512 against the ED-derived table at matched
  rescaling, linf 1.725e-13 (bound 1e-12).

## R3: the σ₁(ω) reconstruction

- `optical_conductivity`: σ₁(ω) from the table with Jackson damping
  (Weisse Eqs. 154 + 157 at T = 0).
- Gates at L = 8, n_cheb = 512: f-sum on the reconstruction 2.25e-3
  (bound 1e-2); binned L1 vs the Lehmann oracle 0.0775 (regression
  bound 0.16, pinned from the bring-up measurement at ~2×; derivation
  in the test_r3 header).

## R4: Drude weight and the clean-metal figure

- `drude_weight`: the kinetic-energy Drude weight D = (π/2N)(−⟨K_x⟩)
  from a stochastic K_x-weighted DOS trace. O(n_cheb) SpMVs per vector,
  not the 2D-table price.
- The clean-metal figure at L = 128, quarter filling, n_cheb = 128,
  stochastic trace R = 32 (mandatory at this L).
- The Drude gate is green: weight within the stated resolution
  Δ = 4πa/n_cheb = 0.397 of ω = 0 matches D at 0.0216 against the
  pinned 0.045.
- Figure: `./build/sigma_clean_metal && python3 tools/plot_sigma_clean_metal.py`
  regenerates `data/sigma_clean_metal_L128.{csv,png}` deterministically.

## The leakage gate, demoted

- Ruled under kill criterion 1, Option A: the leakage assertion is an
  info line, printed in full with its decomposition, never asserted;
  the bound is unmoved.
- Regular-part weight above Δ measures 0.0232 against the
  kernel-leakage bound 0.0026.
- The excess is genuine finite-L physics of the open lattice:
  boundary-assisted absorption plus the near-Fermi transition clusters;
  exact value 0.0179, decaying ~1/L. Not noise, not a reconstruction
  bug. The accounting is exact for the separable spectrum: the
  decomposition prints as true finite-L 0.0179 + kernel 0.0026 +
  stochastic ~0.0007.
- Per kill criterion 1 the DC output is the gated Drude weight; a
  finite DC number exists only as σ(ω) at ω ≥ Δ, stated.

## Build and test

```sh
cmake -S . -B build && cmake --build build && ctest --test-dir build
```

Host-only, needs Eigen3. The build clones the pinned kpm-tightbinding
commit from GitHub into the build directory; point CMake's
`KPM_TIGHTBINDING_REPO` at a local checkout to use one instead.
