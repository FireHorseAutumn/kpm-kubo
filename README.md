# kpm-kubo

Kubo–Greenwood optical and DC conductivity on the verified KPM
tight-binding core: current–current Chebyshev moments, the velocity
operator on the square lattice, DC conductivity as the ω→0 limit.

Boundary, stated once: this is a module, not a model. The Hamiltonian
stays the square-lattice tight-binding core; the observable is what is
new. And a clean metal has no finite DC conductivity: it has a Drude
weight. Any DC number this code ever prints comes with a stated
broadening or disorder scale, or it does not print.

## Relationship to kpm-tightbinding

The planned end state is a module of
[kpm-tightbinding](https://github.com/FireHorseAutumn/kpm-tightbinding),
turning it from a DOS tool into a spectral and transport tool. That repo
is public and untouched by this work: nothing is copied, vendored, or
modified there, and the merge, if the gates earn it, is rung R5, landing
in that repo on its own gates. This repo is the local proving ground.

## Status

- Rungs R1–R4 done: the velocity operator and dense-ED Lehmann oracle
  (R1), the current–current Chebyshev moment table (R2), the σ₁(ω)
  reconstruction with Jackson damping (R3), and the kinetic-energy
  Drude weight with the clean-metal figure (R4).
- CTest 6/6, host-only. kpm-tightbinding is a CMake build dependency
  pinned at `c596c86`; see DESIGN.md §6.

## Gates of record

| Gate | Measured | Bound |
|---|--:|--:|
| f-sum relative error (R1) | 2.5e-15 | 1e-12 |
| <v_x²> diagonal closure (R1) | 3.6e-15 | 1e-12 |
| Moment-table linf vs ED, L = 6 (R2) | 1.725e-13 | 1e-12 |
| f-sum on the σ(ω) reconstruction (R3) | 2.25e-3 | 1e-2 |
| Binned L1 vs the Lehmann oracle (R3) | 0.0775 | 0.16 |
| Drude weight within Δ = 0.397 of ω = 0 (R4) | 0.0216 | 0.045 |

R1 also gates σ(ω) ≥ 0 exactly and an exact L=2 anchor (a single
transition at ω = 2, weight π/8); the exact-DOS cross-gate against
kpm-tightbinding's reference CSVs holds at linf 0 for L = 8..14.

## Honest losses

- The R4 leakage gate fails its original letter: regular-part weight
  above Δ measures 0.0232 against the kernel-leakage bound 0.0026. The
  excess is genuine finite-L physics of the open lattice (exact
  separable-spectrum accounting puts it at 0.0179, decaying ~1/L), not
  noise, not a reconstruction bug. Ruled under kill criterion 1, Option
  A: the gated DC output is the Drude weight, the leakage assertion is
  demoted to an info line printed in full, and the bound is unmoved.
- The R3 L1 bound is a regression bound pinned from the bring-up
  measurement at ~2×, not a derived one; the derivation is in the
  test_r3 header.

## Build and test

```sh
cmake -S . -B build && cmake --build build && ctest --test-dir build
```

Host-only; needs Eigen3. The build clones the pinned kpm-tightbinding
commit from GitHub into the build directory; point CMake's
`KPM_TIGHTBINDING_REPO` at a local checkout to use one instead.

## Layout

- `include/kubo/kubo.hpp`, `src/kubo.cpp`: the velocity operator, the
  moment table, the σ(ω) reconstruction, the Drude weight.
- `src/sigma_clean_metal.cpp`: the clean-metal figure driver.
- `tests/`: `test_r1`–`test_r4`, `test_separable` (the
  separable-spectrum accounting suite), `test_stub`.
- `tools/plot_sigma_clean_metal.py`: the figure plot.

## Deeper reading

- `docs/RUNGS.md`: the full record: every rung, every gate, the
  leakage-gate disposition.
- `DESIGN.md`: the rung decomposition and the pre-committed kill
  criteria.
