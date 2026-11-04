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

R1 done (2026-11-04). Implemented: `velocity_x` — the velocity operator
stored as the real antisymmetric W = [H, x], v_x = iW — and
`kubo_lehmann_ed`, the dense-ED Lehmann oracle for σ(ω) at L ≤ 8. All
pre-committed gates pass at tolerance: f-sum rule relative error
2.5e-15 (bound 1e-12), closure of the <v_x²> diagonal 3.6e-15 (bound
1e-12), σ(ω) ≥ 0 exactly, plus an exact L=2 anchor (single transition
at ω = 2, weight π/8) and the exact-DOS cross-gate against the
kpm-tightbinding reference CSVs (linf 0, L = 8..14). kpm-tightbinding is
a CMake build dependency pinned at commit c596c86; see DESIGN.md §6 for
the decision. The R2–R4 interfaces remain throwing stubs. Build and
test: `cmake -S . -B build && cmake --build build && ctest --test-dir
build` (host-only, needs Eigen3; clones the pinned kpm-tightbinding
commit from GitHub into the build directory — point CMake's
`KPM_TIGHTBINDING_REPO` at a local checkout to use one instead).
