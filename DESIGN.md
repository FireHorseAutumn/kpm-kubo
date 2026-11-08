# DESIGN.md — kpm-kubo

Design record for the Kubo–Greenwood conductivity module on the verified
KPM core. Status: R3 done (2026-11-05). R1: `velocity_x` and
`kubo_lehmann_ed` are implemented and gated (measured: f-sum relative
error 2.5e-15 against the 1e-12 bound; closure 3.6e-15 against 1e-12;
sigma >= 0 exactly; L=2 exact anchor at omega = 2 with weight pi/8;
exact-DOS cross-gate against the kpm-tightbinding reference CSVs at
L = 8..14, linf exactly 0). R2: `current_moments`, the 2D moment table
mu_nm = (1/N) Tr[v_x T_n(H') v_x T_m(H')] of Weisse et al. Eq. 155, with
the deterministic exact trace over the unit basis once n_random >= N and
the stochastic random +/-1 trace below it; gated at L = 6, n_cheb = 512,
a = 4.04 against the ED-derived table at matched rescaling, measured
linf 1.725e-13 against the 1e-12 bound; mu_00 = 2(L-1)/L anchor exact to
2.2e-16. R3: `optical_conductivity`, sigma_1(omega) from the table
(Weisse Eqs. 154 + 157 at T = 0, Jackson kernel applied at
reconstruction, h_nm per Eq. 95); gated at L = 8 with the deterministic
trace and n_cheb = 512: f-sum on the reconstruction 2.25e-3 against the
pre-committed 1e-2 bound, and binned L1 vs the R1 Lehmann oracle 0.0775
against the regression bound 0.16 pinned from the bring-up measurement
at ~2x (pin derivation in the test_r3.cpp header and STATE.md). R4:
`drude_weight` (D = (pi/2N)(-<K_x>) from a stochastic K_x-weighted DOS
trace, O(n_cheb) SpMVs per vector) and the clean-metal figure at
L = 128, quarter filling, n_cheb = 128, R = 32, stated resolution
Delta = 4 pi a / n_cheb = 0.397. Drude gate green: weight within Delta
of omega = 0 matches D at 0.0216 against the pinned 0.045. Leakage gate
RED against the pre-committed kernel bound: regular-part weight above
Delta 0.0232 vs 0.0026 (9x). The excess is measured genuine finite-L
physics of the open lattice — boundary-assisted absorption plus the
near-Fermi transition clusters at omega_d ~ 4.4 d / L, exact value
0.0179, decaying ~1/L (weight above omega = 0.4: 14.1% of total at
L = 32, 7.1% at 64, 3.5% at 128) — not noise, not a reconstruction bug.
Kill criterion 1 stands as written: the DC output is the gated Drude
weight; sigma(omega) is reported at omega >= Delta. Whether the leakage
gate gets a finite-L term in its bound is a DESIGN amendment, not a
bring-up decision. Ruled by Ian 2026-11-08, Option A: Kill 1 applied
test-side — the leakage comparison is printed in full in test_r4 and
never asserted; no bound moved. Option B, the physics-informed finite-L
leakage bound, remains owed: legitimate only after the ~1/L decay model
validates as a prediction at fresh L values before any re-gating — a
bound fitted to the run it gates is a description, not a gate. This note fixes the scope, the oracles, the
rung decomposition with pre-committed gates, and the kill criteria.

Ordering: second of the three KPM extensions and the cheapest — same
moment machinery, one new operator.

## 1. Physics scope

Kubo–Greenwood conductivity of the noninteracting square-lattice
tight-binding core: the real part of the optical conductivity
sigma(omega) from the current–current correlation, evaluated with KPM
moments; DC as the omega -> 0 limit, with the clean-metal caveat of
section 2 stated every time the word "DC" appears.

The velocity operator on the lattice: v_x = i[H, x], a bond operator —
nonzero only between nearest neighbors along x, amplitude set by t and
the lattice constant. It assembles as a CSR with H's sparsity restricted
to the x bonds. No new model: the Hamiltonian stays the verified
square-lattice core.

## 2. Oracles

1. **ED Lehmann on small clusters.** Dense ED (Eigen) at L <= 8 with
   full eigenvectors; sigma(omega) in the Lehmann representation — a
   finite sum over eigenstate pairs. An exact oracle with no KPM
   anywhere near it.
2. **Analytic clean-metal limits.** The clean one-band lattice metal
   puts all its weight in the Drude delta at omega = 0; the regular part
   vanishes for omega > 0; and the f-sum rule ties the total weight to
   the kinetic-energy expectation, computable exactly from the same
   ground state. Exact forms and prefactors to verify from the section 7
   sources before any claim. The consequence is pre-committed here, not
   discovered later: "DC conductivity" of a clean metal is a Drude
   weight, not a finite number. The module reports the Drude weight
   (gated) and sigma(omega) at omega >= a stated resolution. A finite DC
   number exists only behind a stated broadening or disorder scale.
3. **The kpm-tightbinding reference data.** The committed
   `data/exact_dos_L*.csv` files gate that the H entering the moments is
   the verified one. v_x is new; H is not.

## 3. Rungs and gates

Each rung is one sitting. Gates are pre-committed; tolerances are not
adjusted after a run.

- **R1 — velocity builder + ED Lehmann oracle (host).** `velocity_x` as
  CSR; `kubo_lehmann_ed` at L <= 8. Gates: the f-sum rule — the Lehmann
  total weight matches the kinetic-energy expectation from the same
  eigensystem, relative error < 1e-12; sigma(omega) >= 0 at every grid
  point; closure — the sum over final states reconstructs the <v_x^2>
  diagonal to 1e-12. Deliverable: the verified operator and the oracle
  every later rung is gated against.
- **R2 — current–current Chebyshev moments (host).** The 2D moment table
  mu_nm (construction taken from Weisse et al., verified there first).
  Gate at L = 6: moments against the ED-derived deterministic table,
  linf < 1e-12 at matched rescaling; if the stochastic trace is used,
  the kpm-tightbinding noise-envelope discipline applies instead
  (linf <= 4 sigma, sigma = 1/sqrt(2*N*R)).
- **R3 — sigma(omega) reconstruction (host).** From the R2 table at
  moderate L. Gates: reconstruction vs the R1 Lehmann oracle at L = 8,
  binned L1 within a regression bound pinned at R2 when the noise level
  is measured (the kpm-tightbinding L1 discipline: regression bound, not
  accuracy estimate); the f-sum re-gated on the reconstruction —
  Jackson damping redistributes weight, the integral is the conserved
  thing, bound pre-committed at 1%.
- **R4 — the clean-metal figure and the Drude gate (host).** Large L,
  fixed filling. Gates: weight within the stated resolution of omega = 0
  matches the kinetic-energy Drude weight within a pre-committed bound;
  regular-part weight at omega > resolution matches the clean-metal zero
  within the leakage bound computed from the kernel at the chosen
  n_cheb. This rung's figure is the one that turns the parent repo into
  a transport tool.
- **R5 — merge decision.** The validated module merges into
  kpm-tightbinding as its transport module — a commit in that repo, on
  that repo's gates, not before. This repo stays the proving ground.
- **R6 — GPU, only if disorder-averaged sigma studies justify it.**
  Gate: GPU vs CPU moments linf <= 1e-13 at fixed seed; benchmarks keep
  raw per-rep timings plus device/clock/driver/seed/toolchain metadata.
  Not scheduled.

## 4. Kill criteria

- **Kill 1 — fake DC.** If the omega -> 0 limit cannot be stated
  honestly (Drude weight gated, resolution stated), the DC claim is
  dropped to sigma(omega) at omega >= resolution. Pre-committed, so
  nobody negotiates with a plot later.
- **Kill 2 — oracle disagreement.** ED and moments disagree beyond the
  pre-committed tolerance and no bug is found: stop. No tolerance
  tuning.
- **Kill 3 — cost.** The 2D table costs O(n_cheb^2) SpMVs against DOS's
  O(n_cheb). If the target resolution at the target L is unreachable on
  host and no physics question justifies the GPU rung, close with the
  small-L validated module as the record.

## 5. Cost notes

Cheapest of the three extensions: same moment machinery, one more CSR
operator. The only new cost is the 2D moment table — O(n_cheb^2) SpMVs
per random vector where DOS pays O(n_cheb). At the R1–R3 sizes this is
seconds. Cheaper formulations of dynamical correlations exist in the KPM
literature; verify in Weisse et al. before adopting one. Estimates,
labeled as estimates: R1 one sitting; R2 one sitting; R3 one sitting;
R4 one to two.

## 6. Relationship to kpm-tightbinding

Future module of it — the stated end state is that repo becoming a
spectral and transport tool. Today: nothing is copied, vendored, or
modified there; this repo stands alone and touches that repo only as
committed reference data (the exact_dos CSVs) at gate time. Before R1,
`kubo::CsrMatrix` was a placeholder with the same layout so the
interfaces compiled; the merge itself is R5 here.

Decision, taken at R1: the build dependency, not the fallback.
kpm-tightbinding is consumed through CMake (FetchContent) at pinned
commit c596c86a0e343be25e735a8be2433ee163ac0809, and `kubo::CsrMatrix`
is `kpm::CsrMatrix` itself. One wrinkle: that repo's CMake project
requires CUDA at configure time, which a host-only module cannot take
on, so the pinned sources are populated without add_subdirectory and the
host core (`src/lattice.cpp`, headers) is compiled directly here.
Nothing is copied into this tree — the clone lives under `build/` — and
that repo is never modified. The committed exact-DOS CSVs cross-gate
the H entering the module at every test run (linf 0 at L = 8..14).

## 7. Primary sources needed before physics claims

To fetch and verify before any number derived from them is quoted:

- Weisse, Wellein, Alvermann, Fehske, *The kernel polynomial method*,
  Rev. Mod. Phys. 78, 275 (2006), arXiv:cond-mat/0504627 — the KPM
  dynamical-correlation and transport construction: the mu_nm table and
  its reconstruction. The house convention source; re-read, do not quote
  from memory. Verified at R2 against the arXiv text: the moment table
  mu_nm = (1/D) Tr[T_n(H~) J T_m(H~) J] (their Eq. 155) of the current
  matrix-element density j(x, y) (Eq. 153), the stochastic trace over
  random vectors (Eq. 40), and the effort scaling 2RDN to RDN^2; the
  h_nm normalization (Eq. 95) and 2D reconstruction (Eq. 154) verified
  as read at R2 and used at R3. Verified at R3: the noninteracting
  form sigma^reg(omega) = (1/omega) int j(y+omega, y) (f(y) -
  f(y+omega)) dy (their Eq. 157, the Boltzmann weights of Eq. 156
  replaced by Fermi functions) — used at T = 0 as the Fermi-window
  integral, in the house pi-convention of the R1 oracle (Weisse's own
  prefactor differs; the oracle is the authority).
- Kubo (1957), or a standard many-body text (Mahan; Fetter–Walecka), for
  the Kubo–Greenwood formula and its prefactors. Verified at R1: the
  pair-sum form sigma_1(omega) = (pi/N) sum (f_n - f_m) |v_nm|^2 /
  omega_mn delta(omega - omega_mn) was derived from first-order time-
  dependent perturbation theory and checked against Weber and
  Hohenadler, arXiv:1605.02779 (their Eq. 14 is the same pair sum).
- A lattice f-sum-rule / Drude-weight reference (the kinetic-energy
  relation) — the R4 gate numbers come from here. Verified at R1 for the
  relation used in the R1 gate: integral sigma_1 d-omega =
  (pi/2N) (-<K_x>) follows from [x,[H,x]] = -H_x plus cancellation of
  occupied–occupied pairs; the prefactor matches the continuum f-sum
  rule pi n e^2 / 2m (Haule, arXiv:1106.2309) at the band bottom, and
  Weber–Hohenadler state the lattice form integral sigma d-omega =
  -(pi/2) E_kin.

Nothing above is quoted until fetched and read.

## 8. Honest limitations

Noninteracting, clean, mean-free-path-infinite metal unless a broadening
or disorder scale is stated. DC is a Drude weight, not a finite number,
in the clean limit. No vertex corrections, no interactions, no
transport-theory novelty claim: this is a verified observable on a
verified core.
