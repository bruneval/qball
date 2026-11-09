# Traveling pseudopotentials (ETF) for a projectile at constant velocity — norm-conserving KB

Plan for adding the velocity-dependent nonlocal pseudopotential of
M. Stengel, M. Royo, E. Artacho, *Phys. Rev. Lett.* **136**, 196401 (2026)
(and its Supplemental Material, "SM") to QBALL's RT-TDDFT / Ehrenfest code.

Scope: projectile(s) moving at **constant velocity**, **norm-conserving
Kleinman-Bylander** pseudopotentials with diagonal coefficients (the only
nonlocal form QBALL accepts). Hartree atomic units, m_e = ħ = 1.

## 1. Setting

- Each atom κ has a fixed velocity **v**_κ: the projectile has **v**_κ = **v**,
  the target atoms have **v**_κ = 0.
- Positions: **R**_κ(t) = **R**_κ(0) + **v**_κ t.
- QBALL already supports this with `set atoms_dyn IMPULSIVE`
  (`MDIonicStepper`): atoms keep the velocity given in the `atom` command,
  and forces do not move them.

## 2. Traveling nonlocal operator (main text Eq. 5; SM Eqs. 2–9)

    V^nl(r, r') = Σ_κ exp(i v_κ·(r − r')) v_κ^nl(r − R_κ, r' − R_κ)

Local operators are unchanged, because exp(i v·r) commutes with any local
operator: local pseudopotential, Hartree, exchange-correlation, core charge
(NLCC) and Ewald.

## 3. Plane-wave form: a shift of the wavevector, one per atom (main text footnote 36)

Today QBALL builds the projectors of species s at q = k + G − A, where A is
the laser vector potential (A = 0 without a laser). This gives:

- the table `twnl[is][ipr][G]` in `NonLocalPotential::update_twnl`, built
  from |q|, 1/|q| and q_x, q_y, q_z;
- the phase exp(−i q·R) in `NonLocalPotential::energy` (`kpgr`, `ckpgr`,
  `skpgr`).

With traveling pseudopotentials, each atom κ uses its own vector:

    q_κ(G) = k + G − A − v_κ

- Projector: β_κp(G) = (−i)^l t_p(q_κ) exp(−i q_κ·R_κ), where t_p is
  exactly what `update_twnl` computes today, evaluated at q_κ.
- Projection: P^n_κp = Σ_G β*_κp(G) c_n(G)
- Nonlocal energy: E_nl = Σ_n f_n Σ_κp D_p |P^n_κp|², with
  D_p = `wt[is][ipr]`
- Action on a wavefunction:
  (V^nl ψ_n)(G) = Σ_κp D_p β_κp(G) P^n_κp

Why this is cheap here:

1. Using q_κ instead of k+G in the phase only multiplies β_κp by the constant
   exp(i v_κ·R_κ), which cancels between |β⟩ and ⟨β|. So the shift can go
   everywhere, as the vector-potential code already does.
2. With v_κ constant (and no laser), the projector table of a moving atom
   never changes. It is computed once, and only its phase is updated each
   step, as now. The target atoms (v = 0) keep the shared per-species table.
   Only moving atoms need a table of their own.
3. With v_κ = 0 for every atom, everything reduces exactly to the current
   code.

## 4. Conserved energy and stopping power (main text Eqs. 6–7; SM Sec. IV.A)

Velocity derivative (the "gauge field"):

    A^nl_κ = −⟨∂V^nl/∂v_κ⟩ = i⟨[V^nl_κ, r]⟩

Plane-wave form: v_κ enters only through q_κ, so ∂/∂v_κ = −∂/∂q. The
phase's contribution cancels in |P|². This leaves:

    A^nl_κα = Σ_n f_n Σ_p D_p 2 Re[ P^n*_κp ∂P^n_κp/∂q_α ]

where ∂P/∂q_α differentiates only t_p(q), not the phase.

Energy of the whole system:

    K = Σ_κ [ M_κ v_κ²/2 + A^nl_κ·v_κ ] + E_KS

Energy balance: at constant velocity, dK/dt is the power the driving agent
puts in. M v²/2 is constant, so:

    ΔE(t) = [E_KS + A^nl_proj·v](t) − [E_KS + A^nl_proj·v](0)
    S     = (1/|v|) d/dt [E_KS + A^nl_proj·v]        (stopping power)

In practice:

- E_KS must be computed with the traveling projectors.
- The A^nl·v term oscillates without growing, so it only shifts the
  instantaneous S. It does not change the average over many lattice
  periods.
- First implementation: compute A^nl by a central finite difference of
  E_nl,κ with respect to v_κ. This involves only the moving atom and is
  diagnostic only. Later, an analytic version can reuse the derivative code
  the stress calculation already has (`dvnlg`, `dtwnl`).

## 5. Forces

- The Ehrenfest equation becomes d/dt(M Ṙ + A^nl) = −∂E/∂R.
- With `IMPULSIVE`, forces do not move the atoms, so no change is needed for
  the dynamics.
- The force actually exerted on the projectile (the drag) is
  F = −∂E/∂R − dA^nl/dt. If QBALL reports or uses forces in this mode, they
  need the dA^nl/dt term.

## 6. Requirements

- Complex wavefunctions: the traveling operator is not real (it breaks
  time-reversal symmetry). At the Γ point this needs
  `set force_complex_wf ON`, which RT-TDDFT runs at Γ already need.
- Kleinman-Bylander norm-conserving pseudopotentials only (`oncv_`). Exclude
  semilocal pseudopotentials (`nquad > 0`) and ultrasoft ones.

## 7. The initial state (decided: option 1, see section 10)

At t = 0, the projectile is set moving at v.

- Option 1 (plain switch-on, the usual stopping-power setup): start from the
  ground state computed with the projectile at rest and v = 0 projectors,
  then turn on v.
- Option 2 (boosted start): for the projectile's own bound electrons, the
  Galilean-consistent start is the boosted state exp(i v·r) ψ₀ near the
  projectile (SM Eq. 2), and its energy is E₀ + N_e v²/2 (SM Eq. 10).

Option 1 is simpler and is what is usually done. Option 2 removes the
transient caused by the sudden start.

## 8. Validation tests

1. No velocity: with v = 0 everywhere, results are bit-for-bit identical to
   today's.
2. Galilean test: an isolated atom moving through an empty box. With the
   traveling projectors, the energy must stay at E₀ + N_e v²/2. This is
   implemented as the regression test of section 9, which fails with the
   current code.
3. Finite-difference check: the analytic A^nl matches the finite-difference
   A^nl.

## 9. Regression test: moving Ne atom (`testsuite/tddft/02_neon_moving_etf`)

Isolated Ne atom moving at constant velocity in an empty box. Its expected
energies are exact, not taken from a previous run. It fails with the current
code and must pass once the traveling projectors are implemented.

Pseudopotential: `Ne.psp8`, PseudoDojo ONCVPSP PBE v0.4 standard,
downloaded from
https://github.com/PseudoDojo/ONCVPSP-PBE-PDv0.4/tree/master/Ne
(8 valence electrons, 2 s + 2 p projectors, no core charge).

How the boost is done, with no new code: a constant uniform vector
potential A = −u is a pure gauge transformation that gives every electron
the velocity u. QBALL already applies A to both the kinetic energy,
|k+G−A|²/2, and the projectors, q = k+G−A. With the nucleus also moving at u
(`atoms_dyn IMPULSIVE`):

- Exact solution: the at-rest ground state ψ₀ translated with the nucleus,
  ψ(r,t) = ψ₀(r − ut), up to a global phase. The total energy stays exactly
  at E₀ + N_e u²/2.
- With traveling projectors, the moving Ne uses q = k+G−A−u = G, i.e. the
  at-rest projectors, consistent with that exact solution.
- The current code uses q = G + u, which is wrong.

Input (`neon_moving.inp`):

1. Ground state with the atom at rest: 14 bohr cubic box, PBE, 40 Ry,
   `force_complex_wf ON`, converged to 1e-10.
2. `set vector_potential -0.5 0 0`, `atoms_dyn IMPULSIVE` with the atom
   velocity u = 0.5 a.u. along x (given in the `atom` command), ETRS with
   `TD_dt` = `dt` = 0.025 a.u., 20 steps.

Checks (`neon_moving.test`):

| Check | Expected (exact) | Current code |
|---|---|---|
| ground-state energy E₀ | −33.28838802 | −33.28838802 (pass) |
| energy after TD step 1 (tol. 1e-3) | E₀ + 1 = −32.28838802 | −32.10467497 (fail, 0.184 Ha off) |
| energy after TD step 20 (tol. 1e-3) | E₀ + 1 = −32.28838802 | −32.32353897 (fail, 0.035 Ha off) |

Why the setup can be trusted:

- Same run with the nonlocal part of Ne removed (local potential only, for
  which Galilean invariance already holds): the current code gives exactly
  E₀ + N_e u²/2 at the first step (+1.00000000 Ha for u = 0.5, +4.00000000
  for u = 1.0). So the error above comes only from the nonlocal projectors.
- The remaining drift in that local-only run is a time-step error: over
  t = 0.5 a.u. it is 3.9e-4, 9.6e-5 and 2.4e-5 Ha for dt = 0.05, 0.025 and
  0.0125 (scales as dt²). The nonlocal error stays at 0.184 Ha for all three
  time steps, so it is not a time-step artefact.
- Control with the atom at rest and A = 0: the energy stays at E₀ to 1e-7 Ha.

Expected after the fix: the two TD checks within about 1e-4 Ha of
E₀ + 1 Ha (the time-step error of the local-only run).

## 10. Decisions (agreed before implementation)

1. **Input option, default off.** A new variable
   `set traveling_pseudo ON|OFF`, default `OFF`. With `OFF`, results are
   identical to the current code, including all existing tests (e.g. Al₃₂
   MD). The moving Ne test (section 9) gets `set traveling_pseudo ON` once
   the option exists. (Adding it before then would make QBALL stop with
   "no such variable" instead of showing the physics error.)
2. **Scope.** Only Ehrenfest / RT-TDDFT propagation, not ground-state or
   Born-Oppenheimer runs, and only with `atoms_dyn IMPULSIVE` (constant
   velocity). With the option on and any other ionic dynamics (e.g. `MD`),
   QBALL stops with a clear error, because the forces would need the
   missing dA^nl/dt term.
3. **Initial state: plain switch-on** (section 7, option 1). The ground state
   is computed with the atom at rest and v = 0 projectors; the velocity
   starts at t = 0. No extra code.
4. **Energy output.**
   - `<etotal>` stays the Kohn-Sham energy E_KS, so existing scripts and
     tests read the same quantity.
   - A^nl and A^nl·v of each moving atom are printed on separate lines, so
     the deposited energy E_KS + A^nl·v (section 4) can be built from the
     output.
   - A^nl is first computed by central finite difference in v; an analytic
     version may follow.
5. **Forces.** The force code is left unchanged. In `IMPULSIVE` the forces
   do not move the atoms; the reported force on a moving atom lacks the
   dA^nl/dt term (section 5), and this is documented rather than fixed.
6. **Unsupported cases stop with an error** when the option is on:
   - a moving atom with an ultrasoft or semilocal (`nquad > 0`)
     pseudopotential;
   - real wavefunctions (Γ point without `set force_complex_wf ON`).
7. **Tests.**
   - The moving Ne test (section 9) with `set traveling_pseudo ON` must pass.
   - With the option on and all velocities zero, results must be
     bit-for-bit identical to the option off.
   - When the analytic A^nl exists: it must match the finite-difference
     A^nl.

## 11. Implementation status

Done (branch `etf`):

- Input variable `traveling_pseudo` (`src/vars/TravelingPseudo.h`, flag
  `ctrl.traveling_pseudo`), registered in `src/main/qb.cc` and
  `src/qball/qbLink.cc`.
- `NonLocalPotential`:
  - the old body of `update_twnl` is now `compute_twnl(shift, only_is, ...)`,
    which builds the projector tables at q = k+G−A−shift. With shift = 0 it
    runs exactly the old code.
  - `update_twnl_atoms()` builds one table per moving atom
    (`twnl_atom_[is][ia]`, at shift = v); it stops with an error for
    ultrasoft or semilocal pseudopotentials.
  - `energy()` uses the moving atom's own table where the projectors are
    assembled (energy and Hψ) and in the force loop. The phase and the
    force factor k+G−A are unchanged (the constant phase exp(i v·R)
    cancels).
- `EnergyFunctional` switches it on only for RT-TDDFT with
  `atoms_dyn IMPULSIVE`; it stops with an error for other ionic dynamics
  (except `LOCKED`, where it has no effect), for real wavefunctions, and for
  `stress ON` (the stress code does not use the per-atom tables).

Checks:

- Moving Ne test (section 9), now with `set traveling_pseudo ON` and
  `TD_dt` = `dt` = 0.0125 a.u., 40 steps: the first step gives exactly
  E₀ + 1 Ha (−32.28838802); after 40 steps the error is 1.9e-4 Ha.
  With `OFF` the error stays 0.184 Ha.
- The remaining error is a time-step error: at t = 0.5 a.u. it is 3.2e-3,
  7.7e-4, 1.9e-4 and 5.0e-5 Ha for dt = 0.05, 0.025, 0.0125 and 0.00625
  (scales as dt²).
- Option ON with all velocities zero: output identical to OFF (bit for bit).
  Option ON with `atoms_dyn LOCKED`: identical to OFF.
- Errors checked: `atoms_dyn MD`, missing `force_complex_wf ON`.
- `make check-short`: 17/17 pass; the other 16 tests give identical results.

A^nl printing (decision 4), also done:

- Analytic version (used for the output):
  `NonLocalPotential::anl_velocity_analytic(sd, is, ia)` computes
  A_α = Σ_n occ_n Σ_p wt_p/Ω · 2 Re[P*_np ∂P_np/∂q_α] for the moving atom
  only, with the same projection conventions and MPI reductions as the
  norm-conserving part of `energy()`. The q-gradient of each projector table
  comes from `compute_twnl_gradient`: with t_p(q) = v(|q|) Y_lm(n),
  n = q/|q|,
  ∂t/∂q_α = v′ n_α Y + (v/|q|)(∂Y/∂n_α − n_α Σ_β n_β ∂Y/∂n_β),
  where v, v′ come from `Species::dvnlg` and Y_lm, ∂Y/∂n are the same 16
  real harmonics (l = 0–3, same constants and m order) as in
  `compute_twnl` (file-local helper `ylm_gradient`).
- Finite-difference version (kept as a reference):
  `NonLocalPotential::anl_velocity(sd, is, ia)`, A^nl = −∂E_nl/∂v by central
  finite difference (δ = 1e-4 a.u.): for each direction the moving atom's
  table is rebuilt at v ± δ, `energy()` is evaluated with the wavefunctions
  fixed, and the table is restored.
- `EnergyFunctional::anl_velocity(psi, is, ia, analytic = true)`: sum over
  spins and k-points with the same weights as `enl_`.
- `EhrenSampleStepper` prints, after `<etotal>`, for each moving atom:
  `<anl_velocity> name Ax Ay Az </anl_velocity>` and
  `<anl_dot_v> name A·v </anl_dot_v>` (10 decimals).

Checks of A^nl:

- Finite difference vs slope of the first-step energy in v:
  [E(0.51) − E(0.49)]/0.02 = 0.720429 vs −A^nl(0.5) = 0.720440.
- A^nl = 0 (to 1e-4) for the boosted Ne of section 9.
- Analytic vs finite difference (decision 7, third test):
  - moving Ne (s, p projectors), 42 steps: max difference 1.0e-9;
  - moving Au (PseudoDojo `Au.upf`: s, p, d, f projectors), velocity
    (0.3, −0.2, 0.4), vector potential (0.1, 0.05, −0.02): max difference
    1.4e-9 (|A| up to 0.37);
  - sensitivity: multiplying only the l = 3 gradients by 1.01 raises the Au
    difference to 1.2e-3, so the f-projector formulas are really tested.
- Cost on the Ne test (42 steps): 0.501 s/step analytic, 0.502 s/step
  finite difference, 0.497 s/step without A^nl. For a single light atom the
  propagation dominates; the analytic version should matter in large cells,
  where each finite-difference evaluation goes over all atoms (6 per moving
  atom per step). Not measured on a large system.
- `make check-short`: 18/18; all values identical to the finite-difference
  version.

## 12. Regression test: Ne sudden start, A^nl (`testsuite/tddft/03_neon_sudden_anl`)

Same Ne atom as section 9, placed at the centre of the 14 bohr box, but
**without** the vector-potential trick: the ground state is computed at rest,
then the nucleus is suddenly set moving at v = 0.5 a.u. along x
(`IMPULSIVE`, `traveling_pseudo ON`), while the electrons start at rest. Seen
from the atom, the electrons get a kick and slosh around the nucleus.
80 Ry, ETRS, `TD_dt` = `dt` = 0.025 a.u., 42 steps (about 30 s).

Physics: with velocity-dependent projectors, the work done by the agent that
keeps v constant is

    Δ(E_KS + A^nl·v) = v·ΔP,     P = d⟨Σ r⟩/dt = ⟨p⟩ + A^nl

(P = physical electron momentum, ⟨p⟩ = canonical momentum). Also
ΔE_KS = v·Δ⟨p⟩, so E_KS alone misses v·ΔA^nl: a spurious oscillation.

How the quantities are obtained (all printed together for step N, for the
same wavefunctions):

- ⟨p⟩ = −`total_electronic_current` (electron charge −1; canonical part
  only, no nonlocal contribution).
- P = −d(dipole)/dt, central difference of `<electronic_dipole_density>`
  (section 13).
- A^nl: `<anl_velocity>` (section 11).

Findings:

- E_KS jumps at the first step (+0.1837 Ha at 40 Ry): the nonlocal energy
  changes when the projectors switch from at-rest to traveling. Consistent
  with A^nl: the jump is about −½ v A^nl(t=0).
- Afterwards ΔE_KS = v·Δ⟨p⟩ holds within 2.4e-4 Ha per step.
- **The cutoff must be converged.** The relation P = ⟨p⟩ + A^nl assumes a
  complete basis: in a truncated plane-wave basis the Hamiltonian (also its
  local part) does not commute exactly with r. Measured A^nl − (P − ⟨p⟩)
  over the first atomic unit: −0.06 to −0.10 at 40 Ry (independent of the
  time step), +0.003 to +0.008 at 80 Ry and the same at 120 Ry. At 40 Ry
  A^nl itself is unconverged (0.19 instead of 0.27 at t = 0.25).
- Box size matters too: in a 12 bohr box the mismatch grows to 0.024
  (some density ejected by the kick reaches the cell boundaries sooner).
- At 80 Ry, v·A^nl goes from 0.113 Ha (step 10) to 0.067 Ha (step 40).

Why not the Wannier dipole (`wf_diag MLWFC`): it was tried first. It has
isolated bad values (one step off by about 0.1), it is computed after the
propagation (one step later than the energy and current, which biased the
first estimates of A^nl), and it is 5 to 10 times less accurate (P = 4.004
instead of 4 on the boosted Ne of section 9).

Checks (`neon_sudden.test`, computed by `anl_check.sh` from the output):

| Check | Expected | Tolerance | Measured |
|---|---|---|---|
| ground-state energy | −33.35662277 | 1e-4 | −33.35662277 |
| `ks_balance`: ΔE_KS − v·ΔP, steps 10→40 | 0.04330732 (regression value: the missing −v·ΔA^nl) | 1e-2 | 0.04330732 |
| `anl` at step 10: A_x − (P − ⟨p⟩) | 0 (exact; A^nl ≈ 0.22) | 2e-2 | 0.0083 |
| `anl` at step 40 | 0 (exact; A^nl ≈ 0.13) | 2e-2 | 0.0046 |
| `energy_balance`: Δ(E_KS + A^nl·v) − v·ΔP, steps 10→40 | 0 (exact) | 5e-3 | −0.0020 |

Before A^nl was printed, the last three checks failed (no output); with
A^nl = 0 the `anl` checks would be off by about 0.2.

## 13. Density dipole output (`set density_dipole ON`)

Added for the test of section 12, opt-in (default `OFF`):

- `src/vars/DensityDipole.h`, flag `ctrl.density_dipole`, registered in
  `src/main/qb.cc` and `src/qball/qbLink.cc`.
- `ChargeDensity::dipole()`: −∫ r n(r) d³r with r = (i/np0) a + (j/np1) b +
  (k/np2) c on the real-space grid of the cell. Meaningful only for an
  isolated system whose density is negligible at the cell boundaries.
  Spin-unpolarized only (stops with an error otherwise).
- Printed by `EhrenSampleStepper` at every RT-TDDFT step, right after
  `<etotal>`, as `<electronic_dipole_density> Dx Dy Dz </electronic_dipole_density>`
  (10 decimals).

Validation on the boosted Ne of section 9, moved to the box centre:
the dipole starts at −55.99995 (exact −56) in each direction; y and z stay
constant to 6e-9; P = −dD_x/dt is between 3.9992 and 3.9998 at every step
(exact 4, the difference is a time-step error).

## 14. Regression test: Si projectile in Si (`testsuite/tddft/04_silicon_stopping_001`)

After T. Jarrin, Phys. Rev. B 112, 245202 (2025): Si⁰ projectile in bulk Si,
⟨001⟩ channel, "off-center ¼" trajectory, 12-valence-electron
pseudopotential (his Figs. 2 and 4; at 1 MeV he finds 162 eV/Å off-center ¼
vs 75 eV/Å in the channel centre, Table III). Jarrin did not provide input
files; the setup below is ours.

- Projectile: neutral Si at 1 MeV, v = 1.19817 a.u. along z
  (`atoms_dyn IMPULSIVE`, target atoms fixed, `traveling_pseudo ON`).
- Cell: 1×1×1 cubic (8 Si atoms, a = 5.431 Å = 10.2631 bohr) instead of
  Jarrin's 3×3×3; Γ point only.
- Pseudopotential: PseudoDojo ONCVPSP LDA scalar-relativistic, no core
  correction, `Si-sp` (12 valence electrons, 2 projectors each for s, p, d,
  Kleinman-Bylander form), UPF format (UPF2, linear grid from r = 0), from
  https://raw.githubusercontent.com/PseudoDojo/ONCVPSP-LDA-SR-noCCv0.4/master/Si/Si-sp.upf
  (folder: https://github.com/PseudoDojo/ONCVPSP-LDA-SR-noCCv0.4/tree/master/Si).
  (The standard PseudoDojo LDA v0.4 set only has the 4-electron Si; Jarrin's
  Hamann-type pseudopotential is semilocal, which traveling pseudopotentials
  do not support.)
- Geometry: in projection along z the atomic columns are at (m, n)·a/4 with
  m + n even; the channel centre (a/4, 0) has 4 columns at d = a/4.
  Off-center ¼ = d/4 from the nearest column (0, 0): (x, y) = (a/16, 0).
  Start height z₀ = 7a/16, chosen so that the projectile stays at least
  3.66 bohr from all atoms along the path (Jarrin requires ≥ 0.8 Å).
- Ground state of 8 atoms + neutral projectile, LDA, 45 Ry, Fermi smearing
  (8 empty bands), then ETRS, `TD_dt` = `dt` = 0.02 a.u., 20 steps (about
  22 s). This covers only the initial transient: the test is a regression
  test and does not reproduce Jarrin's stopping power (crossing one lattice
  constant takes about 8.6 a.u., i.e. hundreds of steps).

Observations (sudden start, the projectile's own 12 electrons start at
rest):

- The energy jumps by +1.98 Ha at the first step (switch-on of the
  traveling projectors), consistent with −½ A^nl·v (A^nl·v = −3.70 Ha).
- The energy then rises by about 17 Ha over 0.4 a.u. (the bound electrons
  are dragged along: 12 v²/2 ≈ 8.6 Ha alone) and A^nl·v swings between
  −3.7 and +1.3 Ha: the traveling correction is large for a heavy
  projectile.
- Time step: dt = 0.02 vs 0.01 changes the energy after 0.38 a.u. by
  0.018 Ha (0.1 %) and A^nl·v by 0.004 Ha.
- Cutoff: not converged at 45 Ry (energy rise 16.84 / 17.48 / 17.62 Ha at
  45 / 60 / 76 Ry), acceptable for a regression test.

Checks (`silicon_stopping.test`; regression values from our run):

| Check | Value | Tolerance |
|---|---|---|
| ground-state energy | −826.34157529 | 1e-4 |
| energy after step 20 | −809.50031578 | 1e-3 |
| projectile A^nl·v after step 20 | 0.8888255793 | 1e-3 |

UPF vs PSP8: the same pseudopotential is also distributed as `Si-sp.psp8`
(both written by the same ONCVPSP run). The test was first built with the
PSP8 file; with the UPF file:

- the energy rise (step 1: 1.98335902 vs 1.98335922 Ha; step 20:
  16.84125951 vs 16.84126074 Ha) and A^nl·v (step 1: −3.67744073 vs
  −3.67744089; step 20: 0.88882558 vs 0.88882564) agree within about 1e-6 Ha;
- the total energy is shifted by a constant +8.3e-4 Ha (9e-5 Ha per atom),
  the same in the ground state and at steps 1 and 20, almost entirely in the
  local pseudopotential term (`eps` +8.3e-4, all other terms below 1e-5).
  The likely reason is the range of the tabulated local potential: 15.1 bohr
  in the UPF file, 6 bohr in the PSP8 file (QBALL continues it with −Z/r).
  Same kind of offset as found earlier between the UPF and PSP8 versions
  of the PseudoDojo C pseudopotential (3e-4 Ha for diamond).
