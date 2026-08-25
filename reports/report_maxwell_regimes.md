# Report: Maxwell Visco-Elastic Two-Regime Behavior (`fig_maxwell_regimes.jl`)

## Purpose

Reproduces `fig/fig_maxwell_regimes.png` as specified in `tex/paper/paper_v2.tex` (Section "Maxwell case: two dimensionless groups and three regimes", \label{fig:maxwell_regimes}, lines 739-801), panels (a) and (b), stacked vertically to fit a single column of the two-column GJI layout. Illustrates the elastic-rebound-to-viscous-diffusion crossover of the Maxwell step response, in space (panel a) and in time (panel b).

---

## Physics

- Purely viscous velocity: `get_v_step(x,t,b,a)`; purely viscous displacement: `get_u_step(x,t,b,a)`
- Maxwell velocity: `get_v_relax(x,t,b,a,tau_R)`; Maxwell displacement: `get_u_viscoelastic(x,t,b,a,tau_R)`
- Maxwell relaxation time: τ_R := η/μ_m (paper eq. for τ_R), where μ_m is the mantle elastic shear modulus
- Maxwell length: ℓ_M = a·τ_R = EH/[2μ_m(1-ν²)] — a purely elastic quantity, independent of viscosity

---

## Parameters

| Parameter | Value |
|---|---|
| a (Natron CSS viscous best-fit) | 1.180e-04 m/s |
| b | 0.662 m |
| H | 10 km |
| E | 40 GPa |
| ν | 0.25 |
| η = EH/(2a(1-ν²)) | 1.808e+18 Pa·s |
| τ_R,ref = H/[a(1-ν)] | 3.5806 yr |
| ℓ_M = a·τ_R,ref | 13333 m |

τ_R,ref uses the shared-modulus closure already derived in the paper (Section 2.6.3, Eq. tau_R_shared): τ_R = H/[a(1-ν)]. It introduces no new elastic-modulus assumption -- only a, H, ν, already defined for this figure's own STEP/CSS curves -- and is independent of the CSS recurrence interval T, which has no role in the step-response regime this figure illustrates.

---

## Results

Panel (a) near-field (x=ℓ_M) velocity reduction at τ_R/τ_R,ref=5 relative to viscous: **61.0%** (paper caption states "more than 60%").

Panel (b): both Maxwell curves (τ_R/τ_R,ref = 1, 5) lag the purely viscous displacement for t ≲ τ_R (elastic-rebound zone, grey shading) and converge onto the common asymptote u/b = 1 (dotted line) once t ≫ τ_R (x-axis extended to t/τ_R,ref=20 to show full convergence). The orange band marks the displacement deficit of the slowest-relaxing case (τ_R/τ_R,ref=5) relative to the viscous prediction.

---

## Dependencies

- `src/models.jl` — `get_v_step`, `get_v_relax`, `get_u_step`, `get_u_viscoelastic` (no modifications needed)
- `tex/paper/paper_v2.tex` — figure specification (\label{fig:maxwell_regimes}) and physical parameter definitions
- Natron CSS best-fit `a` reproduced from `fig_behavior_natron.jl`

Generated: 2026-08-14T16:58:13.404
