# Report: Maxwell–Viscous Convergence (`fig_maxwell_convergence.jl`)

## Purpose

Demonstrates that the Maxwell visco-elastic step response converges to the purely viscous step response once t ≫ τ_R, from two complementary, fully parameter-free (self-similar) angles: (a) three observation points ξ = x/ℓ_M, each shown as a time series in universal units t̂ = t/τ_R, and (b) the spatial velocity-profile broadening this convergence produces at fixed snapshots of t̂, in the style of Becker & Faccenna's postseismic-relaxation figures. Generates `fig/fig_maxwell_convergence.png` (2 panels).

---

## Physics

- Purely viscous asymptote (dimensional form): `get_u_step(x, tau, b, a) = (b/π)·atan(a·tau/x)`
- Maxwell step response: `get_u_viscoelastic(x, tau, b, a, tau_R)`, solving
  ∂û/∂t = -a|k|/(1+τ_R·a|k|) · û,  û(k,0) = b/|k|, via a segmented frequency-domain
  quadrature (`_k_limit_maxwell` + `_segmented_quadgk` in `src/models.jl`).
- Maxwell velocity-profile response: `get_v_relax(x, tau, b, a, tau_R)`, the companion spatial-velocity solution used in panel (b); purely viscous velocity reference: `get_v_step(x, tau, b, a)`.

**Non-dimensionalization (basis for both panels):** a Maxwell material has exactly one intrinsic timescale, τ_R = η/μ, and one intrinsic lengthscale, ℓ_M = a·τ_R. Non-dimensionalizing with t̂ = t/τ_R and ξ = x/ℓ_M = x/(a·τ_R) eliminates every free parameter (a, b, τ_R all cancel), reducing the whole system to one universal function F(ξ, t̂) = u(x,t)/b — identical for any material. Panel (b) shows F at fixed t̂, varying ξ (velocity version, 5 curves — a spatial-profile slice). Panel (a) shows the displacement version of F at fixed ξ, varying t̂ (3 curves) — a genuine dual/complementary slice of the same universal function, not just a relabeling of axes.

---

## Parameters

This figure is now **fully parameter-free** — both panels use unit values a = b = τ_R = 1 (so t̂ = t and ξ = x directly), and no Natron-specific numbers are used anywhere in the script.

**Historical note on the three ξ values:** panel (a)'s three cases, ξ ∈ {2, 1, 1/3}, originate from the previous revision's fixed-observation-point framing, which used Natron CSS best-fit parameters (a = 1.18e-4 m/s, b = 0.662 m) and a reference relaxation time τ_R,ref = T/4 ≈ 5.75 yr (T = 23.0 yr, the Natron **CSS** recurrence interval, borrowed purely as a source of a round number — it has no physical role in this **step**/single-dike response). That panel fixed x₀ = a·τ_R,ref and swept τ_R/τ_R,ref ∈ {0.5, 1, 3}. In ξ terms, ξ = x₀/ℓ_M(τ_R) = x₀/(a·τ_R) = τ_R,ref/τ_R, so those three cases map exactly onto ξ ∈ {2, 1, 1/3} — which is why panel (a) keeps the same colors (`:firebrick`, `:orangered`, `:sandybrown`) as that earlier revision. No Natron numbers are computed in the current script; this is context only.

---

## Result — Panel (a): self-similar limit at three ξ

At t̂ = t/τ_R = 8.00 (t_hat_max), normalized deviation |u_Maxwell/b − u_viscous/b| from each ξ's own viscous asymptote:

| ξ | |Δ(u/b)| at t̂ = 8 |
|---|---|
| 2 | 0.0225 |
| 1 | 0.0146 |
| 0.333 | 0.0059 |

All three ξ cases converge to within 0.08 (normalized displacement) of their own viscous asymptote by t̂ = 8, confirming that Maxwell visco-elasticity reduces to purely viscous behavior once t̂ ≫ 1, independent of ξ.

**Why three separate dashed viscous curves, not one shared curve (the key subtlety of this consolidation):** in the previous (years-based) revision, all three Maxwell curves shared ONE black viscous asymptote, because `get_u_step` does not depend on τ_R — with x₀ fixed, the viscous limit is one curve in absolute time. Once the x-axis becomes t̂ = t/τ_R (each curve normalized by **its own** τ_R), that shared curve no longer applies: the viscous asymptote becomes ξ-dependent too, `u_step(x=ξ·ℓ_M, t=t̂·τ_R)/b = (1/π)·atan(t̂/ξ)`. So this panel plots three dashed reference curves, one per ξ (same color as its solid Maxwell curve), not a single shared black curve — this is a direct, visible consequence of switching the time axis from years to t̂.

**Note on convergence window:** t_hat_max was widened empirically (same pattern as prior revisions' `grid_mult`/`s_max`) until |Δ(u/b)| < 0.08 for **all three** ξ. Notably, ξ = 2 (the case with the smallest τ_R, formerly the fastest-relaxing curve in years) is actually the **slowest to converge in t̂ units** — since t̂ = t/τ_R, the same physical time corresponds to a larger t̂ for smaller τ_R. t_hat_max = 8 clears the tolerance for all three cases with a comfortable margin and mirrors the "8·τ_R" window convention already established in earlier revisions of this figure.

**Numerical note — unit scale:** `get_u_viscoelastic` in `src/models.jl` has a small-|x| fallback (`abs(x) < 1.0 && return b/2.0`), written for physical x in meters. Calling it literally at x = ξ with a = b = τ_R = 1 makes this guard fire spuriously for ξ = 1/3 (since 1/3 < 1), returning a constant b/2 at every t̂ instead of the real Maxwell curve. Since `src/models.jl` is not modified, this script instead scales a (and hence ℓ_M = a·τ_R, and hence x = ξ·ℓ_M) up by a factor of 10 while keeping τ_R = 1, so t̂ = t still holds exactly and x stays clear of the guard for all three ξ. This is a pure change of arbitrary units, not a physics change — verified numerically that results for ξ = 1, 2 (already unaffected by the guard) are identical to machine precision (< 1e-15) whether computed at unit scale or at this larger scale, confirming `get_u_viscoelastic`'s genuine self-similarity once x is safely above the guard threshold.

---

## Result — Panel (b): velocity-profile broadening

At early normalized times (t̂ = t/τ_R = 0.1, 0.3), the Maxwell velocity profile v̂(x̂) is narrow and peaked close to x̂ = 0 — an elastic, locked-fault-like signature. As t̂ increases (0.8, 1.5, 3.0), the profile progressively broadens and flattens, spreading strain rate over a wider region. By t̂ = 3.0 (firebrick), the Maxwell profile has visibly converged onto the purely-viscous reference profile (black dashed, `get_v_step` at t̂ = 3.0), confirming the same viscous-limit convergence seen in panel (a), now viewed as a diffusive spatial broadening rather than a time series at one point — the Becker & Faccenna-style complementary view requested for this figure.

---

## Dependencies

- `src/models.jl` — `get_u_viscoelastic`, `get_v_step`, `get_v_relax` (no modifications needed; `get_u_step` no longer called directly — panel (a)'s viscous asymptote uses the equivalent non-dimensionalized closed form (1/π)·atan(t̂/ξ) instead)
- Panel (b) color ramp and t̂ snapshot set reused verbatim from `fig_behavior.jl` (project convention)

Generated: 2026-07-02T13:23:32.019
