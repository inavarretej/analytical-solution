# Report: Maxwell Viscoelastic Converging to Purely Viscous (`fig_maxwell_vs_viscous.jl`)

## Purpose

Generates `fig/fig_maxwell_vs_viscous.png`, a two-panel figure showing that the Maxwell viscoelastic model converges to the purely viscous model as the Maxwell relaxation time τ_R decreases, in both the single-event Step form (panel a) and the periodic CSS form (panel b).

---

## τ_R,ref — how it is calculated (for the figure caption)

τ_R,ref is the Maxwell relaxation time of the mantle, τ_R := η/μ_m, evaluated at the Natron best-fit characteristic velocity `a`. It is **not** derived from the CSS recurrence interval T.

**Steps:**
1. Invert the characteristic-velocity relation `a = HE/(2η(1-ν²))` for the substrate viscosity:
   `η = HE/(2a(1-ν²))`
2. Divide by the mantle elastic shear modulus `μ_m` (an independent physical constant, not derived from `a`, `E`, or `H`):
   `τ_R,ref = η/μ_m`

**Numbers used:**

| Symbol | Value | Role |
|---|---|---|
| `a` | 1.18×10⁻⁴ m/s | Natron CSS viscous best-fit (`fig_behavior_natron.jl`) |
| `H` | 10 km | assumed elastic plate thickness |
| `E` | 40 GPa | assumed crustal Young's modulus |
| `ν` | 0.25 | assumed Poisson ratio |
| `μ_m` | 60 GPa | assumed mantle elastic shear modulus (`src/parameters.jl`) |
| `η = HE/(2a(1-ν²))` | 1.808×10¹⁸ Pa·s | derived substrate viscosity |
| **`τ_R,ref = η/μ_m`** | **0.9548 yr** | derived reference Maxwell relaxation time |

This is the same τ_R,ref used in `fig_maxwell_regimes.jl` (independently cross-checked there against the paper's quoted near-field velocity reduction, 61.0% vs. the paper's ">60%").

The four Maxwell curves in both panels use `τ_R = {50, 10, 2, 0.4} × τ_R,ref` — a deliberately wide sweep (chosen after an initial narrower sweep, `{5,1,0.2,0.04}×τ_R,ref`, showed too little visual separation between curves) so that both the early divergence from, and the eventual re-convergence to, the purely viscous curve are visible within the plotted range.

---

## Physics

### Panel (a) — Step model

Single-dike-event functions, already validated in earlier figures this session (`fig_maxwell_regimes.jl`):
- Purely viscous: `get_u_step(x,t,b,a) = (b/π)·atan(at/x)`
- Maxwell: `get_u_viscoelastic(x,t,b,a,τ_R)`, solving `∂û/∂t = -a|k|/(1+τ_R a|k|)·û`

### Panel (b) — CSS model

**This panel required a new function.** The pre-existing `get_v_css_maxwell` in `src/models.jl` sums a *finite* real-space history of individual Maxwell step responses (`Σ_{j=0}^{N} get_v_relax(x, τ-jT, ...)`, `N=floor(τ/T)`) — an approximation intended for systems with few recorded events. It was verified **not** to converge to the closed-form viscous CSS (`get_v_css`) even as `τ_R → 0` and even at N=25 cycles (a diagnostic in this session found the discrepancy unchanged between N=5 and N=25, showing the finite sum was already saturated in N — the mismatch is structural, not a truncation issue).

The paper (`tex/paper/paper_v2.tex`, "Maxwell Cyclo-Steady State" section, eqs. eq:maxwell_css_v/eq:maxwell_css_u) instead derives the Maxwell CSS by imposing the Fourier-space periodicity condition `û(k,τ+T) = û(k,τ) + ib/k` directly and solving for the unique periodic steady state — a fundamentally different (and correct) construction from summing discrete past events. This is now implemented as new functions in `src/models.jl` (added without modifying any existing function):

```
get_v_css_maxwell_closed(x, τ, b, v_plate, a, τ_R)
get_u_css_maxwell_closed(x, τ, b, v_plate, a, τ_R)
```

implementing, with `α(k) := ak/(1+τ_R·ak)`:

```
v_CSS^MAX(x,τ) = (2ba/π) ∫₀^∞ [sin(kx)/(1+τ_R·ak)]·[e^(-α(k)τ)/(1-e^(-α(k)T))] dk
u_CSS^MAX(x,τ) = (2b/π)  ∫₀^∞ [sin(kx)/k]·[(1-e^(-α(k)τ))/(1-e^(-α(k)T))] dk
```

**Verified:** ratio to `get_v_css`/`get_u_css` → 1.0000 exactly as `τ_R → 0`, at every tested phase — confirming this closed form (unlike the finite-sum `get_v_css_maxwell`) is the correct object for this comparison.

`τ` is the phase since the most recent intrusion (`τ ∈ [0,T)`); panel (b) uses the fixed phase `τ = 0.5T`.

---

## Panel descriptions

### Panel (a) — Step displacement at x = aT

Log-scaled time axis (0.1–300 yr) so both the early divergence and the late re-convergence are visible in one panel — a linear axis could not show both (the widest `τ_R` case needs ~200–300 yr to rejoin the viscous curve, which would compress the early-time detail to an unreadable sliver). All four Maxwell curves visibly lag the purely viscous asymptote at early-to-intermediate times, with the lag proportional to `τ_R`, and rejoin it by t ≈ 200–300 yr.

### Panel (b) — CSS velocity at fixed phase τ = 0.5T

x-axis 0–100 km (zoomed to the near-field transition region — beyond ~75 km all curves are already at the far-field plateau). The purely viscous CSS curve rises from 0 at the dike and saturates at the plate-velocity plateau `v·T/b = 1`. The `τ_R = 50·τ_R,ref` case is nearly flat at the plateau value across the whole range (deep elastic regime: near-field stiffening suppresses the sharp viscous rise, leaving a nearly uniform velocity profile instead). Decreasing `τ_R` progressively recovers the viscous curve's characteristic rise from the dike.

Peak amplitudes (v·T/b) at fixed phase τ=0.5T: purely viscous 0.999; τ_R/τ_R,ref = 50 → 0.997; 10 → 0.990; 2 → 0.993; 0.4 → 0.997 (all ≈1, i.e. all curves reach the same far-field plateau — the visible differences are entirely in the near-field rise, not the asymptote).

---

## Dependencies

- `src/models.jl` — `get_u_step`, `get_u_viscoelastic` (existing, unmodified); `get_v_css`, `get_v_css_maxwell_closed` (new, added this session, does not replace or modify `get_v_css_maxwell`)
- Natron best-fit `a` reproduced from `fig_behavior_natron.jl`; `H`, `E`, `μ_m` shared with `fig_maxwell_regimes.jl` and `src/parameters.jl`
