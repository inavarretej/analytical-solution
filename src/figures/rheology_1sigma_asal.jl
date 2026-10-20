using CairoMakie, LaTeXStrings, Statistics

mkpath("fig")

# ── Constants ─────────────────────────────────────────────────────────────────
ν        = 0.25
E_ref    = 40e9
E_min    = 20e9
E_max    = 80e9
H_ref_km = 10.0

# ── Best-fit a and 1σ bounds ─────────────────────────────────────────────────
# From src/figures/asal_css_v_constrained.jl / asal_chi2_landscape_v_constrained.jl:
# EV RSS, joint time series + profile, T = b/V (V = 13 mm/yr); 1σ = Δχ² ≤ 2.30
# projected onto a. (Previously a* = 9.503e-6, 1σ 9.13–9.89e-6, from the older
# fixed-T = 500 yr fit.)
a_opt   = 1.0006e-5  # m/s  — best fit
a_1σ_lo = 9.66e-6    # m/s  — lower 1σ edge → largest η
a_1σ_hi = 1.04e-5    # m/s  — upper 1σ edge → smallest η

# ── Iso-a/E curves: H(η) = 2aη(1−ν²)/E ──────────────────────────────────────
η_log  = collect(range(17.0, 21.0, length = 800))
η_vals = 10 .^ η_log
H_iso(η, a, E) = 2 * a * η * (1 - ν^2) / E / 1e3
clip(v) = clamp.(v, 0.0, 30.0)

# Inner band: 1σ in a, E fixed at E_ref
H_best  = H_iso.(η_vals, a_opt,   E_ref)
H_in_lo = H_iso.(η_vals, a_1σ_lo, E_ref)
H_in_hi = H_iso.(η_vals, a_1σ_hi, E_ref)

# Outer band: combined (1σ in a) + (E ∈ [20,80] GPa)
# max H: largest a, smallest E — min H: smallest a, largest E
H_out_lo = H_iso.(η_vals, a_1σ_lo, E_max)
H_out_hi = H_iso.(η_vals, a_1σ_hi, E_min)

# ── η at H = 10 km ────────────────────────────────────────────────────────────
η_pt(a, E) = E * (H_ref_km * 1e3) / (2 * a * (1 - ν^2))

η_best_pt = η_pt(a_opt,   E_ref)
η_stat_lo = η_pt(a_1σ_hi, E_ref)   # larger a → smaller η
η_stat_hi = η_pt(a_1σ_lo, E_ref)   # smaller a → larger η

σ_a = (log10(a_opt/a_1σ_lo) + log10(a_1σ_hi/a_opt)) / 2
σ_E = log10(E_max / E_ref)
σ_η = sqrt(σ_a^2 + σ_E^2)
η_comb_lo = η_best_pt * 10^(-σ_η)
η_comb_hi = η_best_pt * 10^(+σ_η)

println("At H = 10 km, E = 40 GPa:")
println("  η_best        = $(round(η_best_pt,  sigdigits=3)) Pa·s")
println("  1σ stat (a)   = [$(round(η_stat_lo, sigdigits=3)), $(round(η_stat_hi, sigdigits=3))] Pa·s")
println("  1σ combined   = [$(round(η_comb_lo, sigdigits=3)), $(round(η_comb_hi, sigdigits=3))] Pa·s")
println("  σ(log10 a) = $(round(σ_a, digits=3))   σ(log10 E) = $(round(σ_E, digits=3))   →  E uncertainty dominates")

# ── Figure ────────────────────────────────────────────────────────────────────
fig = Figure(size = (1000, 700), fontsize = 37.5)

ax = Axis(fig[1, 1],
    xlabel  = L"\eta\ \ (\mathrm{Pa\cdot s})",
    ylabel  = "H  (km)",
    title   = "Rheological trade-off — Asal",
    xticks  = (17:21, [L"10^{%$i}" for i in 17:21]),
    limits  = (17.0, 21.0, 0.0, 30.0),
    backgroundcolor = :white)

# Outer band: combined (1σ in a + E range)
band!(ax, η_log, clip(H_out_lo), clip(H_out_hi), color = (:royalblue, 0.10))

# Inner band: 1σ in a, E = E_ref (hairline for Asal — E dominates)
band!(ax, η_log, clip(H_in_lo), clip(H_in_hi), color = (:royalblue, 0.35))

# Best-fit iso-a line
lines!(ax, η_log, clip(H_best),
    color = :royalblue, linewidth = 2.5,
    label = L"a^* = %$(round(a_opt*1e5, digits=2))\times10^{-5}\ \mathrm{m\,s^{-1}}")

# Bounding dashed lines for inner band
lines!(ax, η_log, clip(H_in_lo),
    color = (:royalblue, 0.55), linewidth = 1.2, linestyle = :dash)
lines!(ax, η_log, clip(H_in_hi),
    color = (:royalblue, 0.55), linewidth = 1.2, linestyle = :dash)

# H = 10 km reference
hlines!(ax, [H_ref_km],
    color = (:black, 0.45), linewidth = 1.5, linestyle = :dash)
text!(ax, 17.08, 10.8;
    text = L"H = 10\ \mathrm{km}", fontsize = 26.25, color = :black)

# Best-fit η point
scatter!(ax, [log10(η_best_pt)], [H_ref_km],
    color = :royalblue, markersize = 12, marker = :circle, strokewidth = 0)

# Error bar: 1σ in a (matches inner band at H = 10 km)
errorbars!(ax,
    [log10(η_best_pt)], [H_ref_km],
    [log10(η_best_pt) - log10(η_stat_lo)],
    [log10(η_stat_hi)  - log10(η_best_pt)],
    direction = :x, color = :royalblue, linewidth = 2.2, whiskerwidth = 8)

# Annotation for best-fit
text!(ax, log10(η_best_pt), 11.2;
    text  = L"%$(round(η_best_pt/1e19, digits=1))\times10^{19}",
    color = :royalblue, fontsize = 26.25, align = (:center, :bottom))

# Nooner (2009) — Dabbahu, Afar: η = 5.2×10¹⁸ Pa·s, H = 13.2 km
scatter!(ax, [log10(5.2e18)], [13.2],
    color = :forestgreen, markersize = 12, marker = :diamond, strokewidth = 0)
text!(ax, log10(5.2e18) + 0.06, 14.0;
    text  = "Nooner (2009)\nDabbahu, Afar",
    color = :forestgreen, fontsize = 22.5, align = (:left, :bottom))

# Hamling et al. (2014) — REMOVED pending PDF verification

# Wright et al. (2012) — Dabbahu, Afar: η = 1–3×10¹⁸ Pa·s, H = 10 km (verified PDF)
η_wright_lo  = log10(1.0e18)   # = 18.000
η_wright_hi  = log10(3.0e18)   # ≈ 18.477
η_wright_mid = (η_wright_lo + η_wright_hi) / 2
H_wright = 10.0
scatter!(ax, [η_wright_mid], [H_wright],
    color = :crimson, markersize = 10, marker = :star5, strokewidth = 0)
errorbars!(ax,
    [η_wright_mid], [H_wright],
    [η_wright_mid - η_wright_lo],
    [η_wright_hi - η_wright_mid],
    direction = :x, color = :crimson, linewidth = 2.0, whiskerwidth = 8)
text!(ax, η_wright_hi + 0.06, H_wright;
    text  = "Wright et al. (2012)\nDabbahu, Afar",
    color = :crimson, fontsize = 22.5, align = (:left, :center))

# Cattin et al. (2005) — Asal-Ghoubbet: near-axis effective viscosity from 2D
# thermo-mechanical FEM model of post-1978 relaxation.
# η = 0.5–1 × 10¹⁸ Pa·s (near rift axis, first 6–8 yr post-rifting)
# H = 3–12 km (seismogenic thickness range: ~3 km beneath Fieale, ~12 km off-axis)
η_cat_lo = log10(0.5e18)   # ≈ 17.70
η_cat_hi = log10(1.0e18)   # = 18.00
H_cat_lo = 3.0
H_cat_hi = 12.0
poly!(ax,
    [Point2f(η_cat_lo, H_cat_lo), Point2f(η_cat_hi, H_cat_lo),
     Point2f(η_cat_hi, H_cat_hi), Point2f(η_cat_lo, H_cat_hi)],
    color = (:darkorange, 0.20),
    strokecolor = :darkorange, strokewidth = 1.8)
text!(ax, (η_cat_lo + η_cat_hi) / 2, H_cat_hi + 0.7;
    text  = "Cattin et al. (2005)\nAsal, near-axis",
    color = :darkorange, fontsize = 22.5, align = (:center, :bottom))

# Shading note
text!(ax, 20.7, 1.5;
    text  = "dark: 1σ in a  (E = 40 GPa) |  light: combined (a + E ∈ [20, 80] GPa)",
    color = (:gray, 0.70), fontsize = 22.5, align = (:right, :bottom))

axislegend(ax, position = :lt, framevisible = true, labelsize = 28.12)

save("fig/fig_rheology_1sigma_asal.png", fig)
println("Saved  fig/fig_rheology_1sigma_asal.png")
