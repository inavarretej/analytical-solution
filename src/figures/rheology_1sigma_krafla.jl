using CairoMakie, LaTeXStrings, Statistics

mkpath("fig")

# ── Constants ─────────────────────────────────────────────────────────────────
ν        = 0.25
E_ref    = 40e9
E_min    = 20e9
E_max    = 80e9
H_ref_km = 10.0

# ── Best-fit a and 1σ bounds ─────────────────────────────────────────────────
# From reports/krafla_chi2_landscape_v_constrained.txt (written by
# src/figures/krafla_v_constrained_combined.jl): EV RSS, joint fit to both
# Krafla datasets, T = b/V (V = 1.8 cm/yr), common time origin 1975;
# 1σ = Δχ² ≤ 2.30 projected onto a. (Previously a* = 5.13e-5, 1σ 2.88–8.37e-5,
# from the older free-T fit in fig_chi2_landscape_krafla.jl.)
a_opt   = 4.989e-5    # m/s  — best fit
a_1σ_lo = 4.21e-5     # m/s  — lower 1σ edge → largest η
a_1σ_hi = 5.82e-5     # m/s  — upper 1σ edge → smallest η

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
println("  σ(log10 a) = $(round(σ_a, digits=3))   σ(log10 E) = $(round(σ_E, digits=3))")

# ── Figure ────────────────────────────────────────────────────────────────────
fig = Figure(size = (1600, 700), fontsize = 37.5)

ax = Axis(fig[1, 1],
    xlabel  = L"\eta\ \ (\mathrm{Pa\cdot s})",
    ylabel  = "H  (km)",
    title   = "Rheological trade-off — Krafla",
    xticks  = (17:21, [L"10^{%$i}" for i in 17:21]),
    limits  = (17.0, 21.0, 0.0, 30.0),
    backgroundcolor = :white)

# Outer band: combined (1σ in a + E range)
band!(ax, η_log, clip(H_out_lo), clip(H_out_hi), color = (:firebrick, 0.10))

# Inner band: 1σ in a, E = E_ref
band!(ax, η_log, clip(H_in_lo), clip(H_in_hi), color = (:firebrick, 0.35))

# Best-fit iso-a line
lines!(ax, η_log, clip(H_best),
    color = :firebrick, linewidth = 2.5,
    label = L"a^* = %$(round(a_opt*1e5, digits=2))\times10^{-5}\ \mathrm{m\,s^{-1}}")

lines!(ax, η_log, clip(H_in_lo),
    color = (:firebrick, 0.55), linewidth = 1.2, linestyle = :dash)
lines!(ax, η_log, clip(H_in_hi),
    color = (:firebrick, 0.55), linewidth = 1.2, linestyle = :dash)

# H = 10 km reference
hlines!(ax, [H_ref_km],
    color = (:black, 0.45), linewidth = 1.5, linestyle = :dash)
text!(ax, 17.08, 10.8;
    text = L"H = 10\ \mathrm{km}", fontsize = 26.25, color = :black)

# Best-fit η point
scatter!(ax, [log10(η_best_pt)], [H_ref_km],
    color = :firebrick, markersize = 12, marker = :circle, strokewidth = 0)

errorbars!(ax,
    [log10(η_best_pt)], [H_ref_km],
    [log10(η_best_pt) - log10(η_stat_lo)],
    [log10(η_stat_hi)  - log10(η_best_pt)],
    direction = :x, color = :firebrick, linewidth = 2.2, whiskerwidth = 8)

text!(ax, log10(η_best_pt), 11.2;
    text  = L"%$(round(η_best_pt/1e18, digits=1))\times10^{18}",
    color = :firebrick, fontsize = 26.25, align = (:center, :bottom))

# ── Literature (verified against PDFs, reports/reference_verification.md) ──────
# Labelled through the legend (outside the axis) instead of in-plot text, so
# the overlapping boxes/markers stay readable.

# Foulger et al. (1992) — η = 0.3–2e19 Pa·s for h = 8–30 km (elastic layer),
# b = 5–10 km thin viscous channel, κ = 1.1 m²/s (Nature 358, p. 490).
# Heki et al. (1993) — same layer assumptions, 2-D model with κ ≈ 10 m²/s
# gives the revised η = 0.3–2e18 Pa·s (abstract; p. 14,293).
# Both are channel (Elsasser) viscosities, not half-space viscosities like ours,
# so they are drawn as open boxes rather than points.
for (η_lo, η_hi, col, lab) in ((0.3e19, 2.0e19, :darkorange, "Foulger et al. (1992)"),
                               (0.3e18, 2.0e18, :teal,       "Heki et al. (1993)"))
    poly!(ax, Rect(log10(η_lo), 8.0, log10(η_hi) - log10(η_lo), 22.0),
        color = (col, 0.12), strokecolor = col, strokewidth = 2.0, linestyle = :dash,
        label = lab * "\nchannel model, h = 8–30 km")
end

# Hofton & Foulger (1996) — η = 1.1e18 Pa·s, H = 10 km (Maxwell half-space)
η_hf = 1.1e18
scatter!(ax, [log10(η_hf)], [H_ref_km],
    color = :white, strokecolor = :darkslategray, strokewidth = 2.0,
    markersize = 16, marker = :rect, label = "Hofton & Foulger (1996)")

# Wright et al. (2012) — η = 1–3e18 Pa·s, H ≈ 10 km. Review: these values
# cite Hofton & Foulger (1996) and Heki et al. (1993), not independent.
# Plotted at H = 11 km only to avoid overlapping the Hofton & Foulger marker.
η_wr_lo  = log10(1.0e18)
η_wr_hi  = log10(3.0e18)
η_wr_mid = (η_wr_lo + η_wr_hi) / 2
H_wright = 11.0
errorbars!(ax,
    [η_wr_mid], [H_wright],
    [η_wr_mid - η_wr_lo],
    [η_wr_hi - η_wr_mid],
    direction = :x, color = :darkslategray, linewidth = 2.0, whiskerwidth = 8)
scatter!(ax, [η_wr_mid], [H_wright],
    color = :white, strokecolor = :darkslategray, strokewidth = 2.0,
    markersize = 16, marker = :circle, label = "Wright et al. (2012), review")

# Sigmundsson (1991) — upper bound η ≤ 1e19 Pa·s (post-glacial rebound), H = 10 km
η_sig = 1.0e19
scatter!(ax, [log10(η_sig)], [H_ref_km],
    color = :white, strokecolor = :darkslategray, strokewidth = 2.0,
    markersize = 18, marker = :ltriangle, label = "Sigmundsson (1991), η ≤ 10¹⁹")

# Shading note
text!(ax, 20.9, 1.5;
    text  = "dark: 1σ in a  (E = 40 GPa)\nlight: a + E ∈ [20, 80] GPa",
    color = (:gray, 0.70), fontsize = 22.5, align = (:right, :bottom))

Legend(fig[1, 2], ax, framevisible = true, labelsize = 24, patchsize = (40, 30))

save("fig/fig_rheology_1sigma.png", fig)
println("Saved  fig/fig_rheology_1sigma.png")
