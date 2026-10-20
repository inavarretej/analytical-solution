using CairoMakie, CSV, DataFrames, Printf, ColorSchemes, LaTeXStrings

include("asal_v_constrained_common.jl")   # data, V, predictions (EVE: converged split forms), joint_chi2_raw

# Asal, V-constrained RSS (T = b/V, V = 13 mm/yr): one figure, four panels,
# same layout as src/figures/krafla_v_constrained_combined.jl.
#   (a) EP00–DF00 baseline time series, (b) velocity profile at τ = 19 yr,
#       with the joint EV and EVE RSS fits;
#   (c), (d) joint χ² landscapes in (a, b) for EV and EVE, as Δχ² with
#       1σ/2σ/3σ contours (2 free parameters). The region T < T_min (no repeat
#       intrusion 1978–2016) is masked.
#
# Best fits are read from reports/asal_css_v_constrained_params.csv
# (src/figures/asal_css_v_constrained.jl, run that first). The landscape grid
# minimum is checked against them.
mkpath("fig")
println("Julia threads: $(Threads.nthreads())")

V      = only(V_cases).V
params = CSV.read(params_csv, DataFrame)
best(model) = only(eachrow(filter(r -> r.model == model, params)))
ev_p, eve_p = best("EV"), best("EVE")

const DCHI2 = [2.30, 6.17, 11.83]

function chi2_grid(chi2_fn, a_vec, b_vec)
    G = fill(Inf, length(a_vec), length(b_vec))
    Threads.@threads for j in eachindex(b_vec)
        b_vec[j] / V < T_min_no_repeat && continue
        for i in eachindex(a_vec)
            G[i, j] = chi2_fn(a_vec[i], b_vec[j])
        end
    end
    G
end

# Centre on the grid argmin and widen until the 3σ contour clears every edge
# (the masked T floor does not count as an edge).
function zoomed_grid(chi2_fn, a0, b0; n_search = 41, n_final = 161, max_iter = 8)
    a_half, b_half = 0.2, max(0.15 * b0, 0.01)
    grid(n) = (10 .^ range(log10(a0) - a_half, log10(a0) + a_half, length = n),
               collect(range(max(1e-4, b0 - b_half), b0 + b_half, length = n)))
    for _ in 1:max_iter
        a_vec, b_vec = grid(n_search)
        G = chi2_grid(chi2_fn, a_vec, b_vec)
        imin = argmin(G)
        a0, b0 = a_vec[imin[1]], b_vec[imin[2]]
        m3 = G .≤ minimum(G) + DCHI2[3]
        (any(m3[1, :]) || any(m3[end, :]) || any(m3[:, 1]) || any(m3[:, end])) || break
        a_half *= 1.8
        b_half *= 1.8
    end
    # Crop to the 3σ region plus a margin of its own size on each side, so the
    # contours fill the panel (the Asal minimum is far narrower than 0.2 dex).
    a_vec, b_vec = grid(n_search)
    G  = chi2_grid(chi2_fn, a_vec, b_vec)
    m3 = G .≤ minimum(G) + DCHI2[3]
    la = extrema(log10.(a_vec[vec(any(m3, dims = 2))]))
    lb = extrema(b_vec[vec(any(m3, dims = 1))])
    da = max(la[2] - la[1], 2 * step(range(la..., length = 2)) / n_search)
    db = lb[2] - lb[1]
    a_vec = 10 .^ range(la[1] - 0.6da, la[2] + 0.6da, length = n_final)
    b_vec = collect(range(max(1e-4, lb[1] - 0.6db), lb[2] + 0.6db, length = n_final))
    G = chi2_grid(chi2_fn, a_vec, b_vec)
    (; a_vec, b_vec, G, cmin = minimum(G))
end

models = [(name = "EV",  title = "EV RSS",  p = ev_p,
           chi2 = (a, b) -> joint_chi2_raw(pred_ts_ev,  pred_prof_ev,  a, b, V)),
          (name = "EVE", title = "EVE RSS", p = eve_p,
           chi2 = (a, b) -> joint_chi2_raw(pred_ts_eve, pred_prof_eve, a, b, V))]

panels = map(models) do m
    println("Building $(m.title) landscape …")
    t = @elapsed L = zoomed_grid(m.chi2, m.p.a, m.p.b)
    @printf("  χ²_min = %.3f (grid) vs %.3f (fit)  [%.0f s]\n", L.cmin, m.p.chi2, t)
    (; m.name, m.title, L..., a_fit = m.p.a, b_fit = m.p.b, chi2_fit = m.p.chi2)
end

# ── Figure ───────────────────────────────────────────────────────────────────
fig = Figure(size = (2600, 2200), fontsize = 40)

t_yr_grid = range(0.05, stop = 26.0, length = 400)
t_s_grid  = t_yr_grid .* yr
x_grid_km = collect(range(-12, 12, length = 600))
x_grid_m  = x_grid_km .* 1e3

ax_a = Axis(fig[1, 1], xlabel = "Year", ylabel = "EP00–DF00 distance change (mm)",
            title = "Vigny et al. (2007), baseline EP00–DF00", limits = (t0_year - 1, 2005, -80, 580))
ax_b = Axis(fig[1, 2], xlabel = "Distance from rift axis (km)",
            ylabel = "Velocity, Somalia frame (mm yr⁻¹)",
            title = "Vigny et al. (2007), profile, τ = 19 yr", limits = (-12, 12, -2, 22))

lines!(ax_a, t_yr_grid .+ t0_year, tref(t -> pred_ts_ev(t, ev_p.a, ev_p.b, V), t_s_grid),
       color = :purple, linewidth = 3, label = "EV RSS")
lines!(ax_a, t_yr_grid .+ t0_year, tref(t -> pred_ts_eve(t, eve_p.a, eve_p.b, V), t_s_grid),
       color = :navy, linewidth = 3, label = "EVE RSS")
errorbars!(ax_a, obs_year, obs_dL, obs_sig, color = :black)
scatter!(ax_a, obs_year, obs_dL, color = :black, markersize = 12, label = "Vigny et al. (2007)")
axislegend(ax_a, position = :rb, framevisible = true, labelsize = 30)

lines!(ax_b, x_grid_km, pred_prof_ev(x_grid_m, ev_p.a, ev_p.b, V), color = :purple, linewidth = 3)
lines!(ax_b, x_grid_km, pred_prof_eve(x_grid_m, eve_p.a, eve_p.b, V), color = :navy, linewidth = 3)
vlines!(ax_b, [0.0], color = (:gray, 0.6), linestyle = :dot, linewidth = 2)
vlines!(ax_b, [x_center_m / 1000], color = (:black, 0.5), linestyle = :dash, linewidth = 2)
errorbars!(ax_b, prof_x_km, prof_v, prof_sig, color = :black, whiskerwidth = 8)
scatter!(ax_b, prof_x_km, prof_v, color = :black, markersize = 12)

conf_colors = [:white, :gold, :tomato]
clip_hi     = DCHI2[3] * 1.4
cmap        = cgrad(ColorSchemes.bamako100, 20, categorical = true)
b_floor     = V * T_min_no_repeat

hm = nothing
ax_land = map(enumerate(panels)) do (k, pn)
    ax = Axis(fig[2, k], xscale = log10, xtickformat = vs -> [@sprintf("%.2e", v) for v in vs],
              xlabel = L"a\ \ (\mathrm{m\,s^{-1}})", ylabel = L"b\ \ (\mathrm{m})",
              title = "$(pn.title), Δχ² (T = b/V)")
    cmin = min(pn.cmin, pn.chi2_fit)
    D = replace(pn.G .- cmin, Inf => NaN)
    global hm = heatmap!(ax, pn.a_vec, pn.b_vec, D, colormap = cmap, colorrange = (0, clip_hi),
                         highclip = cmap[end], nan_color = (:gray, 0.25))
    for (lev, c) in zip(DCHI2, conf_colors)
        contour!(ax, pn.a_vec, pn.b_vec, D, levels = [lev], color = c, linewidth = 3)
    end
    b_floor >= pn.b_vec[1] && hlines!(ax, [b_floor], color = :black, linestyle = :dash, linewidth = 2)
    scatter!(ax, [pn.a_fit], [pn.b_fit], color = :magenta, marker = :xcross, markersize = 26)
    axT = Axis(fig[2, k], yaxisposition = :right, ylabel = L"T = b/V\ \ (\mathrm{yr})",
               limits = (extrema(pn.a_vec), (pn.b_vec[1], pn.b_vec[end]) ./ V ./ yr),
               xscale = log10)
    hidexdecorations!(axT); hidespines!(axT)
    ax
end

Colorbar(fig[2, 3], hm, label = L"\Delta\chi^2 = \chi^2 - \chi^2_{\min}",
         ticks = ([0; DCHI2], ["0", "1σ", "2σ", "3σ"]))
leg_elems  = vcat([LineElement(color = c, linewidth = 3) for c in conf_colors],
                  [MarkerElement(color = :magenta, marker = :xcross, markersize = 26)])
leg_labels = ["68.3% (1σ)", "95.4% (2σ)", "99.7% (3σ)", "best fit"]
if any(pn -> b_floor >= pn.b_vec[1], panels)
    push!(leg_elems, LineElement(color = :black, linestyle = :dash, linewidth = 2))
    push!(leg_labels, "T = T_min (no repeat intrusion)")
end
Legend(fig[3, 1:2], leg_elems, leg_labels, orientation = :horizontal, framevisible = true,
       labelsize = 32, backgroundcolor = (:gray60, 1.0))

for (ax, lab) in zip([ax_a, ax_b, ax_land...], ["(a)", "(b)", "(c)", "(d)"])
    text!(ax, 0.02, 0.98, text = lab, space = :relative, align = (:left, :top), fontsize = 44)
end

save("fig/fig_asal_v_constrained_combined.png", fig, px_per_unit = 2)
println("fig/fig_asal_v_constrained_combined.png saved")
