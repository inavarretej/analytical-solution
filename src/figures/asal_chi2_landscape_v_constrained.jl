using CairoMakie, CSV, DataFrames, Printf, ColorSchemes, LaTeXStrings, Dates

include("asal_v_constrained_common.jl")   # data, V cases, predictions, joint_chi2

# Joint (time series + profile) chi^2 landscapes in (a, b) for the
# V-constrained Asal fits (src/figures/asal_css_v_constrained.jl, which must
# be run first -- its best-fit (a, b) seed the grids): EV and EVE RSS, for
# the literature V. T = b/V and V_pl = V/2 are derived at every grid point.
#
# Plotted as Δχ² = χ² − χ²_min so both panels share one colour scale;
# contours at Δχ² = 2.30 / 6.17 / 11.83 (68.3 / 95.4 / 99.7 %, 2 parameters),
# equivalent to the χ²/N levels of src/figures/asal_chi2_landscape_css.jl.
# The region b < V·T_min_no_repeat (T below the no-repeat floor) is masked,
# with the floor drawn as a dashed line.
mkpath("fig")
mkpath("reports")
println("Julia threads: $(Threads.nthreads())")

params = CSV.read(params_csv, DataFrame)
best(Vmm, model) = only(eachrow(filter(r -> r.V_mm_yr == Vmm && r.model == model, params)))

const DCHI2 = [2.30, 6.17, 11.83]

function chi2_grid(chi2_fn, a_vec, b_vec, V)
    G = fill(Inf, length(a_vec), length(b_vec))
    Threads.@threads for j in eachindex(b_vec)
        b = b_vec[j]
        b / V < T_min_no_repeat && continue
        for i in eachindex(a_vec)
            G[i, j] = chi2_fn(a_vec[i], b)
        end
    end
    G
end

# Same self-recentring zoom as src/figures/asal_chi2_landscape_css.jl: centre
# on the grid's own argmin and widen until the 3σ contour no longer touches
# the grid edge (the masked floor does not count as an edge).
function zoomed_grid(chi2_fn, a_opt, b_opt, V; n_search, n_final, max_iter = 8)
    a_half = 0.2
    b_half = max(0.15 * b_opt, 0.01)
    grid(n) = (10 .^ range(log10(a_opt) - a_half, log10(a_opt) + a_half, length = n),
               collect(range(max(1e-4, b_opt - b_half), b_opt + b_half, length = n)))
    for _ in 1:max_iter
        a_vec, b_vec = grid(n_search)
        G = chi2_grid(chi2_fn, a_vec, b_vec, V)
        imin = argmin(G)
        a_opt, b_opt = a_vec[imin[1]], b_vec[imin[2]]
        m3 = G .≤ minimum(G) + DCHI2[3]
        (any(m3[1, :]) || any(m3[end, :]) || any(m3[:, 1]) || any(m3[:, end])) || break
        a_half *= 1.8
        b_half *= 1.8
    end
    a_vec, b_vec = grid(n_final)
    G = chi2_grid(chi2_fn, a_vec, b_vec, V)
    imin = argmin(G)
    return a_vec, b_vec, G, minimum(G), a_vec[imin[1]], b_vec[imin[2]]
end

models = [(name = "EV",  title = "EV RSS",  n_search = 45, n_final = 300,
           chi2 = V -> (a, b) -> joint_chi2_raw(pred_ts_ev,  pred_prof_ev,  a, b, V)),
          (name = "EVE", title = "EVE RSS", n_search = 21, n_final = 71,
           chi2 = V -> (a, b) -> joint_chi2_raw(pred_ts_eve, pred_prof_eve, a, b, V))]

panels = []
for vc in V_cases, m in models
    Vmm = round(Int, vc.V * yr * 1000)
    p   = best(Vmm, m.name)
    println("Building $(m.title) landscape, V = $Vmm mm/yr …")
    t = @elapsed a_vec, b_vec, G, cmin, a_g, b_g = zoomed_grid(m.chi2(vc.V), p.a, p.b, vc.V;
                                                     n_search = m.n_search, n_final = m.n_final)
    @printf("  χ²_min = %.3f (grid) vs %.3f (fit)  at a = %s, b = %.3f  [%.0f s]\n", cmin, p.chi2, sci(a_g), b_g, t)
    push!(panels, (; m.name, m.title, Vmm, V = vc.V, a_vec, b_vec, G, cmin, a_g, b_g, a_fit = p.a, b_fit = p.b))
end

# ── 1σ projected ranges (Δχ² ≤ 2.30) ─────────────────────────────────────────
function one_sigma(pn)
    m1 = pn.G .≤ pn.cmin + DCHI2[1]
    as = pn.a_vec[vec(any(m1, dims = 2))]
    bs = pn.b_vec[vec(any(m1, dims = 1))]
    (a = extrema(as), b = extrema(bs), T = extrema(bs) ./ pn.V ./ yr,
     tauR = reverse(tau_R_from_a.(extrema(as))) ./ yr)
end

open("reports/asal_chi2_landscape_v_constrained.txt", "w") do io
    println(io, "Asal — V-constrained RSS χ² landscapes (joint time series + profile)")
    println(io, "Generated: $(Dates.now())")
    println(io, "1σ ranges: projections of the Δχ² ≤ 2.30 region onto each axis (2 free parameters)")
    println(io)
    for pn in panels
        s = one_sigma(pn)
        println(io, "=== $(pn.title), V = $(pn.Vmm) mm/yr ===")
        @printf(io, "  χ²_min = %.3f  at a = %s m/s, b = %.3f m (grid argmin)\n", pn.cmin, sci(pn.a_g), pn.b_g)
        @printf(io, "  a  ∈ [%s, %s] m/s\n", sci(s.a[1]), sci(s.a[2]))
        @printf(io, "  b  ∈ [%.3f, %.3f] m\n", s.b...)
        @printf(io, "  T  ∈ [%.1f, %.1f] yr   (T = b/V)\n", s.T...)
        pn.name == "EVE" && @printf(io, "  τ_R ∈ [%.2f, %.2f] yr (derived from a)\n", s.tauR...)
        any(isinf, pn.G) && minimum(pn.b_vec[vec(any(pn.G .≤ pn.cmin + DCHI2[1], dims = 1))]) <=
            pn.V * T_min_no_repeat * 1.02 && println(io, "  NOTE: 1σ region reaches the T ≥ T_min_no_repeat floor")
        println(io)
    end
end
println("reports/asal_chi2_landscape_v_constrained.txt  saved")

# ── Figure: rows = V, columns = EV / EVE ─────────────────────────────────────
conf_colors = [:white, :gold, :tomato]
clip_hi     = DCHI2[3] * 1.4
cmap        = cgrad(ColorSchemes.bamako100, 20, categorical = true)

fig = Figure(size = (2000, 900 * length(V_cases) + 200), fontsize = 30)
hms = []
for (k, pn) in enumerate(panels)
    row, col = fldmod1(k, 2)
    ax = Axis(fig[row, col], xscale = log10,
        xlabel = L"a\ \ (\mathrm{m\,s^{-1}})", ylabel = L"b\ \ (\mathrm{m})",
        title = "$(pn.title), V = $(pn.Vmm) mm/yr  (T = b/V)")
    D = replace(pn.G .- pn.cmin, Inf => NaN)
    hm = heatmap!(ax, pn.a_vec, pn.b_vec, D, colormap = cmap, colorrange = (0, clip_hi),
                  highclip = cmap[end], nan_color = (:gray, 0.25))
    push!(hms, hm)
    for (lev, c) in zip(DCHI2, conf_colors)
        contour!(ax, pn.a_vec, pn.b_vec, D, levels = [lev], color = c, linewidth = 2.5)
    end
    b_floor = pn.V * T_min_no_repeat
    b_floor >= pn.b_vec[1] && hlines!(ax, [b_floor], color = :black, linestyle = :dash, linewidth = 1.5)
    scatter!(ax, [pn.a_g], [pn.b_g], color = :cyan, marker = :star5, markersize = 22, strokewidth = 0)
    scatter!(ax, [pn.a_fit], [pn.b_fit], color = :magenta, marker = :xcross, markersize = 16)
end

Colorbar(fig[1:length(V_cases), 3], hms[1], label = L"\Delta\chi^2 = \chi^2 - \chi^2_{\min}",
    ticks = ([0; DCHI2], ["0", "1σ", "2σ", "3σ"]))
# The T-floor entry only appears if some panel actually reaches the floor.
floor_shown = any(pn -> pn.V * T_min_no_repeat >= pn.b_vec[1], panels)
leg_elems = vcat([LineElement(color = c, linewidth = 2.5) for c in conf_colors],
                 [MarkerElement(color = :cyan, marker = :star5, markersize = 22),
                  MarkerElement(color = :magenta, marker = :xcross, markersize = 16)])
leg_labels = ["68.3% (1σ)", "95.4% (2σ)", "99.7% (3σ)", "grid argmin", "optimizer best fit"]
if floor_shown
    push!(leg_elems, LineElement(color = :black, linestyle = :dash, linewidth = 1.5))
    push!(leg_labels, "T = T_min (no repeat intrusion)")
end
Legend(fig[length(V_cases) + 1, 1:2], leg_elems, leg_labels,
    orientation = :horizontal, nbanks = 2, framevisible = true, labelsize = 24)

save("fig/fig_chi2_landscape_asal_v_constrained.png", fig, px_per_unit = 2)
println("fig/fig_chi2_landscape_asal_v_constrained.png saved")
