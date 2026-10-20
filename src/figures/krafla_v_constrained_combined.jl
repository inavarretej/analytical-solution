using CairoMakie, CSV, DataFrames, Printf, ColorSchemes, LaTeXStrings, Dates, Optim

include("../models/viscous.jl")
include("../models/maxwell.jl")
include("../parameters.jl")   # H, E_L, nu_L

# Krafla, V-constrained RSS (T = b/V, V = 1.8 cm/yr, t0 = 1975): one figure,
# four panels.
#   (a), (b): joint EV and EVE RSS fits, Foulger et al. (1992) and
#             Hofton & Foulger (1996) profiles.
#   (c), (d): joint chi^2 landscapes in (a, b) for EV and EVE, as
#             Δχ² = χ² − χ²_min with 1σ/2σ/3σ contours (2 free parameters).
#
# The grids are seeded with the best fits of
# reports/krafla_css_v_constrained_params.csv (src/figures/krafla_css_v_constrained.jl);
# the best fit of each model is then re-derived here as the grid argmin
# refined by Nelder-Mead, and written to
# reports/krafla_v_constrained_combined_params.csv. The Hofton
# offset c is linear in χ², so it is solved analytically at every (a, b), as
# in the fit itself. The landscapes are therefore profiles over c.
mkpath("fig")
mkpath("reports")
println("Julia threads: $(Threads.nthreads())")
sci(x) = @sprintf("%.2e", x)

yr    = 3600.0 * 24 * 365
V_lit = 0.018 / yr   # full opening rate (DeMets et al. 1994); far field ±V/2
T0    = 1975

tau_R_from_a(a) = H / (a * (1 - nu_L))

# ── Data ─────────────────────────────────────────────────────────────────────
foulger = CSV.read("data/extracted_data.csv", DataFrame)
hofton  = CSV.read("data/hofton.csv", DataFrame)
xf, uf, sf = foulger.X .* 1e3, foulger.Y ./ 100, foulger.SD ./ 100
xh, uh, sh = hofton.X  .* 1e3, hofton.Y  ./ 100, hofton.SD  ./ 100
tf = ((1987 - T0) * yr, (1990 - T0) * yr)
th = ((1987 - T0) * yr, (1992 - T0) * yr)

params = CSV.read("reports/krafla_css_v_constrained_params.csv", DataFrame)
best(model) = only(eachrow(filter(r -> r.model == model, params)))

# ── Forward models (same as krafla_css_v_constrained.jl) ─────────────────────
pred_ev(x, ti, tf_, a, b) = diff_css(x, ti, tf_, b, V_lit, a)
# EVE uses the converged split form. get_u_css_maxwell_closed at its default
# k_limit_scale = 1 (used by krafla_css_v_constrained.jl) under-truncates the
# k-integral at these T and biases χ² low by ~9 at the old best fit
# (58.5 vs 67.3); the split form agrees with k_limit_scale = 100 to ~0.1 in χ².
pred_eve(x, ti, tf_, a, b) = begin
    tR = tau_R_from_a(a)
    get_u_css_maxwell_split.(x, tf_, b, V_lit, a, tR) .-
    get_u_css_maxwell_split.(x, ti,  b, V_lit, a, tR)
end

const wh = 1 ./ sh .^ 2
function joint_chi2(pred, a, b)
    pf = pred(xf, tf..., a, b)
    ph = pred(xh, th..., a, b)
    c  = sum(wh .* (uh .- ph)) / sum(wh)
    (sum(((uf .- pf) ./ sf) .^ 2) + sum(((uh .- ph .- c) ./ sh) .^ 2), c)
end

# ── Landscapes ───────────────────────────────────────────────────────────────
const DCHI2 = [2.30, 6.17, 11.83]

function chi2_grid(pred, a_vec, b_vec)
    G = Matrix{Float64}(undef, length(a_vec), length(b_vec))
    Threads.@threads for j in eachindex(b_vec)
        for i in eachindex(a_vec)
            G[i, j] = joint_chi2(pred, a_vec[i], b_vec[j])[1]
        end
    end
    G
end

# Centre on the grid argmin and widen until the 3σ contour clears every edge.
function zoomed_grid(pred, a0, b0; n_search = 41, n_final = 201, max_iter = 8)
    a_half, b_half = 0.2, 0.1 * b0
    grid(n) = (10 .^ range(log10(a0) - a_half, log10(a0) + a_half, length = n),
               collect(range(max(0.05, b0 - b_half), b0 + b_half, length = n)))
    for _ in 1:max_iter
        a_vec, b_vec = grid(n_search)
        G = chi2_grid(pred, a_vec, b_vec)
        imin = argmin(G)
        a0, b0 = a_vec[imin[1]], b_vec[imin[2]]
        m3 = G .≤ minimum(G) + DCHI2[3]
        (any(m3[1, :]) || any(m3[end, :]) || any(m3[:, 1]) || any(m3[:, end])) || break
        a_half *= 1.6
        b_half *= 1.6
    end
    a_vec, b_vec = grid(n_final)
    G = chi2_grid(pred, a_vec, b_vec)
    imin = argmin(G)
    (; a_vec, b_vec, G, cmin = minimum(G), a_g = a_vec[imin[1]], b_g = b_vec[imin[2]])
end

models = [(name = "EV",  title = "EV RSS",  pred = pred_ev),
          (name = "EVE", title = "EVE RSS", pred = pred_eve)]

N = length(xf) + length(xh)
panels = map(models) do m
    p = best(m.name)
    println("Building $(m.title) landscape …")
    t = @elapsed L = zoomed_grid(m.pred, p.a, p.b)
    @printf("  χ²_min = %.3f (grid) vs %.3f (fit)  at a = %s, b = %.3f  [%.0f s]\n",
            L.cmin, p.chi2, sci(L.a_g), L.b_g, t)
    res = optimize(q -> joint_chi2(m.pred, 10^q[1], q[2])[1], [log10(L.a_g), L.b_g], NelderMead())
    a_fit, b_fit = 10^Optim.minimizer(res)[1], Optim.minimizer(res)[2]
    chi2_fit, c_fit = joint_chi2(m.pred, a_fit, b_fit)
    chi2_fit > L.cmin && ((a_fit, b_fit, chi2_fit) = (L.a_g, L.b_g, L.cmin); c_fit = joint_chi2(m.pred, a_fit, b_fit)[2])
    @printf("  best fit: a = %s, b = %.3f, c = %.3f, χ² = %.3f (χ²_r = %.3f)  [old fit χ² = %.3f]\n",
            sci(a_fit), b_fit, c_fit, chi2_fit, chi2_fit / (N - 3), p.chi2)
    (; m.name, m.title, L..., a_fit, b_fit, c_fit, chi2_fit)
end
CSV.write("reports/krafla_v_constrained_combined_params.csv", DataFrame(
    model = [pn.name for pn in panels], a = [pn.a_fit for pn in panels], b = [pn.b_fit for pn in panels],
    c_hofton = [pn.c_fit for pn in panels], chi2 = [pn.chi2_fit for pn in panels],
    chi2_r = [pn.chi2_fit / (N - 3) for pn in panels],
    T_yr = [pn.b_fit / V_lit / yr for pn in panels],
    tau_R_yr = [pn.name == "EVE" ? tau_R_from_a(pn.a_fit) / yr : NaN for pn in panels]))

# ── 1σ projected ranges ──────────────────────────────────────────────────────
open("reports/krafla_chi2_landscape_v_constrained.txt", "w") do io
    println(io, "Krafla — V-constrained RSS χ² landscapes (joint Foulger 1992 + Hofton & Foulger 1996)")
    println(io, "Generated: $(Dates.now())")
    println(io, "V = 1.8 cm/yr (full rate), T = b/V, t0 = $T0, N = $N, Hofton offset c solved analytically")
    println(io, "1σ ranges: projections of the Δχ² ≤ 2.30 region onto each axis (2 parameters in the landscape)")
    println(io)
    for pn in panels
        m1 = pn.G .≤ min(pn.cmin, pn.chi2_fit) + DCHI2[1]
        as = extrema(pn.a_vec[vec(any(m1, dims = 2))])
        bs = extrema(pn.b_vec[vec(any(m1, dims = 1))])
        println(io, "=== $(pn.title) ===")
        @printf(io, "  best fit: a = %s m/s, b = %.3f m, c = %.3f m, χ² = %.3f, χ²_r = %.3f (N - 3 dof)\n",
                sci(pn.a_fit), pn.b_fit, pn.c_fit, pn.chi2_fit, pn.chi2_fit / (N - 3))
        @printf(io, "  T = %.1f yr", pn.b_fit / V_lit / yr)
        pn.name == "EVE" && @printf(io, "   τ_R = %.2f yr", tau_R_from_a(pn.a_fit) / yr)
        println(io)
        @printf(io, "  a  ∈ [%s, %s] m/s\n", sci(as[1]), sci(as[2]))
        @printf(io, "  b  ∈ [%.3f, %.3f] m\n", bs...)
        @printf(io, "  T  ∈ [%.1f, %.1f] yr   (T = b/V)\n", (bs ./ V_lit ./ yr)...)
        pn.name == "EVE" && @printf(io, "  τ_R ∈ [%.2f, %.2f] yr (derived from a)\n",
                                    (reverse(tau_R_from_a.(as)) ./ yr)...)
        println(io)
    end
end
println("reports/krafla_chi2_landscape_v_constrained.txt  saved")

# ── Figure ───────────────────────────────────────────────────────────────────
ev, eve = panels
fig = Figure(size = (2600, 2200), fontsize = 40)

x_km = range(-180, stop = 180, length = 1000)
x_m  = x_km .* 1e3

function fit_panel!(pos, title, xo, uo, so, ti, tf_, c_ev, c_eve, data_label)
    ax = Axis(pos, xlabel = "Distance from Krafla volcanic system (km)",
              ylabel = "Displacement (cm)", limits = (-180, 180, -15, 20), title = title)
    lines!(ax, x_km, (pred_ev(x_m, ti, tf_, ev.a_fit, ev.b_fit) .+ c_ev) .* 100,
           color = :purple, linewidth = 3, label = "EV RSS")
    lines!(ax, x_km, (pred_eve(x_m, ti, tf_, eve.a_fit, eve.b_fit) .+ c_eve) .* 100,
           color = :navy, linewidth = 3, label = "EVE RSS")
    errorbars!(ax, xo ./ 1e3, uo .* 100, so .* 100, color = :black)
    scatter!(ax, xo ./ 1e3, uo .* 100, color = :black, markersize = 12, label = data_label)
    axislegend(ax, position = :rb, framevisible = true, labelsize = 30)
    ax
end

ax_a = fit_panel!(fig[1, 1], "Foulger et al. (1992), τ = 12–15 yr", xf, uf, sf, tf..., 0.0, 0.0,
                  "Foulger et al. (1992)")
ax_b = fit_panel!(fig[1, 2], "Hofton & Foulger (1996), τ = 12–17 yr", xh, uh, sh, th..., ev.c_fit, eve.c_fit,
                  "Hofton & Foulger (1996)")

conf_colors = [:white, :gold, :tomato]
clip_hi     = DCHI2[3] * 1.4
cmap        = cgrad(ColorSchemes.bamako100, 20, categorical = true)

hm = nothing
ax_land = map(enumerate(panels)) do (k, pn)
    ax = Axis(fig[2, k], xscale = log10, xtickformat = vs -> [@sprintf("%.1e", v) for v in vs],
              xlabel = L"a\ \ (\mathrm{m\,s^{-1}})", ylabel = L"b\ \ (\mathrm{m})",
              title = "$(pn.title), Δχ² (T = b/V)")
    D = pn.G .- min(pn.cmin, pn.chi2_fit)
    global hm = heatmap!(ax, pn.a_vec, pn.b_vec, D, colormap = cmap, colorrange = (0, clip_hi),
                         highclip = cmap[end])
    for (lev, c) in zip(DCHI2, conf_colors)
        contour!(ax, pn.a_vec, pn.b_vec, D, levels = [lev], color = c, linewidth = 3)
    end
    scatter!(ax, [pn.a_fit], [pn.b_fit], color = :magenta, marker = :xcross, markersize = 26)
    # Right axis: T = b/V for the same b range
    axT = Axis(fig[2, k], yaxisposition = :right, ylabel = L"T = b/V\ \ (\mathrm{yr})",
               limits = (extrema(pn.a_vec), (pn.b_vec[1], pn.b_vec[end]) ./ V_lit ./ yr),
               xscale = log10)
    hidexdecorations!(axT); hidespines!(axT)
    ax
end

Colorbar(fig[2, 3], hm, label = L"\Delta\chi^2 = \chi^2 - \chi^2_{\min}",
         ticks = ([0; DCHI2], ["0", "1σ", "2σ", "3σ"]))
Legend(fig[3, 1:2],
       vcat([LineElement(color = c, linewidth = 3) for c in conf_colors],
            [MarkerElement(color = :magenta, marker = :xcross, markersize = 26)]),
       ["68.3% (1σ)", "95.4% (2σ)", "99.7% (3σ)", "best fit"],
       orientation = :horizontal, framevisible = true, labelsize = 32,
       backgroundcolor = (:gray60, 1.0))

for (ax, lab) in zip([ax_a, ax_b, ax_land...], ["(a)", "(b)", "(c)", "(d)"])
    text!(ax, 0.02, 0.98, text = lab, space = :relative, align = (:left, :top), fontsize = 44)
end

save("fig/fig_krafla_v_constrained_combined.png", fig, px_per_unit = 2)
println("fig/fig_krafla_v_constrained_combined.png saved")
