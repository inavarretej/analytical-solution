using CairoMakie, CSV, DataFrames, Printf, Optim, QuasiMonteCarlo, Statistics, Random, Dates, ColorSchemes, LaTeXStrings

# Asal: EV RSS and paper EVE RSS (best fits read from
# reports/asal_css_v_constrained_params.csv, not refit) against the
# alternative EVE RSS of src/models/maxwell_alt.jl (initial state from the
# visco-elastic kernel integrated across each intrusion; vault note "EVE
# initial condition - missing instantaneous elastic response"), fitted here
# with the same procedure, data, T = b/V and T floor as asal_css_v_constrained.jl.
# Comparison only: no existing script, figure or report is modified.
#   fig/fig_asal_eve_alt_compare.png
#   reports/asal_eve_alt_report.txt, reports/asal_eve_alt_params.csv
mkpath("fig"); mkpath("reports")
Random.seed!(20261020)

include("asal_v_constrained_common.jl")   # data, predictions, joint_chi2 (unchanged)
include("../models/maxwell_alt.jl")

V = V_cases[1].V

pred_ts_alt(t, a, b, V; k_limit_scale = nothing) =
    2 .* get_u_css_maxwell_alt.(x0_m, mod.(t, b / V), b, V, a, tau_R_from_a(a)) .* 1000
pred_prof_alt(x, a, b, V; k_limit_scale = nothing) =
    get_v_css_maxwell_alt.(x_rel(x), mod(t_epoch_profile, b / V), b, V, a, tau_R_from_a(a)) .* 1000 .* yr .+
    V * yr * 1000 / 2

# ── Existing EV and paper-EVE best fits ──────────────────────────────────────
P = CSV.read(params_csv, DataFrame)
row(m) = only(eachrow(filter(r -> r.model == m, P)))
ev, eve = row("EV"), row("EVE")

# ── Fit the alternative EVE (same helper as asal_css_v_constrained.jl) ───────
function lhs_log(lb, ub, n)
    plan = QuasiMonteCarlo.sample(n, log10.(lb), log10.(ub), LatinHypercubeSample())
    [10 .^ plan[:, i] for i in 1:n]
end
function fit_np(f, lb, ub; n = 4000, n_top = 3)
    samps  = lhs_log(lb, ub, n)
    scores = [f(s) for s in samps]
    order  = sortperm(scores)
    opts   = Optim.Options(iterations = 8000)
    cands  = [(p = Optim.minimizer(r), chi2 = Optim.minimum(r))
              for r in (optimize(f, samps[i], NelderMead(), opts) for i in order[1:n_top])]
    sort!(cands, by = c -> c.chi2)
    println("    top-$n_top refinements: chi2 = ", round.([c.chi2 for c in cands], digits = 3))
    cands[1].p, cands[1].chi2
end

println("Fitting alternative EVE RSS (T = b/V) to Asal time series + profile …")
lb, ub = [1e-7, V * T_min_no_repeat], [1e-2, 3.0]
t_fit = @elapsed p_alt, chi2_alt = fit_np(p -> joint_chi2(pred_ts_alt, pred_prof_alt, p[1], p[2], V), lb, ub)
a_alt, b_alt = p_alt
@printf("  a = %.4e m/s, b = %.4f m, χ² = %.3f  [%.0f s]\n", a_alt, b_alt, chi2_alt, t_fit)

models = [(name = "EV RSS",               a = ev.a,  b = ev.b,  ts = pred_ts_ev,  pr = pred_prof_ev,  color = :purple,     style = :solid),
          (name = "EVE RSS (paper)",      a = eve.a, b = eve.b, ts = pred_ts_eve, pr = pred_prof_eve, color = :navy,       style = :dash),
          (name = "EVE RSS (alternative)", a = a_alt, b = b_alt, ts = pred_ts_alt, pr = pred_prof_alt, color = :darkorange, style = :solid)]

split_chi2(m) = (weighted_chi2(tref(t -> m.ts(t, m.a, m.b, V), obs_t), obs_dL, obs_sig),
                 weighted_chi2(m.pr(prof_x_m, m.a, m.b, V), prof_v, prof_sig))
stats = map(models) do m
    cts, cpr = split_chi2(m)
    (; m.name, m.a, m.b, T = m.b / V, tauR = m.name == "EV RSS" ? NaN : tau_R_from_a(m.a),
       chi2 = cts + cpr, chi2_ts = cts, chi2_prof = cpr)
end

# ── Alternative-EVE χ² landscape in (a, b) ───────────────────────────────────
const DCHI2 = [2.30, 6.17, 11.83]
function chi2_grid(a_vec, b_vec)
    G = Matrix{Float64}(undef, length(a_vec), length(b_vec))
    Threads.@threads for j in eachindex(b_vec)
        for i in eachindex(a_vec)
            G[i, j] = b_vec[j] / V < T_min_no_repeat ? NaN :
                      joint_chi2_raw(pred_ts_alt, pred_prof_alt, a_vec[i], b_vec[j], V)
        end
    end
    G
end
function zoomed_grid(a0, b0; n_search = 31, n_final = 121, max_iter = 8)
    a_half, b_half = 0.05, 0.02 * b0
    grid(n) = (10 .^ range(log10(a0) - a_half, log10(a0) + a_half, length = n),
               collect(range(b0 - b_half, b0 + b_half, length = n)))
    for _ in 1:max_iter
        a_vec, b_vec = grid(n_search)
        G = chi2_grid(a_vec, b_vec)
        m3 = replace(G, NaN => Inf) .≤ minimum(filter(!isnan, G)) + DCHI2[3]
        (any(m3[1, :]) || any(m3[end, :]) || any(m3[:, 1]) || any(m3[:, end])) || break
        a_half *= 1.6; b_half *= 1.6
    end
    a_vec, b_vec = grid(n_final)
    (a_vec, b_vec, chi2_grid(a_vec, b_vec))
end
println("Building alternative-EVE landscape …")
t_land = @elapsed a_vec, b_vec, G = zoomed_grid(a_alt, b_alt)
cmin = min(minimum(filter(!isnan, G)), chi2_alt)
@printf("  done [%.0f s], grid χ²_min = %.3f\n", t_land, minimum(filter(!isnan, G)))
in1 = replace(G, NaN => Inf) .≤ cmin + DCHI2[1]
a1 = extrema(a_vec[[i for i in axes(G, 1) if any(in1[i, :])]])
b1 = extrema(b_vec[[j for j in axes(G, 2) if any(in1[:, j])]])

# ── Report ───────────────────────────────────────────────────────────────────
k = 2; nu = N_total - k; aic(c) = c + 2k
open("reports/asal_eve_alt_report.txt", "w") do io
    println(io, "Asal — EV RSS vs paper EVE RSS vs alternative EVE RSS (T = b/V, V = 13 mm/yr)")
    println(io, "Generated: $(Dates.now())")
    println(io, "EV and paper-EVE best fits read from $params_csv (not refit); alternative EVE fitted here.")
    println(io, "Alternative EVE: src/models/maxwell_alt.jl (initial state from the kernel integrated across each intrusion).")
    println(io, "Same data, T = b/V, V_pl = V/2, T floor and tau_R = H/(a(1-nu)) as asal_css_v_constrained.jl.\n")
    @printf(io, "  %-22s %10s %8s %7s %7s %9s %9s %9s %7s %8s\n",
            "Model", "a (m/s)", "b (m)", "T (yr)", "τR (yr)", "χ²", "χ²_ts", "χ²_prof", "χ²_r", "ΔAIC")
    for s in stats
        @printf(io, "  %-22s %10.3e %8.4f %7.1f %7.2f %9.2f %9.2f %9.2f %7.3f %+8.2f\n",
                s.name, s.a, s.b, s.T / yr, s.tauR / yr, s.chi2, s.chi2_ts, s.chi2_prof, s.chi2 / nu,
                aic(s.chi2) - aic(stats[1].chi2))
    end
    @printf(io, "\n  χ²/N: time series (N=%d) and profile (N=%d)\n", N_ts, N_prof)
    for s in stats
        @printf(io, "  %-22s  ts %.3f   profile %.3f\n", s.name, s.chi2_ts / N_ts, s.chi2_prof / N_prof)
    end
    @printf(io, "\n  Alternative EVE 1σ (Δχ² ≤ 2.30 projections): a ∈ [%.3e, %.3e] m/s, b ∈ [%.4f, %.4f] m, T ∈ [%.1f, %.1f] yr\n",
            a1..., b1..., (b1 ./ V ./ yr)...)
    @printf(io, "  η (H = 10 km, E = 40 GPa, ν = 0.25): EV %.2e, EVE paper %.2e, EVE alternative %.2e Pa s\n",
            [E_L * H / (2 * s.a * (1 - nu_L^2)) for s in stats]...)
end
CSV.write("reports/asal_eve_alt_params.csv",
          DataFrame(model = [s.name for s in stats], a = [s.a for s in stats], b = [s.b for s in stats],
                    T_yr = [s.T / yr for s in stats], tau_R_yr = [s.tauR / yr for s in stats],
                    chi2 = [s.chi2 for s in stats], chi2_ts = [s.chi2_ts for s in stats],
                    chi2_prof = [s.chi2_prof for s in stats]))
print(read("reports/asal_eve_alt_report.txt", String))

# ── Figure ───────────────────────────────────────────────────────────────────
fig = Figure(size = (2600, 2000), fontsize = 36)
t_grid  = collect(range(t_ref_s, (2004 - t0_year) * yr, length = 400))
x_grid  = collect(range(-12e3, 12e3, length = 801))
lab(m, s) = @sprintf("%s  (χ² = %.0f)", m.name, s.chi2)

ax_a = Axis(fig[1, 1], title = "EP00–DF00 baseline change", xlabel = "Year", ylabel = "ΔL (mm)",
            limits = (t0_year, 2005, -80, 580))
ax_b = Axis(fig[1, 2], title = "Velocity profile, τ = 19 yr (Somalia frame)", xlabel = "Distance from rift axis (km)",
            ylabel = "Velocity (mm/yr)", limits = (-12, 12, -2, 22))
ax_c = Axis(fig[2, 1], title = "Time-series residuals (data − model)", xlabel = "Year", ylabel = "Residual (mm)",
            limits = (t0_year, 2005, nothing, nothing))
for (m, s) in zip(models, stats)
    lines!(ax_a, t0_year .+ t_grid ./ yr, tref(t -> m.ts(t, m.a, m.b, V), t_grid), color = m.color,
           linestyle = m.style, linewidth = 4, label = lab(m, s))
    lines!(ax_b, x_grid ./ 1e3, m.pr(x_grid, m.a, m.b, V), color = m.color, linestyle = m.style, linewidth = 4)
    res = obs_dL .- tref(t -> m.ts(t, m.a, m.b, V), obs_t)
    scatterlines!(ax_c, obs_year, res, color = m.color, linestyle = m.style, markersize = 12, linewidth = 2)
end
errorbars!(ax_a, obs_year, obs_dL, obs_sig, color = :black)
scatter!(ax_a, obs_year, obs_dL, color = :black, markersize = 12)
errorbars!(ax_b, prof_x_km, prof_v, prof_sig, color = :black)
scatter!(ax_b, prof_x_km, prof_v, color = :black, markersize = 12)
vlines!(ax_b, [x_center_m / 1e3], color = :gray, linestyle = :dash)
errorbars!(ax_c, obs_year, zeros(N_ts), obs_sig, color = (:black, 0.5))
hlines!(ax_c, [0.0], color = :black)
axislegend(ax_a, position = :rb, labelsize = 30)

ax_d = Axis(fig[2, 2], xscale = log10, xtickformat = vs -> [@sprintf("%.2e", v) for v in vs],
            xlabel = L"a\ \ (\mathrm{m\,s^{-1}})", ylabel = L"b\ \ (\mathrm{m})",
            title = "EVE RSS (alternative), Δχ²")
cmap = cgrad(ColorSchemes.bamako100, 20, categorical = true)
hm = heatmap!(ax_d, a_vec, b_vec, G .- cmin, colormap = cmap, colorrange = (0, DCHI2[3] * 1.4),
              highclip = cmap[end], nan_color = :gray80)
for (lev, c) in zip(DCHI2, [:white, :gold, :tomato])
    contour!(ax_d, a_vec, b_vec, G .- cmin, levels = [lev], color = c, linewidth = 3)
end
scatter!(ax_d, [a_alt], [b_alt], color = :magenta, marker = :xcross, markersize = 26)
Colorbar(fig[2, 3], hm, label = L"\Delta\chi^2", ticks = ([0; DCHI2], ["0", "1σ", "2σ", "3σ"]))
for (ax, l) in zip((ax_a, ax_b, ax_c, ax_d), ("(a)", "(b)", "(c)", "(d)"))
    text!(ax, 0.02, 0.98, text = l, space = :relative, align = (:left, :top), fontsize = 44)
end
save("fig/fig_asal_eve_alt_compare.png", fig, px_per_unit = 2)
println("Saved fig/fig_asal_eve_alt_compare.png")
