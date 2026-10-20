using CairoMakie, CSV, DataFrames, Printf, Optim, QuasiMonteCarlo, Statistics, Random, Dates

mkpath("fig")
mkpath("reports")
Random.seed!(20260421)

include("asal_v_constrained_common.jl")   # data, V cases, predictions, joint_chi2

# ── Fitting helper (mirrors src/fits/asal.jl's fit_np) ────────────────────────
# Top-K LHS refinement with plain NelderMead; box bounds enforced by the
# 1e12 penalty inside joint_chi2 (see the Fminbox note in src/fits/asal.jl).
function lhs_log(lb, ub, n)
    plan = QuasiMonteCarlo.sample(n, log10.(lb), log10.(ub), LatinHypercubeSample())
    [10 .^ plan[:, i] for i in 1:n]
end

function fit_np(f, lb, ub; n = 8000, n_top = 5, tol = 1e-3)
    samps      = lhs_log(lb, ub, n)
    scores     = [f(s) for s in samps]
    order      = sortperm(scores)
    opts       = Optim.Options(iterations = 8000)
    candidates = [(p = Optim.minimizer(res), chi2 = Optim.minimum(res))
                  for res in (optimize(f, samps[i], NelderMead(), opts) for i in order[1:n_top])]
    sort!(candidates, by = cand -> cand.chi2)
    best    = candidates[1]
    n_agree = count(cand -> abs(cand.chi2 - best.chi2) < tol * best.chi2, candidates)
    println("    top-$n_top LHS refinements: chi2 = ", round.([cand.chi2 for cand in candidates], digits=3),
            "  ($n_agree/$n_top agree within $(tol*100)%)")
    n_agree < n_top && @warn "top-$n_top LHS refinements did not all agree -- possible local-minimum sensitivity" best_chi2=best.chi2
    best.p, best.chi2
end

# ── Fits: EV and EVE, for each entry of V_cases (currently V = 13 mm/yr) ─────
results = map(V_cases) do vc
    V    = vc.V
    lb   = [1e-7, V * T_min_no_repeat]
    ub   = [1e-2, 3.0]

    println("Fitting EV RSS (T = b/V, $(vc.label)) jointly to Asal time series + profile …")
    p_ev, chi2_ev = fit_np(p -> joint_chi2(pred_ts_ev, pred_prof_ev, p[1], p[2], V), lb, ub; n = 8000, n_top = 5)

    println("Fitting EVE RSS (T = b/V, $(vc.label)) jointly to Asal time series + profile …")
    p_eve, chi2_eve = fit_np(p -> joint_chi2(pred_ts_eve, pred_prof_eve, p[1], p[2], V), lb, ub; n = 10000, n_top = 3)

    (; vc.label, V,
       a_ev = p_ev[1],  b_ev = p_ev[2],  T_ev = p_ev[2] / V,  chi2_ev,
       a_eve = p_eve[1], b_eve = p_eve[2], T_eve = p_eve[2] / V, tauR_eve = tau_R_from_a(p_eve[1]), chi2_eve)
end

# ── Report ───────────────────────────────────────────────────────────────────
k  = 2                    # a, b free; T, V_pl (and tau_R for EVE) derived
nu = N_total - k
aic(c) = c + 2k
bic(c) = c + k * log(N_total)

function write_report(io)
    println(io, "Asal rift — V-constrained RSS fits (T = b/V, V_pl = V/2)")
    println(io, "Generated: $(Dates.now())")
    println(io, "Data: Vigny et al. (2007) EP00–DF00 time series (N=$N_ts) + velocity profile (N=$N_prof), jointly fit")
    println(io, "Free: a, b.  Derived: T = b/V, V_pl = V/2, tau_R = H/(a(1-nu)) for EVE (H=10km, nu=0.25)")
    @printf(io, "Floor: T >= %.1f yr (no repeat intrusion 1978-2016, Smittarello et al. 2016)\n", T_min_no_repeat / yr)
    println(io)
    sep = "──────────────────────────────────────────────────────────────────────────────"
    for r in results
        println(io, "=== $(r.label) ===")
        @printf(io, "  EV RSS:   a = %.3e m/s   b = %.4f m   T = %.1f yr   V_pl = %.1f mm/yr\n",
                r.a_ev, r.b_ev, r.T_ev / yr, r.V * yr * 1000 / 2)
        @printf(io, "  EVE RSS:  a = %.3e m/s   b = %.4f m   T = %.1f yr   τ_R = %.2f yr   V_pl = %.1f mm/yr\n",
                r.a_eve, r.b_eve, r.T_eve / yr, r.tauR_eve / yr, r.V * yr * 1000 / 2)
        for (name, T) in [("EV", r.T_ev), ("EVE", r.T_eve)]
            T <= T_min_no_repeat * 1.001 && println(io, "  NOTE: $name fit sits on the T >= T_min_no_repeat floor")
        end
        println(io, sep)
        @printf(io, "  %-10s  k   χ²        χ²_r     AIC       BIC\n", "Model")
        println(io, sep)
        for (name, c) in [("EV RSS", r.chi2_ev), ("EVE RSS", r.chi2_eve)]
            @printf(io, "  %-10s  %d   %8.3f  %7.3f  %8.3f  %8.3f\n", name, k, c, c / nu, aic(c), bic(c))
        end
        println(io, sep)
        @printf(io, "  ΔAIC (EVE − EV) = %+.2f\n\n", aic(r.chi2_eve) - aic(r.chi2_ev))
    end
end

write_report(stdout)
open(write_report, "reports/asal_css_v_constrained_report.txt", "w")
println("reports/asal_css_v_constrained_report.txt  saved")

# Best-fit parameters, read back by src/figures/asal_chi2_landscape_v_constrained.jl
CSV.write(params_csv, DataFrame(
    V_mm_yr = [round(Int, r.V * yr * 1000) for r in results for _ in 1:2],
    model   = [m for _ in results for m in ("EV", "EVE")],
    a       = [x for r in results for x in (r.a_ev, r.a_eve)],
    b       = [x for r in results for x in (r.b_ev, r.b_eve)],
    chi2    = [x for r in results for x in (r.chi2_ev, r.chi2_eve)]))
println("$params_csv  saved")

# ── Figure: time series (left) and velocity profile (right) ──────────────────
t_yr_grid = range(0.05, stop = 26.0, length = 600)
t_s_grid  = t_yr_grid .* yr
x_grid_km = collect(range(-12, 12, length = 800))
x_grid_m  = x_grid_km .* 1e3

styles = [(ev = (color = :purple, linestyle = :solid), eve = (color = :navy,  linestyle = :solid)),
          (ev = (color = :orange, linestyle = :dash),  eve = (color = :brown, linestyle = :dash))]

fig = Figure(size = (2600, 1000), fontsize = 32)
ax1 = Axis(fig[1, 1], xlabel = "Year", ylabel = "EP00–DF00 distance change (mm)",
    title = "Asal time series", limits = (t0_year - 1, 2005, -80, 580))
ax2 = Axis(fig[1, 2], xlabel = "Distance from rift axis (km)",
    ylabel = "Rift-perpendicular velocity, Somalia frame (mm yr⁻¹)",
    title = "Asal velocity profile, t = 19 yr", limits = (-12, 12, -2, 22))

for (r, st) in zip(results, styles)
    Vmm = round(Int, r.V * yr * 1000)
    lab_ev  = "EV RSS, V=$Vmm | a=$(sci(r.a_ev)) m/s  b=$(round(r.b_ev, digits=2)) m  T=$(round(Int, r.T_ev / yr)) yr"
    lab_eve = "EVE RSS, V=$Vmm | a=$(sci(r.a_eve)) m/s  b=$(round(r.b_eve, digits=2)) m  T=$(round(Int, r.T_eve / yr)) yr  τ_R=$(round(r.tauR_eve / yr, digits=2)) yr"

    lines!(ax1, t_yr_grid .+ t0_year, tref(t -> pred_ts_ev(t, r.a_ev, r.b_ev, r.V), t_s_grid);
        st.ev..., linewidth = 2.5, label = lab_ev)
    lines!(ax1, t_yr_grid .+ t0_year, tref(t -> pred_ts_eve(t, r.a_eve, r.b_eve, r.V; k_limit_scale=30.0), t_s_grid);
        st.eve..., linewidth = 2.5, label = lab_eve)

    lines!(ax2, x_grid_km, pred_prof_ev(x_grid_m, r.a_ev, r.b_ev, r.V); st.ev..., linewidth = 2.5)
    lines!(ax2, x_grid_km, pred_prof_eve(x_grid_m, r.a_eve, r.b_eve, r.V; k_limit_scale=30.0); st.eve..., linewidth = 2.5)
end

scatter!(ax1, obs_year, obs_dL, color = :black, markersize = 9, label = "Vigny et al. (2007)")
errorbars!(ax1, obs_year, obs_dL, obs_sig, color = :black)
vlines!(ax2, [0.0], color = (:gray, 0.4), linestyle = :dot, linewidth = 1.0)
vlines!(ax2, [x_center_m / 1000], color = (:black, 0.3), linestyle = :dash, linewidth = 1.0)
errorbars!(ax2, prof_x_km, prof_v, prof_sig, color = :black, whiskerwidth = 5)
scatter!(ax2, prof_x_km, prof_v, color = :black, markersize = 9)

axislegend(ax1, position = :rb, framevisible = true, labelsize = 18)

save("fig/fig_asal_css_v_constrained.png", fig, px_per_unit = 2)
println("fig/fig_asal_css_v_constrained.png saved")
