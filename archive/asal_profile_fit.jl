using CairoMakie, CSV, DataFrames, QuadGK, SpecialFunctions
using LaTeXStrings, Printf, Optim, Dates, Statistics

include("src/models.jl")

sci(x) = @sprintf("%.2e", x)

# ── Constants ──────────────────────────────────────────────────────────────────
yr      = 3600.0 * 24 * 365.25
t0_year = 1978.85

# Profile epoch: Vigny uses "mean date [1997] of GPS campaigns, 19 years from 1978"
t_epoch = 19.0 * yr   # s post-rifting

# Model center: x = +1 km NE of topographic rift axis (Petit Rift, Vigny §3.2)
x_center_m = 1000.0

# ── Load profile data (Vigny et al. 2007, Fig. 5, digitized) ─────────────────
df_prof  = CSV.read("data/asal_profile.csv", DataFrame)
obs_x_km = df_prof.X     # distance from topographic rift axis (km)
obs_v    = df_prof.Y     # rift-perp velocity, Somalia frame (mm/yr)
obs_sig  = df_prof.SD    # 1σ (mm/yr)
obs_x_m  = obs_x_km .* 1e3
N        = length(obs_x_km)

# ── Load time series (same as asal_fit.jl) ────────────────────────────────────
df_ts    = CSV.read("data/asal_timeseries.csv", DataFrame)
obs_year = df_ts.year
obs_dL   = df_ts.dL_mm
obs_sig_ts = df_ts.sig_mm
obs_t    = (obs_year .- t0_year) .* yr
N_ts     = length(obs_year)
x0_m     = 2000.0   # half EP00–DF00 perp. baseline (Vigny 2007 convention)

# ── Time-series forward models ΔL(t) = 2×u(x₀,t) [m→mm] ─────────────────────
predict_foulger_ts(t_s, p) = 2 .* get_u_foulger.(x0_m, t_s, p[1], p[2]) .* 1000
predict_step_ts(t_s, p)    = 2 .* get_u_step.(x0_m, t_s, p[2], p[1]) .* 1000
predict_maxwell_ts(t_s, p) = [2 * get_u_viscoelastic(x0_m, t, p[2], p[1], p[3]) * 1000
                               for t in t_s]

chi2_ts(pred, obs, sig) = sum(((obs .- pred) ./ sig) .^ 2)

function fit_ts(predict_fn, p0_list; n_iter=15000)
    best_p, best_c = p0_list[1], Inf
    for p0 in p0_list
        cost(p) = any(p .<= 0) ? 1e12 :
                  chi2_ts(predict_fn(obs_t, p), obs_dL, obs_sig_ts)
        res = optimize(cost, p0, NelderMead(), Optim.Options(iterations=n_iter))
        if Optim.minimum(res) < best_c
            best_p = Optim.minimizer(res)
            best_c = Optim.minimum(res)
        end
    end
    return best_p, best_c
end

p0_foulger_ts = [[κ, U] for κ in [0.001, 0.015, 0.1, 1.0]
                          for U in [0.1,   0.3,   0.5, 0.8]]
p0_step_ts    = [[a, ω] for a in [1e-6, 1e-5, 5e-5, 1e-4]
                         for ω in [0.2,  0.4,  0.6,  0.8]]
p0_maxwell_ts = [[a, ω, τR * yr] for a  in [1e-6, 1e-5, 5e-5]
                                   for ω  in [0.3,  0.5]
                                   for τR in [0.5,  2.0, 10.0]]

println("Fitting time series to get profile prediction parameters …")
p_foulger, chi2_foulger = fit_ts(predict_foulger_ts, p0_foulger_ts)
println("  Foulger 2p  χ² = $(round(chi2_foulger, digits=3))")
p_step,    chi2_step    = fit_ts(predict_step_ts,    p0_step_ts)
println("  Step    2p  χ² = $(round(chi2_step,    digits=3))")
p_maxwell, chi2_maxwell = fit_ts(predict_maxwell_ts, p0_maxwell_ts)
println("  Maxwell 3p  χ² = $(round(chi2_maxwell, digits=3))")

# ── Profile prediction: v(x, t_epoch) from time-series parameters ─────────────
# Models return m/s → convert to mm/yr (* 1000 * yr).
# x_rel = x_obs − x_center_m  (distance from Petit Rift center)
function pred_v_foulger(x_obs_m)
    x_rel = x_obs_m .- x_center_m
    get_v_foulger.(x_rel, t_epoch, p_foulger[1], p_foulger[2]) .* 1000 .* yr
end

function pred_v_step(x_obs_m)
    x_rel = x_obs_m .- x_center_m
    get_v_step.(x_rel, t_epoch, p_step[2], p_step[1]) .* 1000 .* yr
end

function pred_v_maxwell(x_obs_m)
    x_rel = x_obs_m .- x_center_m
    [get_v_relax(x, t_epoch, p_maxwell[2], p_maxwell[1], p_maxwell[3]) * 1000 * yr
     for x in x_rel]
end

# Vigny reference: κ = 0.015 m²/s, U₀ = 0.4 m (their Figure 3 fit)
pred_v_vigny(x_obs_m) =
    get_v_foulger.(x_obs_m .- x_center_m, t_epoch, 0.015, 0.4) .* 1000 .* yr

# ── Background plate rate: uniform offset ─────────────────────────────────────
# v_total(x) = v_transient(x) + V_plate
# V_plate is the Somalia–Danakil background rate; transient params fixed from
# time-series fit, only V_plate is fitted to the profile data (1 free param).
function fit_vplate(pred_v_fn)
    cost(p) = sum(((obs_v .- (pred_v_fn(obs_x_m) .+ p[1])) ./ obs_sig).^2)
    res = optimize(cost, [8.0], NelderMead(), Optim.Options(iterations=5000))
    Optim.minimizer(res)[1]
end

println("\nFitting Danakil background rate (step at x = 0) …")
Vp_vigny   = fit_vplate(pred_v_vigny)
Vp_foulger = fit_vplate(pred_v_foulger)
Vp_step    = fit_vplate(pred_v_step)
Vp_maxwell = fit_vplate(pred_v_maxwell)
for (name, Vp) in [("Vigny ref",  Vp_vigny),   ("Foulger 2p", Vp_foulger),
                    ("Step    2p", Vp_step),     ("Maxwell 3p", Vp_maxwell)]
    @printf("  %-12s  V_plate = %.2f mm/yr\n", name, Vp)
end

# ── Chi² helpers ─────────────────────────────────────────────────────────────
chi2_prof(pred_fn) = sum(((obs_v .- pred_fn(obs_x_m)) ./ obs_sig) .^ 2)
chi2_r_prof(pred_fn, k) = chi2_prof(pred_fn) / (N - k)

# k_eff = 1: only V_plate is fitted to profile; transient params from independent data
chi2_bg(pred_fn, Vp) = sum(((obs_v .- (pred_fn(obs_x_m) .+ Vp)) ./ obs_sig).^2)
chi2_r_bg(pred_fn, Vp) = chi2_bg(pred_fn, Vp) / (N - 1)

println("\nProfile chi² without background correction:")
@printf("  Vigny ref   χ²=%.2f  χ²_r=%.2f\n",
    chi2_prof(pred_v_vigny),   chi2_r_prof(pred_v_vigny, 2))
@printf("  Foulger 2p  χ²=%.2f  χ²_r=%.2f\n",
    chi2_prof(pred_v_foulger), chi2_r_prof(pred_v_foulger, 2))
@printf("  Step    2p  χ²=%.2f  χ²_r=%.2f\n",
    chi2_prof(pred_v_step),    chi2_r_prof(pred_v_step, 2))
@printf("  Maxwell 3p  χ²=%.2f  χ²_r=%.2f\n",
    chi2_prof(pred_v_maxwell), chi2_r_prof(pred_v_maxwell, 3))

println("\nProfile chi² with Danakil background step (1 free param: V_plate):")
@printf("  Vigny ref   V_plate=%5.2f  χ²=%.2f  χ²_r=%.2f\n",
    Vp_vigny,   chi2_bg(pred_v_vigny,   Vp_vigny),   chi2_r_bg(pred_v_vigny,   Vp_vigny))
@printf("  Foulger 2p  V_plate=%5.2f  χ²=%.2f  χ²_r=%.2f\n",
    Vp_foulger, chi2_bg(pred_v_foulger, Vp_foulger), chi2_r_bg(pred_v_foulger, Vp_foulger))
@printf("  Step    2p  V_plate=%5.2f  χ²=%.2f  χ²_r=%.2f\n",
    Vp_step,    chi2_bg(pred_v_step,    Vp_step),    chi2_r_bg(pred_v_step,    Vp_step))
@printf("  Maxwell 3p  V_plate=%5.2f  χ²=%.2f  χ²_r=%.2f\n",
    Vp_maxwell, chi2_bg(pred_v_maxwell, Vp_maxwell), chi2_r_bg(pred_v_maxwell, Vp_maxwell))

# ── Write report ───────────────────────────────────────────────────────────────
mkpath("reports")
open("reports/asal_profile_report.txt", "w") do io
    println(io, "Asal rift — velocity profile comparison (Fig. 5 equivalent)")
    println(io, "Generated: $(Dates.now())")
    println(io, "Data: Vigny et al. (2007) Fig. 5, N=$N stations")
    println(io, "Approach: FORWARD PREDICTION from time-series parameters (not an independent fit)")
    println(io, "Background: uniform V_plate offset fitted to profile (1 free param)")
    println(io, "Model center: x = +1 km NE of topographic axis (Petit Rift, Vigny §3.2)")
    println(io, "Profile epoch: t₀ = 19 yr (mean GPS epoch 1997, Vigny §3.2)")
    println(io)

    println(io, "=== Time-series best-fit parameters (from Fig. 3 data) ===")
    @printf(io, "  Foulger 2p:  κ = %.3e m²/s   U₀ = %.4f m   (χ² = %.3f)\n",
        p_foulger[1], p_foulger[2], chi2_foulger)
    @printf(io, "  Step    2p:  a = %.3e m/s    ω  = %.4f m   (χ² = %.3f)\n",
        p_step[1],    p_step[2],    chi2_step)
    @printf(io, "  Maxwell 3p:  a = %.3e m/s    ω  = %.4f m   τ_R = %.2f yr  (χ² = %.3f)\n",
        p_maxwell[1], p_maxwell[2], p_maxwell[3]/yr, chi2_maxwell)
    @printf(io, "  Vigny (2007): κ = 0.015 m²/s   U₀ = 0.4 m\n")
    println(io)

    println(io, "=== Profile residual without background (N=$N) ===")
    println(io, "  (χ²_r >> 1 reflects unmodeled Danakil background motion)")
    @printf(io, "  %-12s  χ²=%7.2f  χ²_r=%.2f\n", "Vigny ref",
        chi2_prof(pred_v_vigny),   chi2_r_prof(pred_v_vigny, 2))
    @printf(io, "  %-12s  χ²=%7.2f  χ²_r=%.2f\n", "Foulger 2p",
        chi2_prof(pred_v_foulger), chi2_r_prof(pred_v_foulger, 2))
    @printf(io, "  %-12s  χ²=%7.2f  χ²_r=%.2f\n", "Step 2p",
        chi2_prof(pred_v_step),    chi2_r_prof(pred_v_step, 2))
    @printf(io, "  %-12s  χ²=%7.2f  χ²_r=%.2f\n", "Maxwell 3p",
        chi2_prof(pred_v_maxwell), chi2_r_prof(pred_v_maxwell, 3))
    println(io)

    println(io, "=== Profile residual with Danakil background step (k_eff=1) ===")
    println(io, "  V_plate: best-fit Somalia–Danakil rate from profile alone")
    println(io, "─────────────────────────────────────────────────────────────")
    @printf(io, "  %-12s  V_plate (mm/yr)   χ²      χ²_r\n", "Model")
    println(io, "─────────────────────────────────────────────────────────────")
    for (name, pred_fn, Vp) in [("Vigny ref",  pred_v_vigny,   Vp_vigny),
                                  ("Foulger 2p", pred_v_foulger, Vp_foulger),
                                  ("Step 2p",    pred_v_step,    Vp_step),
                                  ("Maxwell 3p", pred_v_maxwell, Vp_maxwell)]
        @printf(io, "  %-12s  %14.2f   %7.2f  %5.2f\n",
            name, Vp, chi2_bg(pred_fn, Vp), chi2_r_bg(pred_fn, Vp))
    end
    println(io, "─────────────────────────────────────────────────────────────")
    println(io, "  Reference Somalia–Danakil rate: ~16 mm/yr (Vigny et al. 2007)")
end
println("reports/asal_profile_report.txt  saved")

# ── Figure ─────────────────────────────────────────────────────────────────────
mkpath("fig")
x_grid_km = range(-12, stop=12, length=800)
x_grid_m  = x_grid_km .* 1e3

fig = Figure(size=(1200, 700), fontsize=26)
ax  = Axis(fig[1, 1],
    xlabel = "Distance from rift axis (km)",
    ylabel = "Rift-perpendicular velocity (mm yr⁻¹)\nrelative to Somalia",
    limits = (-12, 12, -2, 22))

# Rift zone shading
vspan!(ax, -3, 3, color=(:lightgray, 0.3))
vspan!(ax, -1, 1, color=(:lightgray, 0.4))
vlines!(ax, [0.0], color=(:gray, 0.4), linestyle=:dot, linewidth=1)
vlines!(ax, [1.0], color=(:black, 0.3), linestyle=:dash, linewidth=1)

# Data
scatter!(ax, obs_x_km, obs_v, color=:black, markersize=10,
    label="Vigny et al. (2007)")
errorbars!(ax, obs_x_km, obs_v, obs_sig, color=:black)

# Model curves: transient + background step
lines!(ax, x_grid_km, pred_v_vigny(x_grid_m)   .+ Vp_vigny,
    color=(:red, 0.4), linewidth=2, linestyle=:dash,
    label="Foulger (Vigny 2007) | V_plate=$(round(Vp_vigny,digits=1)) mm/yr")
lines!(ax, x_grid_km, pred_v_foulger(x_grid_m) .+ Vp_foulger,
    color=:red, linewidth=2.5,
    label="Foulger 2p | κ=$(sci(p_foulger[1])) m²/s  V_plate=$(round(Vp_foulger,digits=1)) mm/yr")
lines!(ax, x_grid_km, pred_v_step(x_grid_m)    .+ Vp_step,
    color=:forestgreen, linewidth=2.5,
    label="Step 2p    | a=$(sci(p_step[1])) m/s  V_plate=$(round(Vp_step,digits=1)) mm/yr")

axislegend(ax, position=:rt, framevisible=true, fontsize=17)

text!(ax, -11, 1.0,
    text="Shaded: area activated in 1978 (light)\nand Petit Rift (dark)\nStep at x = 0: Danakil background",
    fontsize=13, color=(:gray, 0.7))

save("fig/fig_asal_profile.png", fig)
println("fig/fig_asal_profile.png  saved")
