using CairoMakie, CSV, DataFrames, QuadGK, SpecialFunctions
using LaTeXStrings, Printf, Optim, QuasiMonteCarlo, Dates, ColorSchemes, Statistics

include("../models/viscous.jl")
include("../models/maxwell.jl")
include("../parameters.jl")   # for H, E_L, nu_L (shared plate/half-space modulus assumption)

sci(x) = @sprintf("%.2e", x)   # shared by src/figures/asal_timeseries_css.jl

include("../figures/asal_timeseries_css.jl")
include("../figures/asal_velocity_profile_css.jl")
include("../figures/asal_chi2_landscape_css.jl")

# ── Constants ──────────────────────────────────────────────────────────────────
yr      = 3600.0 * 24 * 365.25
t0_year = 1978.85               # rifting onset (Vigny et al. 2007)

# ── Data ──────────────────────────────────────────────────────────────────────
df       = CSV.read("data/asal_timeseries_merged.csv", DataFrame)
obs_year = df.year
obs_dL   = df.dL_mm
obs_sig  = df.sig_mm
obs_t    = (obs_year .- t0_year) .* yr
N        = length(obs_year)
x0_m     = 2000.0
t_ref_s  = (obs_year[1] - t0_year) * yr

weighted_chi2(pred, obs, sig) = sum(((obs .- pred) ./ sig) .^ 2)

# ═══════════════════════ Linear (non-CSS) fits — was asal_fit.jl ═══════════════

predict_foulger(t, p) = 2 .* get_u_foulger.(x0_m, t, p[1], p[2]) .* 1000
predict_step(t, p)    = 2 .* get_u_step.(x0_m, t, p[2], p[1])    .* 1000
predict_foulger_tref(t, p) = predict_foulger(t, p) .- predict_foulger([t_ref_s], p)[1]
predict_step_tref(t, p)    = predict_step(t, p)    .- predict_step([t_ref_s],    p)[1]

function fit_model(predict_fn, p0_list; n_iter=15000)
    best_p, best_c = p0_list[1], Inf
    for p0 in p0_list
        cost(p) = any(p .<= 0) ? 1e12 : weighted_chi2(predict_fn(obs_t, p), obs_dL, obs_sig)
        res = optimize(cost, p0, NelderMead(), Optim.Options(iterations=n_iter))
        if Optim.minimum(res) < best_c
            best_p, best_c = Optim.minimizer(res), Optim.minimum(res)
        end
    end
    return best_p, best_c
end

p0_foulger_list = [[κ, U0] for κ  in [5e-3, 1e-2, 1.5e-2, 2e-2] for U0 in [0.2, 0.4, 0.6, 0.8]]
p0_step_list    = [[a, ω]  for a  in [1e-6, 1e-5, 5e-5, 1e-4]  for ω  in [0.2, 0.4, 0.6, 0.8]]

println("Fitting Asal time series (N = $N observations) …")
p_foulger, chi2_foulger = fit_model(predict_foulger_tref, p0_foulger_list)
println("  Foulger (fitted)  done  χ² = $(round(chi2_foulger; digits=3))")
p_step, chi2_step       = fit_model(predict_step_tref,    p0_step_list)
println("  Viscous           done  χ² = $(round(chi2_step;    digits=3))")

OFFSET_MM = 42.77
p_foulger_paper = [0.015, 0.4]
predict_foulger_fixoff(t, p) = predict_foulger(t, p) .- OFFSET_MM
chi2_foulger_paper = weighted_chi2(predict_foulger_fixoff(obs_t, p_foulger_paper), obs_dL, obs_sig)
println("  Foulger (paper)   done  χ² = $(round(chi2_foulger_paper; digits=3))")

aic(c, k) = c + 2k

function loocv(predict_fn, p0; n_iter=5000)
    mspe = 0.0
    for i in 1:N
        idx  = filter(j -> j != i, 1:N)
        t_tr = obs_t[idx];  u_tr = obs_dL[idx];  s_tr = obs_sig[idx]
        cost(p) = any(p .<= 0) ? 1e12 : sum(((u_tr .- predict_fn(t_tr, p)) ./ s_tr) .^ 2)
        p_hat = Optim.minimizer(optimize(cost, p0, NelderMead(), Optim.Options(iterations=n_iter)))
        pred  = predict_fn([obs_t[i]], p_hat)[1]
        mspe += ((obs_dL[i] - pred) / obs_sig[i])^2
    end
    return mspe / N
end

println("Running LOOCV ($N folds × 2 models) …")
mspe_foulger = loocv(predict_foulger_tref, p_foulger)
println("  Foulger (fitted)  done  MSPE = $(round(mspe_foulger; digits=4))")
mspe_step    = loocv(predict_step_tref,    p_step)
println("  Viscous           done  MSPE = $(round(mspe_step;    digits=4))")
mspe_foulger_paper = sum(((obs_dL .- predict_foulger_fixoff(obs_t, p_foulger_paper)) ./ obs_sig).^2) / N
println("  Foulger (paper)   done  MSPE = $(round(mspe_foulger_paper; digits=4))  [forward only, no refit]")

H_m = 10e3;  E_Pa = 40e9;  nu = 0.25
eta = E_Pa * H_m / (2 * p_step[1] * (1 - nu^2))
off_foulger       = predict_foulger([t_ref_s], p_foulger)[1]
off_foulger_paper = OFFSET_MM
off_step          = predict_step([t_ref_s],    p_step)[1]

mkpath("reports")
open("reports/asal_report.txt", "w") do io
    println(io, "Asal rift — Foulger vs Viscous (tref-off)")
    println(io, "Generated: $(Dates.now())")
    println(io, "Data: Vigny et al. (2007), EP00–DF00 time series, N=$N")
    println(io, "x₀ = $(x0_m) m  |  t_ref = $(round((obs_year[1]-t0_year); digits=3)) yr post-rifting")
    println(io, "Offset: model(t_ref) subtracted — no free parameter")
    println(io)
    println(io, "=== Parameters ===")
    @printf(io, "  Foulger (fitted):  κ = %.3e m²/s   U₀ = %.4f m   (off = %.2f mm)\n",
            p_foulger[1], p_foulger[2], off_foulger)
    @printf(io, "  Foulger (paper):   κ = %.3e m²/s   U₀ = %.4f m   (off = %.2f mm)  [fixed]\n",
            p_foulger_paper[1], p_foulger_paper[2], off_foulger_paper)
    @printf(io, "  Viscous:           a = %.3e m/s    ω  = %.4f m   (off = %.2f mm)\n",
            p_step[1], p_step[2], off_step)
    println(io)
    sep = "──────────────────────────────────────────────────────────────"
    println(io, "=== Model comparison (N = $N) ===")
    println(io, sep)
    @printf(io, "  %-20s  k   χ²       χ²_r    AIC      MSPE\n", "Model")
    println(io, sep)
    @printf(io, "  %-20s  2   %7.3f   %6.3f   %7.3f   %.4f\n",
            "Foulger (fitted)", chi2_foulger, chi2_foulger/(N-2), aic(chi2_foulger,2), mspe_foulger)
    @printf(io, "  %-20s  —   %7.3f   %6.3f   %s      %.4f  [forward]\n",
            "Foulger (paper)", chi2_foulger_paper, chi2_foulger_paper/(N-2), "  —    ", mspe_foulger_paper)
    @printf(io, "  %-20s  2   %7.3f   %6.3f   %7.3f   %.4f\n",
            "Viscous", chi2_step, chi2_step/(N-2), aic(chi2_step,2), mspe_step)
    println(io, sep)
    @printf(io, "  ΔAIC Viscous−Foulger(fitted) = %+.2f\n", aic(chi2_step,2) - aic(chi2_foulger,2))
    println(io, "  (|Δ|<2: indistinguishable; |Δ|>10: decisive)")
    println(io)
    println(io, "=== Physical checks ===")
    @printf(io, "  Foulger (fitted):  2U₀ = %.3f m   off = %.1f mm   (field: 1.9 m)\n", 2p_foulger[2], off_foulger)
    @printf(io, "  Foulger (paper):   2U₀ = %.3f m   off = %.1f mm   (field: 1.9 m)\n", 2p_foulger_paper[2], off_foulger_paper)
    @printf(io, "  Viscous:           2ω  = %.3f m   off = %.1f mm   (field: 1.9 m)\n", 2p_step[2], off_step)
    @printf(io, "  Viscous:           η   = %.2e Pa·s  (H=10 km, E=40 GPa, ν=0.25)\n", eta)
end
println("reports/asal_report.txt  saved")

# ═══════════════════════ CSS fits — was asal_fit_css.jl ════════════════════════

predict_css(t, p)      = 2 .* get_u_css.(x0_m, t, p[2], p[2] / p[3], p[1]) .* 1000
predict_css_tref(t, p) = predict_css(t, p) .- predict_css([t_ref_s], p)[1]

function lhs_log(lb, ub, n)
    plan = QuasiMonteCarlo.sample(n, log10.(lb), log10.(ub), LatinHypercubeSample())
    [10 .^ plan[:, i] for i in 1:n]
end

function fit_np(f, lb, ub; n = 8000)
    samps = lhs_log(lb, ub, n)
    p0    = samps[argmin([f(s) for s in samps])]
    res   = optimize(f, lb, ub, p0, Fminbox(NelderMead()))
    Optim.minimizer(res), Optim.minimum(res)
end

chi2_css_free = p -> begin
    a, b, T = p
    (a <= 0 || b <= 0 || T <= 0) && return 1e12
    weighted_chi2(predict_css_tref(obs_t, [a, b, T]), obs_dL, obs_sig)
end
lb_free = [1e-7, 0.1, 100.0 * yr]
ub_free = [1e-2, 3.0, 5000.0 * yr]

println("Fitting CSS (T free) to Asal time series …")
p_css_free, chi2_css_free_val = fit_np(chi2_css_free, lb_free, ub_free)
a_css_free, b_css_free, T_css_free = p_css_free
println("  CSS (T free)   done  a=$(sci(a_css_free)) m/s  b=$(round(b_css_free,digits=3)) m  T=$(round(T_css_free/yr,digits=1)) yr  χ²=$(round(chi2_css_free_val,digits=3))")

T_fixed = 500.0 * yr
chi2_css_fixed = p -> begin
    a, b = p
    (a <= 0 || b <= 0) && return 1e12
    weighted_chi2(predict_css_tref(obs_t, [a, b, T_fixed]), obs_dL, obs_sig)
end
lb_fixed = [1e-7, 0.1]
ub_fixed = [1e-2, 3.0]

println("Fitting CSS (T = 500 yr fixed) to Asal time series …")
p_css_fixed, chi2_css_fixed_val = fit_np(chi2_css_fixed, lb_fixed, ub_fixed)
a_css_fixed, b_css_fixed = p_css_fixed
println("  CSS (T=500yr)  done  a=$(sci(a_css_fixed)) m/s  b=$(round(b_css_fixed,digits=3)) m  χ²=$(round(chi2_css_fixed_val,digits=3))")

tau_R_from_a(a) = H / (a * (1 - nu_L))

predict_css_maxwell(t, p; k_limit_scale=1.0) = begin
    a, b, T = p
    tau_R = tau_R_from_a(a)
    2 .* get_u_css_maxwell_closed.(x0_m, t, b, b / T, a, tau_R; k_limit_scale=k_limit_scale) .* 1000
end
predict_css_maxwell_tref(t, p; k_limit_scale=1.0) =
    predict_css_maxwell(t, p; k_limit_scale=k_limit_scale) .-
    predict_css_maxwell([t_ref_s], p; k_limit_scale=k_limit_scale)[1]

chi2_maxcss_free = p -> begin
    a, b, T = p
    (a <= 0 || b <= 0 || T <= 0) && return 1e12
    weighted_chi2(predict_css_maxwell_tref(obs_t, [a, b, T]), obs_dL, obs_sig)
end

println("Fitting Maxwell CSS (T free) to Asal time series …")
p_maxcss_free, chi2_maxcss_free_val = fit_np(chi2_maxcss_free, lb_free, ub_free; n=10000)
a_maxcss_free, b_maxcss_free, T_maxcss_free = p_maxcss_free
tauR_maxcss_free = tau_R_from_a(a_maxcss_free)
println("  MaxCSS (T free)   done  a=$(sci(a_maxcss_free)) m/s  b=$(round(b_maxcss_free,digits=3)) m  T=$(round(T_maxcss_free/yr,digits=1)) yr  τ_R=$(round(tauR_maxcss_free/yr,digits=2)) yr (derived)  χ²=$(round(chi2_maxcss_free_val,digits=3))")

chi2_maxcss_fixed = p -> begin
    a, b = p
    (a <= 0 || b <= 0) && return 1e12
    weighted_chi2(predict_css_maxwell_tref(obs_t, [a, b, T_fixed]), obs_dL, obs_sig)
end

println("Fitting Maxwell CSS (T = 500 yr fixed) to Asal time series …")
p_maxcss_fixed, chi2_maxcss_fixed_val = fit_np(chi2_maxcss_fixed, lb_fixed, ub_fixed; n=10000)
a_maxcss_fixed, b_maxcss_fixed = p_maxcss_fixed
tauR_maxcss_fixed = tau_R_from_a(a_maxcss_fixed)
println("  MaxCSS (T=500yr)  done  a=$(sci(a_maxcss_fixed)) m/s  b=$(round(b_maxcss_fixed,digits=3)) m  τ_R=$(round(tauR_maxcss_fixed/yr,digits=2)) yr (derived)  χ²=$(round(chi2_maxcss_fixed_val,digits=3))")

bic(c, k) = c + k * log(N)
model_rows = [("CSS 3p (T free)",     chi2_css_free_val,     3),
              ("CSS 2p (T=500yr)",    chi2_css_fixed_val,    2),
              ("MaxCSS 3p (T free)",  chi2_maxcss_free_val,  3),
              ("MaxCSS 2p (T=500yr)", chi2_maxcss_fixed_val, 2)]

println("\n=== Model comparison (N = $N) ===")
println("─────────────────────────────────────────────────────────────")
@printf("  %-20s  k   chi²     chi²_r   AIC      BIC\n", "Model")
println("─────────────────────────────────────────────────────────────")
for (name, c, k) in model_rows
    nu2 = N - k
    @printf("  %-20s  %d   %6.3f   %6.3f   %6.3f   %6.3f\n", name, k, c, c/nu2, aic(c,k), bic(c,k))
end
println("─────────────────────────────────────────────────────────────")
println("  ΔAIC (MaxCSS T-free − CSS T-free)   = $(round(aic(chi2_maxcss_free_val,3) - aic(chi2_css_free_val,3), digits=2))  (same k=3)")
println("  ΔAIC (MaxCSS T=500yr − CSS T=500yr) = $(round(aic(chi2_maxcss_fixed_val,2) - aic(chi2_css_fixed_val,2), digits=2))  (same k=2)")
println("  (negative = Maxwell variant preferred at equal parameter count; |Δ|>2 positive evidence, |Δ|>6 strong)")

open("reports/asal_css_report.txt", "w") do io
    println(io, "Asal rift — CSS viscous vs Maxwell comparison (exploratory)")
    println(io, "Generated: $(Dates.now())")
    println(io, "Data: Vigny et al. (2007), EP00–DF00 time series, N=$N")
    println(io, "MaxCSS tau_R derived from a via shared plate/half-space modulus (H=10km, E=40GPa, nu=0.25)")
    println(io)
    println(io, "=== Parameters ===")
    @printf(io, "  CSS (T free):      a = %.3e m/s    b  = %.4f m   T = %.1f yr\n", a_css_free, b_css_free, T_css_free/yr)
    @printf(io, "  CSS (T=500yr):     a = %.3e m/s    b  = %.4f m   [T fixed]\n", a_css_fixed, b_css_fixed)
    @printf(io, "  MaxCSS (T free):   a = %.3e m/s    b  = %.4f m   T = %.1f yr   τ_R = %.2f yr (derived)\n",
            a_maxcss_free, b_maxcss_free, T_maxcss_free/yr, tauR_maxcss_free/yr)
    @printf(io, "  MaxCSS (T=500yr):  a = %.3e m/s    b  = %.4f m   τ_R = %.2f yr (derived)  [T fixed]\n",
            a_maxcss_fixed, b_maxcss_fixed, tauR_maxcss_fixed/yr)
    println(io)
    sep = "──────────────────────────────────────────────────────────────"
    println(io, "=== Model comparison (N = $N) ===")
    println(io, sep)
    @printf(io, "  %-20s  k   χ²       χ²_r    AIC      BIC\n", "Model")
    println(io, sep)
    for (name, c, k) in model_rows
        nu2 = N - k
        @printf(io, "  %-20s  %d   %7.3f   %6.3f   %7.3f   %7.3f\n", name, k, c, c/nu2, aic(c,k), bic(c,k))
    end
    println(io, sep)
end
println("reports/asal_css_report.txt  saved")

# ═══════════════════════ Figures (paper_v2.tex figures only) ═══════════════════

plot_asal_timeseries_css(; obs_year, obs_dL, obs_sig, x0_m, t_ref_s, t0_year, yr,
    a_css_free, b_css_free, T_css_free, a_css_fixed, b_css_fixed, T_fixed,
    a_maxcss_free, b_maxcss_free, T_maxcss_free, tauR_maxcss_free,
    a_maxcss_fixed, b_maxcss_fixed, tauR_maxcss_fixed)

plot_asal_velocity_profile_css(; a_css_free, b_css_free, T_css_free,
    a_css_fixed, b_css_fixed, T_fixed, a_maxcss_free, b_maxcss_free, T_maxcss_free,
    tauR_maxcss_free, a_maxcss_fixed, b_maxcss_fixed, tauR_maxcss_fixed, yr)

plot_asal_chi2_landscape_css(; obs_year, obs_dL, obs_sig, x0_m, t_ref_s, t0_year, yr,
    a_css_fixed, b_css_fixed, T_fixed)
