using CairoMakie, CSV, DataFrames, QuadGK, SpecialFunctions
using LaTeXStrings, Printf, Optim, QuasiMonteCarlo, ColorSchemes, Statistics

include("../models/viscous.jl")
include("../models/maxwell.jl")
include("../parameters.jl")   # for H, E_L, nu_L (shared plate/half-space modulus assumption)

sci(x) = @sprintf("%.2e", x)   # shared by src/figures/krafla_fits_css.jl and krafla_hofton_diffusion.jl

include("../figures/krafla_fits_css.jl")
include("../figures/krafla_chi2_landscape_css.jl")
include("../figures/krafla_hofton_diffusion.jl")

# ── Data ──────────────────────────────────────────────────────────────────────
df       = CSV.read("data/extracted_data.csv", DataFrame)
obs_x_m  = df.X .* 1000.0
obs_u    = df.Y  ./ 100.0
obs_sig  = df.SD ./ 100.0
obs_x_km = obs_x_m ./ 1e3
N        = length(obs_x_m)

yr    = 3600.0 * 24 * 365
tau_8 = 8.0  * yr
tau   = 11.0 * yr

x_km = range(-500, stop=500, length=6000)
x_m  = x_km .* 1e3

function weighted_chi2(preds, obs, sig)
    sum(((obs .- preds) ./ sig) .^ 2)
end
aic(c, k) = c + 2k
bic(c, k) = c + k * log(N)

# ═══════════════════ Linear (2p/3p) — was krafla_fit.jl ════════════════════════
# NOTE: preserved as-is — these are hardcoded best-fit values transcribed from
# reports/relaxation_report.txt (2026-04-21), not fit inline here. Only the
# LOOCV loop below runs an actual optimizer. See Global Constraints.
k_opt         = 6.1858e-01
b_foulger_opt = 0.8099
a_step_opt    = 5.1257e-05
b_step_opt    = 1.1437
a_3p          = 6.6107e-05
b_3p          = 1.1495
tauR_3p       = 1.4841 * yr

chi2_foulger = weighted_chi2(diff_foulger(obs_x_m, tau_8, tau, b_foulger_opt, k_opt), obs_u, obs_sig)
chi2_step    = weighted_chi2(diff_step(obs_x_m, tau_8, tau, b_step_opt, a_step_opt), obs_u, obs_sig)
chi2_maxwell = weighted_chi2(diff_viscoelastic(obs_x_m, tau_8, tau, b_3p, a_3p, tauR_3p), obs_u, obs_sig)

println("\n=== Model comparison (N = $N) ===")
println("─────────────────────────────────────────────────────")
@printf("  %-14s  k   chi²    chi²_r   AIC     BIC\n", "Model")
println("─────────────────────────────────────────────────────")
for (name, c, k) in [("Foulger 2p", chi2_foulger, 2), ("Step 2p", chi2_step, 2), ("Maxwell 3p", chi2_maxwell, 3)]
    nu2 = N - k
    @printf("  %-14s  %d   %6.3f  %6.3f    %6.3f  %6.3f\n", name, k, c, c/nu2, aic(c,k), bic(c,k))
end
println("─────────────────────────────────────────────────────")
println("  ΔAIC (Step 2p − Foulger 2p) = $(round(aic(chi2_step,2) - aic(chi2_foulger,2), digits=2))")
println("  ΔBIC (Step 2p − Foulger 2p) = $(round(bic(chi2_step,2) - bic(chi2_foulger,2), digits=2))")

function loocv(predict_fn, p0; n_iter=3000)
    mspe = 0.0
    for i in 1:N
        idx   = filter(j -> j != i, 1:N)
        x_tr  = obs_x_m[idx];  u_tr = obs_u[idx];  s_tr = obs_sig[idx]
        cost(p) = any(p .<= 0) ? 1e10 : sum(((u_tr .- predict_fn(x_tr, p)) ./ s_tr) .^ 2)
        p_hat = Optim.minimizer(optimize(cost, p0, NelderMead(), Optim.Options(iterations=n_iter)))
        pred  = predict_fn([obs_x_m[i]], p_hat)[1]
        mspe += ((obs_u[i] - pred) / obs_sig[i])^2
    end
    return mspe / N
end

predict_foulger(x, p) = diff_foulger(x, tau_8, tau, p[2], p[1])
predict_step(x, p)    = diff_step(x, tau_8, tau, p[2], p[1])
predict_maxwell(x, p) = diff_viscoelastic(x, tau_8, tau, p[2], p[1], p[3])

println("\nRunning LOOCV ($(N) folds × 3 models) …")
mspe_foulger = loocv(predict_foulger, [k_opt,     b_foulger_opt])
mspe_step    = loocv(predict_step,    [a_step_opt, b_step_opt])
mspe_maxwell = loocv(predict_maxwell, [a_3p, b_3p, tauR_3p])

println("\n=== LOOCV results ===")
for (name, mspe) in [("Foulger 2p", mspe_foulger), ("Step 2p", mspe_step), ("Maxwell 3p", mspe_maxwell)]
    @printf("  %-14s  %.4f\n", name, mspe)
end
println("  ΔMSPE (Step 2p − Foulger 2p) = $(round(mspe_step - mspe_foulger, digits=4))")

# ═══════════════════════ CSS fits — was krafla_fit_css.jl ══════════════════════

pred_css(x, a, b, T) = diff_css(x, tau_8, tau, b, b / T, a)

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
    weighted_chi2(pred_css(obs_x_m, a, b, T), obs_u, obs_sig)
end
lb_free = [1e-7, 0.1, 100.0 * yr]
ub_free = [1e-2, 3.0, 5000.0 * yr]

println("Fitting CSS (T free) to Krafla data …")
p_css_free, chi2_css_free_val = fit_np(chi2_css_free, lb_free, ub_free)
a_css_free, b_css_free, T_css_free = p_css_free

T_fixed = 200.0 * yr
chi2_css_fixed = p -> begin
    a, b = p
    weighted_chi2(pred_css(obs_x_m, a, b, T_fixed), obs_u, obs_sig)
end
lb_fixed = [1e-7, 0.1]
ub_fixed = [1e-2, 3.0]

println("Fitting CSS (T = 200 yr fixed) to Krafla data …")
p_css_fixed, chi2_css_fixed_val = fit_np(chi2_css_fixed, lb_fixed, ub_fixed)
a_css_fixed, b_css_fixed = p_css_fixed

tau_R_from_a(a) = H / (a * (1 - nu_L))

pred_css_maxwell(x, a, b, T; k_limit_scale=1.0) = begin
    tau_R = tau_R_from_a(a)
    get_u_css_maxwell_closed.(x, tau,   b, b / T, a, tau_R; k_limit_scale=k_limit_scale) .-
    get_u_css_maxwell_closed.(x, tau_8, b, b / T, a, tau_R; k_limit_scale=k_limit_scale)
end

chi2_maxcss_free = p -> begin
    a, b, T = p
    (a <= 0 || b <= 0 || T <= 0) && return 1e12
    weighted_chi2(pred_css_maxwell(obs_x_m, a, b, T), obs_u, obs_sig)
end

println("Fitting Maxwell CSS (T free) to Krafla data …")
p_maxcss_free, chi2_maxcss_free_val = fit_np(chi2_maxcss_free, lb_free, ub_free; n=10000)
a_maxcss_free, b_maxcss_free, T_maxcss_free = p_maxcss_free
tauR_maxcss_free = tau_R_from_a(a_maxcss_free)

chi2_maxcss_fixed = p -> begin
    a, b = p
    (a <= 0 || b <= 0) && return 1e12
    weighted_chi2(pred_css_maxwell(obs_x_m, a, b, T_fixed), obs_u, obs_sig)
end

println("Fitting Maxwell CSS (T = 200 yr fixed) to Krafla data …")
p_maxcss_fixed, chi2_maxcss_fixed_val = fit_np(chi2_maxcss_fixed, lb_fixed, ub_fixed; n=10000)
a_maxcss_fixed, b_maxcss_fixed = p_maxcss_fixed
tauR_maxcss_fixed = tau_R_from_a(a_maxcss_fixed)

println("\n=== Model comparison (N = $N) ===")
for (name, c, k) in [("CSS 3p (T free)", chi2_css_free_val, 3), ("CSS 2p (T=200yr)", chi2_css_fixed_val, 2),
                      ("MaxCSS 3p (T free)", chi2_maxcss_free_val, 3), ("MaxCSS 2p (T=200yr)", chi2_maxcss_fixed_val, 2)]
    nu2 = N - k
    @printf("  %-20s  %d   %6.3f   %6.3f   %6.3f   %6.3f\n", name, k, c, c/nu2, aic(c,k), bic(c,k))
end

plot_krafla_fits_css(; obs_x_km, obs_u, obs_sig, x_km, x_m,
    a_css_free, b_css_free, T_css_free, a_css_fixed, b_css_fixed, T_fixed,
    a_maxcss_free, b_maxcss_free, T_maxcss_free, tauR_maxcss_free,
    a_maxcss_fixed, b_maxcss_fixed, tauR_maxcss_fixed, tau_8, tau, yr)

# ═══════════════════════ Remaining figures ══════════════════════════════════

plot_krafla_chi2_landscape_css(; obs_x_m, obs_u, obs_sig, tau_8, tau, a_css_fixed, b_css_fixed, T_fixed)

plot_krafla_hofton_diffusion()
