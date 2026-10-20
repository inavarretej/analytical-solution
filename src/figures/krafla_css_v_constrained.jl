using CairoMakie, CSV, DataFrames, Printf, Optim, QuasiMonteCarlo, Statistics, Random

mkpath("fig")
Random.seed!(20260421)

include("../models/viscous.jl")
include("../models/maxwell.jl")
include("../parameters.jl")   # H, E_L, nu_L -- shared plate/half-space modulus assumption
sci(x) = @sprintf("%.2e", x)

# CSS/MaxCSS joint fit to both Krafla datasets with the recurrence interval T
# NOT free: T is set by the plate-rate constraint V = b/(2T), where b is the
# dike opening (the same convention used everywhere else -- get_u_step's
# discontinuity at the source is exactly b, eq:opening_conserved) and V is a
# literature far-field spreading rate (not fit).
#
# The model functions (get_u_css, get_u_css_maxwell_closed, ...) take a
# "v_plate" argument and set T = b/v_plate internally, with far-field velocity
# v_plate/2, not v_plate (see vault note "Factor-of-2 error in the EV RSS
# velocity formula" / "From STEP to RSS velocity (corrected derivation)").
#
# V_lit = 1.8 cm/yr is DeMets et al.'s (1994) NUVEL-1A relative Eurasia-North
# America plate velocity -- global plate-motion models report RELATIVE
# (full, two-plate) rates by convention, and our EV/EVE model's own geometry
# is symmetric spreading about a fixed rift axis (v_RSS is odd in x). So the
# physically correct one-sided far-field rate for a single flank is V_lit/2,
# and passing v_plate_arg = V_lit directly below already gives exactly that
# (far-field velocity = v_plate_arg/2 = V_lit/2) -- no extra factor needed.
# (An earlier version of this comment argued for v_plate_arg = 2*V_lit,
# reading Foulger et al. 1992's informal "expected displacement" arithmetic
# as evidence V_lit is used as a one-sided rate; on closer reading that
# passage is about the peak/near-field amplitude of the transient relative
# to background, not a rigorous statement about the half/full convention,
# so it isn't good evidence either way -- the standard NUVEL-1A relative-rate
# convention is the safer default.)
yr    = 3600.0 * 24 * 365
V_lit = 0.018 / yr   # m/s -- DeMets et al. (1994) measured NE Iceland spreading rate, 1.8 cm/yr, as used by Hofton & Foulger (1996)
v_plate_arg = V_lit   # passed to get_*_css* directly; far-field velocity = v_plate_arg/2 = V_lit/2 (one-sided)

tau_R_from_a(a) = H / (a * (1 - nu_L))

# Time origin (dike-injection date, t = 0) for BOTH datasets, set with the env
# variable KRAFLA_T0 (e.g. 1975 = start of the episode, 1979 = end of the dike
# intrusions, Foulger et al. 1992). Default 1975, the common origin that fits
# best (same default as krafla_v_constrained_1sigma.jl).
T0_env = get(ENV, "KRAFLA_T0", "")
T0_foulger = T0_hofton = isempty(T0_env) ? 1975 : parse(Int, T0_env)
fig_suffix = isempty(T0_env) ? "" : "_t0_$(T0_env)"

# ── Data ─────────────────────────────────────────────────────────────────────
foulger          = CSV.read("data/extracted_data.csv", DataFrame)
obs_x_m_foulger  = foulger.X .* 1000.0
obs_u_foulger    = foulger.Y  ./ 100.0
obs_sig_foulger  = foulger.SD ./ 100.0
N_foulger        = length(obs_x_m_foulger)
tau_i_foulger    = (1987 - T0_foulger) * yr
tau_f_foulger    = (1990 - T0_foulger) * yr

hofton           = CSV.read("data/hofton.csv", DataFrame)
obs_x_m_hofton   = hofton.X .* 1000.0
obs_u_hofton     = hofton.Y  ./ 100.0
obs_sig_hofton   = hofton.SD ./ 100.0
N_hofton         = length(obs_x_m_hofton)
tau_i_hofton     = (1987 - T0_hofton) * yr
tau_f_hofton     = (1992 - T0_hofton) * yr

N_total = N_foulger + N_hofton

# ── Fitting helpers (mirrors src/fits/krafla.jl's lhs_log / fit_np) ───────────
function lhs_log(lb, ub, n)
    plan = QuasiMonteCarlo.sample(n, log10.(lb), log10.(ub), LatinHypercubeSample())
    [10 .^ plan[:, i] for i in 1:n]
end

weighted_chi2(preds, obs, sig) = sum(((obs .- preds) ./ sig) .^ 2)

# Hofton offset c, profiled out analytically. Hofton & Foulger's (1996) zero
# level is arbitrary (their Fig. 5 caption), so the model is shifted by a
# constant c on that dataset only. chi2 is a parabola in c, and
# d(chi2)/dc = 0 gives c* = sum(w (u - m)) / sum(w), w = 1/sigma^2: the
# weighted mean residual. Every fit below therefore searches only its two
# model parameters; c is still a free parameter (nu = N - 3), it is just
# solved exactly at each trial point instead of searched for.
const w_hofton = 1 ./ obs_sig_hofton .^ 2
best_c(pred_h) = sum(w_hofton .* (obs_u_hofton .- pred_h)) / sum(w_hofton)

# Joint chi2 of a two-parameter model pred(x, tau_i, tau_f, p1, p2) over both
# datasets, with the Hofton offset at its optimum. Returns (chi2, c).
function joint_chi2(pred, p1, p2)
    (p1 <= 0 || p2 <= 0) && return (1e12, 0.0)
    pf = pred(obs_x_m_foulger, tau_i_foulger, tau_f_foulger, p1, p2)
    ph = pred(obs_x_m_hofton,  tau_i_hofton,  tau_f_hofton,  p1, p2)
    c  = best_c(ph)
    (weighted_chi2(pf, obs_u_foulger, obs_sig_foulger) +
     weighted_chi2(ph .+ c, obs_u_hofton, obs_sig_hofton), c)
end

# LHS over (p1, p2) (both positive, log-uniform), then Nelder-Mead from the
# top n_top LHS candidates, keeping the best. Warns if the refinements
# disagree by more than `tol` in chi2 (possible local-minimum sensitivity;
# the EVE/MaxCSS objective in particular has shallow local minima,
# consistent with the k_limit_scale=1 integration noted in maxwell.jl).
function fit_np(pred, lb, ub; n = 8000, n_top = 5, tol = 1e-3)
    f(q)    = joint_chi2(pred, q[1], q[2])[1]
    samps   = lhs_log(lb, ub, n)
    order   = sortperm([f(s) for s in samps])
    candidates = [(p = Optim.minimizer(res), chi2 = Optim.minimum(res))
                  for res in (optimize(f, lb, ub, samps[i], Fminbox(NelderMead())) for i in order[1:n_top])]
    sort!(candidates, by = cand -> cand.chi2)
    best    = candidates[1]
    n_agree = count(cand -> abs(cand.chi2 - best.chi2) < tol * best.chi2, candidates)
    println("    top-$n_top LHS refinements: chi2 = ", round.([cand.chi2 for cand in candidates], digits=3),
            "  ($n_agree/$n_top agree within $(tol*100)%)")
    n_agree < n_top && @warn "top-$n_top LHS refinements did not all agree -- possible local-minimum sensitivity" best_chi2=best.chi2
    best.p, best.chi2
end

# Grid polish after the LHS fit: the LHS + Nelder-Mead search can stop in a
# shallow local minimum (for EV at T0 = 1975 an earlier version returned
# chi2 = 57.3 where the true minimum is 54.7, see
# src/figures/krafla_v_constrained_1sigma.jl). Scans chi2 on a (p1, p2) grid
# around the LHS result (p1 log-spaced, p2 linear) and restarts Nelder-Mead
# from the grid minimum. Returns ([p1, p2, c], chi2) for whichever of the
# LHS / polished solutions is lower.
function polish(pred, p_lhs, chi2_lhs; n = 41, p1_half = 0.4, p2_frac = 0.3)
    p10, p20 = p_lhs[1], p_lhs[2]
    p1_vec = 10 .^ range(log10(p10) - p1_half, log10(p10) + p1_half, length = n)
    p2_vec = range(p20 * (1 - p2_frac), p20 * (1 + p2_frac), length = n)
    G = [joint_chi2(pred, p1, p2)[1] for p1 in p1_vec, p2 in p2_vec]
    imin = argmin(G)
    any(i -> i in (1, n), Tuple(imin)) && @warn "grid minimum on the grid edge -- widen the polish grid" p1=p1_vec[imin[1]] p2=p2_vec[imin[2]]
    res = optimize(q -> joint_chi2(pred, q[1] * p10, q[2])[1], [p1_vec[imin[1]] / p10, p2_vec[imin[2]]], NelderMead())
    p1_p, p2_p = Optim.minimizer(res)[1] * p10, Optim.minimizer(res)[2]
    chi2_p = joint_chi2(pred, p1_p, p2_p)[1]
    @printf("    grid polish: chi2 = %.3f (LHS: %.3f)\n", chi2_p, chi2_lhs)
    p1, p2, chi2 = chi2_p < chi2_lhs ? (p1_p, p2_p, chi2_p) : (p_lhs[1], p_lhs[2], chi2_lhs)
    [p1, p2, joint_chi2(pred, p1, p2)[2]], chi2
end

# ── EV (viscous CSS) ───────────────────────────────────────────────────────────
pred_css(x, tau_i, tau_f, a, b) = diff_css(x, tau_i, tau_f, b, v_plate_arg, a)

println("Fitting EV RSS (T = b/V, V = $(V_lit*yr*100) cm/yr) jointly to Foulger (1992) + Hofton & Foulger (1996) …")
p_ev, chi2_ev_val = fit_np(pred_css, [1e-7, 0.1], [1e-2, 3.0]; n = 8000, n_top = 5)
p_ev, chi2_ev_val = polish(pred_css, p_ev, chi2_ev_val)
a_ev, b_ev, c_ev = p_ev
T_ev = b_ev / v_plate_arg

# ── EVE (Maxwell CSS / MaxCSS), tau_R derived from a ──────────────────────────
pred_maxcss(x, tau_i, tau_f, a, b; k_limit_scale=1.0) = begin
    tau_R = tau_R_from_a(a)
    get_u_css_maxwell_closed.(x, tau_f, b, v_plate_arg, a, tau_R; k_limit_scale=k_limit_scale) .-
    get_u_css_maxwell_closed.(x, tau_i, b, v_plate_arg, a, tau_R; k_limit_scale=k_limit_scale)
end

println("Fitting EVE RSS (T = b/V) jointly to Foulger (1992) + Hofton & Foulger (1996) …")
p_eve, chi2_eve_val = fit_np(pred_maxcss, [1e-7, 0.1], [1e-2, 3.0]; n = 10000, n_top = 3)
p_eve, chi2_eve_val = polish(pred_maxcss, p_eve, chi2_eve_val; n = 21)
a_eve, b_eve, c_eve = p_eve
T_eve   = b_eve / v_plate_arg
tauR_eve = tau_R_from_a(a_eve)

# ── Foulger diffusion model: published (per-paper) and our joint re-fit ───────
# Published: Foulger et al. (1992) kappa=1.1, U0=1.0 (no offset -- their own
# near-rift reference points); Hofton & Foulger (1996) kappa=3.11, U0=1.27.
# Both papers report the dike half-width U0; our b is the full opening, b = 2 U0
# (their own zero is arbitrary, so only the offset is solved, kappa/b fixed).
kappa_pub_foulger, b_pub_foulger, c_pub_foulger = 1.1, 2 * 1.0, 0.0
kappa_pub_hofton,  b_pub_hofton                 = 3.11, 2 * 1.27
c_pub_hofton = best_c(diff_foulger(obs_x_m_hofton, tau_i_hofton, tau_f_hofton, b_pub_hofton, kappa_pub_hofton))

# Joint re-fit: one shared (kappa, b) across both datasets, Hofton offset
# solved analytically -- same search as EV/EVE so the comparison is fair.
pred_diff(x, tau_i, tau_f, kappa, b) = diff_foulger(x, tau_i, tau_f, b, kappa)
println("Fitting joint Foulger diffusion model (kappa, b; c analytic) …")
res_diff = optimize(q -> joint_chi2(pred_diff, q[1], q[2])[1], [0.05, 0.2], [10.0, 10.0], [1.86, 2.08], Fminbox(NelderMead()))
p_diff, chi2_diff_val = polish(pred_diff, Optim.minimizer(res_diff), Optim.minimum(res_diff))
kappa_new, b_new, c_new_hofton = p_diff
c_new_foulger = 0.0

# ── Report ───────────────────────────────────────────────────────────────────
nu = N_total - 3   # a, b, c free; T (and tau_R for EVE) derived, not fit

@printf("\n=== V-constrained joint RSS fit (N = %d, V = %.1f cm/yr) ===\n", N_total, V_lit * yr * 100)
@printf("  EV RSS:   a = %.3e m/s   b = %.3f m   T = %.1f yr (derived)   chi2/nu = %.3f\n",
        a_ev, b_ev, T_ev/yr, chi2_ev_val/nu)
@printf("            c (Hofton offset, analytic) = %.3f m\n", c_ev)
@printf("  EVE RSS:  a = %.3e m/s   b = %.3f m   T = %.1f yr (derived)   tau_R = %.2f yr (derived)   chi2/nu = %.3f\n",
        a_eve, b_eve, T_eve/yr, tauR_eve/yr, chi2_eve_val/nu)
@printf("            c (Hofton offset, analytic) = %.3f m\n", c_eve)
@printf("  Foulger diffusion, joint re-fit: kappa = %.3f m^2/s   b = %.3f m   c (Hofton) = %.3f m   chi2/nu = %.3f\n",
        kappa_new, b_new, c_new_hofton, chi2_diff_val/nu)

# Best-fit parameters, read back by src/figures/krafla_v_constrained_farfield.jl
mkpath("reports")
params_csv = "reports/krafla_css_v_constrained_params$(fig_suffix).csv"
CSV.write(params_csv, DataFrame(
    V_cm_yr = V_lit * yr * 100, model = ["EV", "EVE"],
    a = [a_ev, a_eve], b = [b_ev, b_eve], c_hofton = [c_ev, c_eve],
    chi2 = [chi2_ev_val, chi2_eve_val]))
println("$params_csv  saved")

# ── Figure: one panel per dataset ─────────────────────────────────────────────
x_km = range(-180, stop=180, length=2000)
x_m  = x_km .* 1e3

function plot_panel!(ax; obs_x_km, obs_u, obs_sig, tau_i, tau_f, data_label,
        c_ev_panel, c_eve_panel, kappa_pub, b_pub, c_pub, c_new, show_diffusion=true)
    lines!(ax, x_km, (pred_css(x_m, tau_i, tau_f, a_ev, b_ev) .+ c_ev_panel) .* 100,
        color=:purple, linewidth=2.5,
        label="EV RSS | a=$(sci(a_ev)) m/s  b=$(round(b_ev,digits=2)) m  T=$(round(T_ev/yr,digits=0)) yr")
    lines!(ax, x_km, (pred_maxcss(x_m, tau_i, tau_f, a_eve, b_eve; k_limit_scale=30.0) .+ c_eve_panel) .* 100,
        color=:navy, linewidth=2.5,
        label="EVE RSS | a=$(sci(a_eve)) m/s  b=$(round(b_eve,digits=2)) m  T=$(round(T_eve/yr,digits=0)) yr  τ_R=$(round(tauR_eve/yr,digits=2)) yr")
    if show_diffusion
        lines!(ax, x_km, (diff_foulger(x_m, tau_i, tau_f, b_pub, kappa_pub) .+ c_pub) .* 100,
            color=:red, linewidth=2, linestyle=:dash,
            label="Foulger diffusion (published) | κ=$(sci(kappa_pub)) m²/s  b=$(round(b_pub,digits=2)) m")
        lines!(ax, x_km, (diff_foulger(x_m, tau_i, tau_f, b_new, kappa_new) .+ c_new) .* 100,
            color=:darkorange, linewidth=2, linestyle=:dashdot,
            label="Foulger diffusion (joint re-fit) | κ=$(sci(kappa_new)) m²/s  b=$(round(b_new,digits=2)) m")
    end

    scatter!(ax, obs_x_km, obs_u .* 100, color=:black, marker=:circle, markersize=10, label=data_label)
    errorbars!(ax, obs_x_km, obs_u .* 100, obs_sig .* 100, color=:black)
    return ax
end

function make_figure(path; show_diffusion=true)
    fig = Figure(size=(2600, 1000), fontsize=46.43)
    ax1 = Axis(fig[1, 1], xlabel = "Distance from Krafla volcanic system (km)",
        ylabel = "Displacement (cm)", limits = (-180, 180, -15, 20),
        title = "Foulger et al. (1992), τ=$(round(Int, tau_i_foulger/yr))–$(round(Int, tau_f_foulger/yr)) yr")
    ax2 = Axis(fig[1, 2], xlabel = "Distance from Krafla volcanic system (km)",
        ylabel = "Displacement (cm)", limits = (-180, 180, -15, 20),
        title = "Hofton & Foulger (1996), τ=$(round(Int, tau_i_hofton/yr))–$(round(Int, tau_f_hofton/yr)) yr")

    plot_panel!(ax1; obs_x_km=obs_x_m_foulger ./ 1e3, obs_u=obs_u_foulger, obs_sig=obs_sig_foulger,
        tau_i=tau_i_foulger, tau_f=tau_f_foulger, data_label="Foulger et al. (1992)",
        c_ev_panel=0.0, c_eve_panel=0.0, kappa_pub=kappa_pub_foulger, b_pub=b_pub_foulger, c_pub=c_pub_foulger, c_new=c_new_foulger,
        show_diffusion=show_diffusion)

    plot_panel!(ax2; obs_x_km=obs_x_m_hofton ./ 1e3, obs_u=obs_u_hofton, obs_sig=obs_sig_hofton,
        tau_i=tau_i_hofton, tau_f=tau_f_hofton, data_label="Hofton & Foulger (1996), Fig. 5",
        c_ev_panel=c_ev, c_eve_panel=c_eve, kappa_pub=kappa_pub_hofton, b_pub=b_pub_hofton, c_pub=c_pub_hofton, c_new=c_new_hofton,
        show_diffusion=show_diffusion)

    axislegend(ax1, position=:rb, framevisible=true, labelsize=25)
    axislegend(ax2, position=:rb, framevisible=true, labelsize=25)

    save(path, fig, px_per_unit=2)
    println("\n$path saved")
end

make_figure("fig/fig_krafla_css_v_constrained$(fig_suffix).png")
# Same fits, EV/EVE only (no Foulger diffusion curves)
make_figure("fig/fig_krafla_css_v_constrained_models_only$(fig_suffix).png"; show_diffusion=false)
