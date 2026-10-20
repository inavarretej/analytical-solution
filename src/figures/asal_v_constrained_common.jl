using CSV, DataFrames, Printf

include("../models/viscous.jl")
include("../models/maxwell.jl")
include("../parameters.jl")   # H, E_L, nu_L -- shared plate/half-space modulus assumption
sci(x) = @sprintf("%.2e", x)

# Shared by src/figures/asal_css_v_constrained.jl (fits) and
# src/figures/asal_chi2_landscape_v_constrained.jl (landscapes).
#
# CSS/MaxCSS joint fit (EP00-DF00 time series + rift-perpendicular velocity
# profile) to Asal with the recurrence interval T NOT free: T = b/V, with V
# a literature spreading rate -- same approach as the Krafla default
# (src/figures/krafla_css_v_constrained.jl).
#
# The model functions take a "v_plate" argument and set T = b/v_plate
# internally, with far-field velocity +-v_plate/2 on each flank (v_RSS is odd
# in x). V is a full (two-flank) opening rate, so it is passed directly.
#
# Velocity: V = 13 +- 1 mm/yr, the long-term Arabia-Somalia rate of Vigny
# et al. (2007) -- the analogue of the NUVEL-1A plate rate used for Krafla.
# Not their 16 +- 1 mm/yr opening across the Asal rift (GPS 1991-2003): the
# authors attribute that excess to post-1978 transient relaxation, which the
# RSS model already produces itself, so using it for V would count the
# transient twice. (V = 16 was tried on 2026-09-23 as a sensitivity test:
# better chi^2, since it is measured on the same profile being fit, T ~ 40 yr
# instead of ~50 yr, same EV-over-EVE preference.)
#
# Background offset of the profile: Vigny's profile is in the Somalia frame,
# where the SW flank is ~fixed and the NE flank moves at the full rate. The
# model's far-field asymptotes are -V/2 and +V/2, so the frame shift is
# V_pl = V/2 by construction -- no longer a fitted nuisance parameter (it was
# in src/fits/asal.jl, where it came out at ~8 mm/yr, i.e. ~16/2).
yr      = 3600.0 * 24 * 365.25
t0_year = 1978.85               # rifting onset (Vigny et al. 2007)

V_cases = [(label = "V = 13 mm/yr (Arabia–Somalia plate rate)", V = 0.013 / yr)]

tau_R_from_a(a) = H / (a * (1 - nu_L))

# ── Data ─────────────────────────────────────────────────────────────────────
df       = CSV.read("data/asal_timeseries_merged.csv", DataFrame)
obs_year = df.year
obs_dL   = df.dL_mm
obs_sig  = df.sig_mm
obs_t    = (obs_year .- t0_year) .* yr
N_ts     = length(obs_year)
x0_m     = 2000.0
t_ref_s  = (obs_year[1] - t0_year) * yr

df_prof   = CSV.read("data/asal_profile.csv", DataFrame)
prof_x_km = df_prof.X
prof_v    = df_prof.Y
prof_sig  = df_prof.SD
prof_x_m  = prof_x_km .* 1e3
N_prof    = length(prof_x_km)

t_epoch_profile = 19.0 * yr   # Vigny et al. (2007) profile mean epoch
x_center_m      = 1000.0      # model centre, +1 km NE of the topographic rift axis
x_rel(x_m)      = x_m .- x_center_m

N_total = N_ts + N_prof

# No second intrusion recorded since 1978 through 2016 (Smittarello et al.
# 2016) -- same floor as src/fits/asal.jl. With T = b/V this becomes a lower
# bound on b: b >= V * T_min_no_repeat.
T_min_no_repeat = (2016.0 - t0_year) * yr

weighted_chi2(pred, obs, sig) = sum(((obs .- pred) ./ sig) .^ 2)

# ── Predictions (mod T: the closed forms are valid only for tau in [0, T)) ────
pred_ts_ev(t, a, b, V)  = 2 .* get_u_css.(x0_m, mod.(t, b / V), b, V, a) .* 1000
pred_prof_ev(x, a, b, V) = get_v_css.(x_rel(x), mod(t_epoch_profile, b / V), b, V, a) .* 1000 .* yr .+ V * yr * 1000 / 2

# EVE uses the converged split forms (2026-10-16). The closed forms at the
# default k_limit_scale = 1 under-truncate the k-integral and biased the old
# best-fit χ² low by 0.8 (434.3 vs 435.1), almost all from the time series.
# k_limit_scale is accepted and ignored so the plotting calls still work.
pred_ts_eve(t, a, b, V; k_limit_scale=nothing) =
    2 .* get_u_css_maxwell_split.(x0_m, mod.(t, b / V), b, V, a, tau_R_from_a(a)) .* 1000
pred_prof_eve(x, a, b, V; k_limit_scale=nothing) =
    get_v_css_maxwell_split.(x_rel(x), mod(t_epoch_profile, b / V), b, V, a, tau_R_from_a(a)) .* 1000 .* yr .+
    V * yr * 1000 / 2

tref(f, t) = f(t) .- f([t_ref_s])[1]

# Unpenalized joint chi^2 (used directly by the landscapes, which mask the
# T floor themselves rather than plotting a 1e12 wall).
joint_chi2_raw(pred_ts, pred_prof, a, b, V) =
    weighted_chi2(tref(t -> pred_ts(t, a, b, V), obs_t), obs_dL, obs_sig) +
    weighted_chi2(pred_prof(prof_x_m, a, b, V), prof_v, prof_sig)

function joint_chi2(pred_ts, pred_prof, a, b, V)
    (a <= 0 || b <= 0 || a > 1e-2 || b > 3.0 || b / V < T_min_no_repeat) && return 1e12
    joint_chi2_raw(pred_ts, pred_prof, a, b, V)
end

params_csv = "reports/asal_css_v_constrained_params.csv"
