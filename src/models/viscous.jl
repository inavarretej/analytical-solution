# Viscous forward models: absolute displacements, velocities, and
# differenced displacements between two observation epochs.
#
using SpecialFunctions

# ── Absolute displacement profiles ────────────────────────────────────────────

# RSS (Recurrent Steady State): infinite periodic sum of STEP responses,
# elapsed times tau, tau+T, tau+2T, ... since each past intrusion. Reducing
# that sum via the Digamma/Gamma functions gives prefactor b/pi (NOT 2b/pi --
# an earlier derivation mis-indexed the sum and doubled the amplitude; see
# vault note "Factor-of-2 error in the EV RSS velocity formula" / "From STEP
# to RSS velocity (corrected derivation)"). With this prefactor,
# int_0^T v_RSS dtau = b/2 per flank, matching the STEP far-field plateau
# u^STEP(x,inf)=b/2 exactly -- the RSS never contradicts the single-event
# limit, it reproduces it once per cycle.
#
# T = b/v is kept as the defining relation (unchanged): callers that already
# construct v = b/T_target to select a target recurrence interval continue to
# get exactly T_target here, unaffected by this fix. What changes is that the
# model's own far-field asymptotic velocity is now v/2, not v -- so a caller
# that wants the true physical plate rate V to be reproduced at x->infinity
# must pass v = 2V, not v = V.
function get_u_css(x, tau, b, v, a)
    T  = b / v
    z  = (tau / T) + im * (abs(x) / (a * T))
    z0 = im * (abs(x) / (a * T))
    return sign(x) * (b / pi) * imag(loggamma(z) - loggamma(z0))
end

function get_u_step(x, tau, b, a)
    return (b / pi) * atan.(a * tau ./ x)
end

# Foulger et al. (1992) diffusion model, written in terms of the full dike
# opening b (their dike of width 2 U0, so U0 = b/2): u = (b/2) erfc(...).
function get_u_foulger(x, tau, kappa, b)
    z = x / (2 * sqrt(kappa * tau))
    return (b / 2) * erfc(z)
end

# ── Velocity profiles ──────────────────────────────────────────────────────────

function get_v_css(x, tau, b, v_plate, a)
    T = b / v_plate
    z = (tau / T) + im * (abs(x) / (a * T))
    return sign(x) * (b / (pi * T)) * imag(digamma(z))
end

function get_v_step(x, tau, b, a)
    return (b * a / pi) * (x / (x^2 + (a * tau)^2))
end

function get_v_foulger(x, tau, kappa, b)
    (kappa <= 0 || tau <= 0) && return 0.0
    exp1 = -x^2 / (4 * kappa * tau)
    return (b / (2 * tau * sqrt(pi))) * (x / (2 * sqrt(kappa * tau))) * exp(exp1)
end

# ── Differenced displacements (prediction between two epochs) ─────────────────

diff_step(x, tau_i, tau_f, b, a) =
    get_u_step.(x, tau_f, b, a) .- get_u_step.(x, tau_i, b, a)

diff_foulger(x, tau_i, tau_f, b, kappa) =
    get_u_foulger.(x, tau_f, kappa, b) .- get_u_foulger.(x, tau_i, kappa, b)

diff_css(x, tau_i, tau_f, b, v_plate, a) =
    get_u_css.(x, tau_f, b, v_plate, a) .-
    get_u_css.(x, tau_i, b, v_plate, a)
