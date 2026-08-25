# Viscous forward models: absolute displacements, velocities, and
# differenced displacements between two observation epochs.
#
# Dependencies (must be loaded before include):
#   SpecialFunctions (erfc, loggamma, digamma)

# ── Absolute displacement profiles ────────────────────────────────────────────

function get_u_css(x, tau, b, v, a)
    T  = b / v
    z  = (tau / T) + im * (abs(x) / (a * T))
    z0 = im * (abs(x) / (a * T))
    return sign(x) * (2 * b / pi) * imag(loggamma(z) - loggamma(z0))
end

function get_u_step(x, tau, b, a)
    return (b / pi) * atan.(a * tau ./ x)
end

function get_u_foulger(x, tau, kappa, b)
    z = x / (2 * sqrt(kappa * tau))
    return b * erfc(z)
end

# ── Velocity profiles ──────────────────────────────────────────────────────────

function get_v_css(x, tau, b, v_plate, a)
    T = b / v_plate
    z = (tau / T) + im * (abs(x) / (a * T))
    return sign(x) * (2 * b / (pi * T)) * imag(digamma(z))
end

function get_v_step(x, tau, b, a)
    return (b * a / pi) * (x / (x^2 + (a * tau)^2))
end

function get_v_foulger(x, tau, kappa, b)
    (kappa <= 0 || tau <= 0) && return 0.0
    exp1 = -x^2 / (4 * kappa * tau)
    return (b / (tau * sqrt(pi))) * (x / (2 * sqrt(kappa * tau))) * exp(exp1)
end

# ── Differenced displacements (prediction between two epochs) ─────────────────

diff_step(x, tau_i, tau_f, b, a) =
    get_u_step.(x, tau_f, b, a) .- get_u_step.(x, tau_i, b, a)

diff_foulger(x, tau_i, tau_f, b, kappa) =
    get_u_foulger.(x, tau_f, kappa, b) .- get_u_foulger.(x, tau_i, kappa, b)

diff_css(x, tau_i, tau_f, b, v_plate, a) =
    get_u_css.(x, tau_f, b, v_plate, a) .-
    get_u_css.(x, tau_i, b, v_plate, a)
