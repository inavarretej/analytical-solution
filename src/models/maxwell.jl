# Maxwell viscoelastic forward models: absolute displacements, velocities,
# and differenced displacements between two observation epochs.
#
# Dependencies (must be loaded before include):
#   QuadGK (quadgk)

function _k_limit_maxwell(x, tau, a, tau_R; epsilon=1e-4)
    k_visc = 20.0 / (a * tau + 1e-12)
    k_near = 100.0 / max(abs(x), 1.0)
    k_tail = exp(-tau / tau_R) / (epsilon * max(abs(x), 1.0))
    return max(k_visc, k_near, k_tail)
end

function _segmented_quadgk(f, k_limit, x; rtol=1e-10, atol=1e-12)
    n_max = min(floor(Int, k_limit * abs(x) / pi), 300)
    bpts  = [n * pi / abs(x) for n in 1:n_max]
    segs  = [0.0; bpts; k_limit]
    val   = 0.0
    for i in 1:(length(segs) - 1)
        v, _ = quadgk(f, segs[i], segs[i+1]; rtol=rtol, atol=atol)
        val += v
    end
    return val
end

function get_u_viscoelastic(x, tau, b, a, tau_R)
    (a <= 0 || tau < 0 || tau_R <= 0) && return 0.0
    abs(x) < 1.0 && return sign(x) * b / 2.0

    k_limit   = _k_limit_maxwell(x, tau, a, tau_R)
    integrand(k) = sin(k * x) / k * exp(-(a * k * tau) / (1 + tau_R * a * k))
    val = _segmented_quadgk(integrand, k_limit, x)
    return sign(x) * (b / 2) - (b / pi) * val
end

function get_v_relax(x, tau, b, a, tau_R)
    (a <= 0 || tau < 0 || tau_R <= 0) && return 0.0

    k_limit   = _k_limit_maxwell(x, tau, a, tau_R)
    integrand(k) = sin(k * x) / (1 + tau_R * a * k) *
                   exp(-(a * k * tau) / (1 + tau_R * a * k))
    val = _segmented_quadgk(integrand, k_limit, x)
    return (b * a / pi) * val
end

# ── CSS Maxwell: sum of all past Maxwell responses (u replaces t with τ−jT) ──

function get_u_css_maxwell(x, tau, b, v_plate, a, tau_R)
    T = b / v_plate
    N = floor(Int, tau / T)
    total = 0.0
    for j in 0:N
        total += get_u_viscoelastic(x, tau - j * T, b, a, tau_R)
    end
    return total
end

function get_v_css_maxwell(x, tau, b, v_plate, a, tau_R)
    T = b / v_plate
    N = floor(Int, tau / T)
    total = 0.0
    for j in 0:N
        total += get_v_relax(x, tau - j * T, b, a, tau_R)
    end
    return total
end

# ── CSS Maxwell, closed periodic form (paper eq:maxwell_css_v / eq:maxwell_css_u) ──
# Distinct from get_u/v_css_maxwell above: those sum a FINITE real-space history
# starting at t=0 (an approximation for "few events on record"); these instead
# impose the Fourier periodicity condition u_hat(k,tau+T) = u_hat(k,tau) + ib/k
# directly (infinite periodic history), matching how get_u_css/get_v_css (viscous)
# are themselves derived. Verified to converge to get_u_css/get_v_css as tau_R -> 0.
# tau is the phase since the most recent intrusion, tau in [0,T).

_alpha_css_maxwell(k, a, tau_R) = a * k / (1 + tau_R * a * k)

function get_v_css_maxwell_closed(x, tau, b, v_plate, a, tau_R; k_limit_scale=1.0)
    T = b / v_plate
    k_limit = k_limit_scale * max(_k_limit_maxwell(x, tau, a, tau_R), 50.0 / (a * T))
    integrand(k) = begin
        α = _alpha_css_maxwell(k, a, tau_R)
        sin(k * x) / (1 + tau_R * a * k) * exp(-α * tau) / (1 - exp(-α * T))
    end
    val = _segmented_quadgk(integrand, k_limit, x)
    return (2 * b * a / pi) * val
end

# k_limit_scale: the default _k_limit_maxwell cutoff, tuned for get_u_viscoelastic,
# under-truncates this closed-periodic integrand when T is large (e.g. a
# multi-century CSS recurrence interval fit to a short observation window):
# the (1-exp(-alpha*tau))/(1-exp(-alpha*T)) factor has a longer effective
# k-tail than the plain exponential decay the default cutoff assumes. At the
# default scale=1 this produces a spurious non-monotonic ripple of a few
# tenths of a mm (visible when plotted on a fine time grid) -- confirmed by
# checking monotonicity of u(tau) at fixed x, which requires scale >~30 to
# eliminate entirely. Since raising the scale costs ~40-60x per evaluation,
# leave the default at 1 (fast, adequate for chi^2 optimization, where this
# noise does not affect the fitted parameters) and pass a higher scale only
# when generating the smooth curve for a figure.
function get_u_css_maxwell_closed(x, tau, b, v_plate, a, tau_R; k_limit_scale=1.0)
    T = b / v_plate
    k_limit = k_limit_scale * max(_k_limit_maxwell(x, tau, a, tau_R), 50.0 / (a * T))
    integrand(k) = begin
        α = _alpha_css_maxwell(k, a, tau_R)
        sin(k * x) / k * (1 - exp(-α * tau)) / (1 - exp(-α * T))
    end
    val = _segmented_quadgk(integrand, k_limit, x)
    return (2 * b / pi) * val
end

# ── Differenced displacements (prediction between two epochs) ─────────────────

diff_viscoelastic(x, tau_i, tau_f, b, a, tau_R) =
    get_u_viscoelastic.(x, tau_f, b, a, tau_R) .-
    get_u_viscoelastic.(x, tau_i, b, a, tau_R)

diff_css_maxwell(x, tau_i, tau_f, b, v_plate, a, tau_R) =
    get_u_css_maxwell.(x, tau_f, b, v_plate, a, tau_R) .-
    get_u_css_maxwell.(x, tau_i, b, v_plate, a, tau_R)
