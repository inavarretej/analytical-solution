# Maxwell viscoelastic forward models: absolute displacements, velocities,
# and differenced displacements between two observation epochs.
#
using QuadGK
using SpecialFunctions  # sinint, cosint (tau -> 0+ limits)

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

# ── Limits tau -> 0+ (right after an intrusion) ──────────────────────────────
# Used for the first row of the EVE contour figures (Figs. 7 and 8): at very
# small tau the exponential factor exp(-alpha tau) barely decays with k, the
# adaptive k cutoff becomes huge, and get_v_relax / get_v_css_maxwell_closed
# cost seconds per call (e.g. ~7 s per interior-kernel evaluation at
# tau = 0.001 tau_R). The 0+ limits below are cheap and exact. Full
# documentation (derivation, checks, timings): Obsidian note
# 30-models/analytical/EVE contour figures - tau = 0+ limit.md

# Auxiliary sine/cosine-integral function f(z) = int_0^inf sin(t)/(t + z) dt
# = Ci(z) sin z - (Si(z) - pi/2) cos z   (Abramowitz & Stegun 5.2.6/5.2.12).
_f_aux(z) = cosint(z) * sin(z) - (sinint(z) - pi / 2) * cos(z)

# EVE STEP velocity at tau = 0+: with exp(-alpha tau) = 1,
#   v = (b a/pi) int_0^inf sin(kx)/(1 + tau_R a k) dk
#     = sign(x) (b/(pi tau_R)) f(|x|/ell_M),   ell_M = a tau_R.
# Checked against get_v_relax at tau = 0.001 tau_R: relative difference
# <= 1e-3 (x = 100 m .. 300 km), i.e. O(tau/tau_R) as expected.
function get_v_relax_0plus(x, b, a, tau_R)
    x == 0 && return 0.0
    return sign(x) * (b / (pi * tau_R)) * _f_aux(abs(x) / (a * tau_R))
end

# EVE RSS (closed periodic form) velocity at tau = 0+. With
# G(k) = 1/(1 - exp(-alpha T)) and alpha -> 1/tau_R as k -> inf,
#   G_inf = 1/(1 - exp(-T/tau_R)),
#   v = G_inf v_STEP(x, 0+) + (b a/pi) int_0^inf sin(kx) [G(k) - G_inf]/(1 + tau_R a k) dk.
# G - G_inf = O(1/k), so the remainder integrand decays as sin(kx)/k^2 and
# converges absolutely. Checked against get_v_css_maxwell_closed at
# tau = 1e-3 T and 2e-4 T: relative differences shrink ~linearly with tau
# (e.g. 4.3e-3 -> 9.5e-4 at x = 1 km), i.e. convergence to this limit.
# Slow (~1-4 s per call, oscillatory infinite k-integral): tabulate it (see
# the EVE RSS contour script) rather than calling it inside the kernels.
function get_v_css_maxwell_0plus(x, b, v_plate, a, tau_R)
    x == 0 && return 0.0
    T     = b / v_plate
    G_inf = 1 / (1 - exp(-T / tau_R))
    remainder(k) = begin
        α = _alpha_css_maxwell(k, a, tau_R)
        G = -1 / expm1(-α * T)
        sin(k * abs(x)) * (G - G_inf) / (1 + tau_R * a * k)
    end
    rem, _ = quadgk(remainder, 0, Inf; rtol = 1e-8)
    return G_inf * get_v_relax_0plus(x, b, a, tau_R) + sign(x) * (b * a / pi) * rem
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
#
# Prefactor is b (not 2b), mirroring the same correction applied to
# get_u_css/get_v_css in viscous.jl -- this file was calibrated to converge to
# those functions as tau_R -> 0, so it inherited their factor-of-2 indexing
# bug and needs the same fix. T = b/v_plate is kept unchanged (see viscous.jl
# for why); see vault note "Factor-of-2 error in the EV RSS velocity formula".

_alpha_css_maxwell(k, a, tau_R) = a * k / (1 + tau_R * a * k)

function get_v_css_maxwell_closed(x, tau, b, v_plate, a, tau_R; k_limit_scale=1.0)
    T = b / v_plate
    k_limit = k_limit_scale * max(_k_limit_maxwell(x, tau, a, tau_R), 50.0 / (a * T))
    integrand(k) = begin
        α = _alpha_css_maxwell(k, a, tau_R)
        sin(k * x) / (1 + tau_R * a * k) * exp(-α * tau) / (1 - exp(-α * T))
    end
    val = _segmented_quadgk(integrand, k_limit, x)
    return (b * a / pi) * val
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
    return (b / pi) * val
end

# Same displacement as get_u_css_maxwell_closed, valid down to x -> 0+. The
# factor F(k) = (1 - exp(-alpha tau))/(1 - exp(-alpha T)) tends to
# F_inf = (1 - exp(-tau/tau_R))/(1 - exp(-T/tau_R)) as k -> inf, so the plain
# k-integral needs k >> 1/|x| and is under-truncated below ~1 m. Splitting
#   u = sign(x) (b/2) F_inf + (b/pi) int_0^inf sin(kx)/k [F(k) - F_inf] dk
# leaves a remainder that decays as 1/k^2 on the scale k ~ 1/(a tau_R), so a
# fixed cutoff 1e4/(a tau_R) converges for any x. Checked against
# get_u_css_maxwell_closed (k_limit_scale = 100) at T = 5 tau_R,
# x = 1e-4 .. 3 aT: differences <= 3e-5 b; at x = 1e-12 aT it equals the x = 0+
# closed form. Also ~10-100x faster.
function get_u_css_maxwell_split(x, tau, b, v_plate, a, tau_R)
    T = b / v_plate
    F_inf = expm1(-tau / tau_R) / expm1(-T / tau_R)
    remainder(k) = begin
        α = _alpha_css_maxwell(k, a, tau_R)
        sin(k * x) / k * (expm1(-α * tau) / expm1(-α * T) - F_inf)
    end
    val = _segmented_quadgk(remainder, 1e4 / (a * tau_R), x)
    return sign(x) * (b / 2) * F_inf + (b / pi) * val
end

# Velocity counterpart of get_u_css_maxwell_split. As k -> inf the integrand of
# get_v_css_maxwell_closed tends to sin(kx) G_inf / (tau_R a k), with
# G_inf = exp(-tau/tau_R)/(1 - exp(-T/tau_R)); its Dirichlet integral gives the
# axis jump sign(x) b G_inf/(2 tau_R) (paper eq:maxwell_rss_jump). Subtracting
# it leaves a remainder that decays as 1/k^2 on the scale k ~ 1/(a tau_R), so
# the same fixed cutoff 1e4/(a tau_R) converges for any x.
function get_v_css_maxwell_split(x, tau, b, v_plate, a, tau_R)
    T = b / v_plate
    G_inf = exp(-tau / tau_R) / -expm1(-T / tau_R)
    remainder(k) = begin
        α = _alpha_css_maxwell(k, a, tau_R)
        sin(k * x) * (exp(-α * tau) / -expm1(-α * T) / (1 + tau_R * a * k) -
                      G_inf / (tau_R * a * k))
    end
    val = _segmented_quadgk(remainder, 1e4 / (a * tau_R), x)
    return sign(x) * b * G_inf / (2 * tau_R) + (b * a / pi) * val
end

# ── Differenced displacements (prediction between two epochs) ─────────────────

diff_viscoelastic(x, tau_i, tau_f, b, a, tau_R) =
    get_u_viscoelastic.(x, tau_f, b, a, tau_R) .-
    get_u_viscoelastic.(x, tau_i, b, a, tau_R)

diff_css_maxwell(x, tau_i, tau_f, b, v_plate, a, tau_R) =
    get_u_css_maxwell.(x, tau_f, b, v_plate, a, tau_R) .-
    get_u_css_maxwell.(x, tau_i, b, v_plate, a, tau_R)
