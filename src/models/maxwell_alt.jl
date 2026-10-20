# Alternative EVE RSS (vault note "EVE initial condition - missing
# instantaneous elastic response", Step 8): the state at t = 0+ is obtained by
# integrating the visco-elastic kernel across each intrusion, which gives
#   u_hat(k, 0+) = u_hat(k, T-) + (ib/k) / (1 + l_M k),   l_M = a tau_R.
# Compared with the paper's EVE RSS (get_*_css_maxwell_split in maxwell.jl),
# the relaxing amplitude carries one extra factor 1/(1 + l_M k):
#   v_alt = (b a / pi) int sin(kx) / (1 + l_M k)^2 * G(k) dk,
#   u_alt = (b / pi)   int sin(kx) / (k (1 + l_M k)) * F(k) dk,
# G = exp(-alpha tau)/(1 - exp(-alpha T)), F = (1 - exp(-alpha tau))/(1 - exp(-alpha T)).
# u_alt is the post-diking displacement since the most recent intrusion (it
# excludes the co-diking elastic step), so it is what differenced geodetic
# epochs within one cycle measure. This file does not modify maxwell.jl; it
# only reuses its helpers.
#
# As k -> inf, G -> G_inf and F -> F_inf. The parts with G_inf, F_inf are
# integrated in closed form with the auxiliary functions of Ci and Si,
#   int_0^inf sin(z s)/(1+z)   dz = f(s) = Ci(s) sin s - (Si(s) - pi/2) cos s,
#   int_0^inf sin(z s)/(1+z)^2 dz = s g(s),  g(s) = -Ci(s) cos s - (Si(s) - pi/2) sin s,
# so the remainders decay as 1/k^3 and a fixed cutoff converges for any x.
using QuadGK, SpecialFunctions

_f_alt(s) = cosint(s) * sin(s) - (sinint(s) - pi / 2) * cos(s)
_g_alt(s) = -cosint(s) * cos(s) - (sinint(s) - pi / 2) * sin(s)

# int_0^inf sin(kx) / (k (1 + l k)) dk = sign(x) [pi/2 - f(|x|/l)]
_I1_alt(x, l) = x == 0 ? 0.0 : sign(x) * (pi / 2 - _f_alt(abs(x) / l))
# int_0^inf sin(kx) / (1 + l k)^2 dk = sign(x) (|x|/l) g(|x|/l) / l
_I2_alt(x, l) = x == 0 ? 0.0 : sign(x) * (abs(x) / l) * _g_alt(abs(x) / l) / l

function get_v_css_maxwell_alt(x, tau, b, v_plate, a, tau_R)
    T = b / v_plate
    l = a * tau_R
    G_inf = exp(-tau / tau_R) / -expm1(-T / tau_R)
    remainder(k) = begin
        α = _alpha_css_maxwell(k, a, tau_R)
        sin(k * x) / (1 + l * k)^2 * (exp(-α * tau) / -expm1(-α * T) - G_inf)
    end
    val = _segmented_quadgk(remainder, 1e4 / l, x)
    return (b * a / pi) * (G_inf * _I2_alt(x, l) + val)
end

function get_u_css_maxwell_alt(x, tau, b, v_plate, a, tau_R)
    T = b / v_plate
    l = a * tau_R
    F_inf = expm1(-tau / tau_R) / expm1(-T / tau_R)
    remainder(k) = begin
        α = _alpha_css_maxwell(k, a, tau_R)
        sin(k * x) / (k * (1 + l * k)) * (expm1(-α * tau) / expm1(-α * T) - F_inf)
    end
    val = _segmented_quadgk(remainder, 1e4 / l, x)
    return (b / pi) * (F_inf * _I1_alt(x, l) + val)
end
