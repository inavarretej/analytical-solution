using CairoMakie, SpecialFunctions, QuadGK
using LaTeXStrings

include("../models/viscous.jl")
include("../models/maxwell.jl")

# ── Dimensionless units: a = T = b = v_plate = 1 ──────────────────────────────
# x̂ = x/(aT),  τ̂ = τ/T ∈ [0,1),  u normalised by b,  v normalised by b/T
a_d  = 1.0;  T_d = 1.0;  b_d = 1.0;  vp_d = b_d / T_d

x_hat = collect(range(0.05, 5.0, length = 400))

# ── Panel (a): CSS velocity at several inter-rifting phases ───────────────────
tau_hats   = [0.1, 0.3, 0.5, 0.9]
phase_cols = [:royalblue, :darkorange, :forestgreen, :firebrick]
phase_labels = [L"\hat{\tau} = 0.1", L"\hat{\tau} = 0.3",
                L"\hat{\tau} = 0.5", L"\hat{\tau} = 0.9"]

v_css_norm(x̂, τ̂) =
    get_v_css(x̂ * a_d * T_d, τ̂ * T_d, b_d, vp_d, a_d) * T_d / b_d

# ── Panel (b): convergence of finite velocity sum to CSS at τ̂ = 0.5 ──────────
# Velocity (not displacement) converges: each past dike's velocity decays as 1/t.
τ̂_conv = 0.5

function v_finite_norm(x̂, τ̂, N)
    # v normalised by ω/T (= 1 in dimensionless units)
    sum(get_v_step(x̂, τ̂ + j, b_d, a_d) for j in 0:(N-1))
end

v_css_norm_conv(x̂, τ̂) = get_v_css(x̂, τ̂, b_d, vp_d, a_d)  # already normalised

Ns     = [1, 3, 10]
conv_cols   = [:royalblue, :darkorange, :forestgreen]
conv_labels = ["N = 1", "N = 3", "N = 10"]

mkpath("fig")

# ── Figure layout ─────────────────────────────────────────────────────────────
fig2 = Figure(size = (800, 1680), fontsize = 20)

ax21 = Axis(fig2[1, 1],
    xlabel = L"\hat{x} = x/(aT)",
    ylabel = L"v_{\mathrm{CSS}} \cdot T / \omega",
    title  = "(a)  CSS velocity — inter-rifting evolution",
    limits = (0, 5, 0, nothing))

for (τ̂, col, lbl) in zip(tau_hats, phase_cols, phase_labels)
    v_vals = v_css_norm.(x_hat, τ̂)
    lines!(ax21, x_hat, v_vals, color = col, linewidth = 2.5, label = lbl)
end

# Reference: time-averaged CSS velocity = ω/(aT²) → normalised = 1/aT = 1 here
# (not plotted, for context only)
axislegend(ax21, position = :rt, framevisible = true, labelsize = 18)

ax22 = Axis(fig2[2, 1],
    xlabel = L"\hat{x} = x/(aT)",
    ylabel = L"v \cdot T / \omega",
    title  = L"(b)  Convergence of velocity sum to CSS ($\hat{\tau} = 0.5$)",
    limits = (0, 5, 0, nothing))

for (N, col, lbl) in zip(Ns, conv_cols, conv_labels)
    v_vals = v_finite_norm.(x_hat, τ̂_conv, N)
    lines!(ax22, x_hat, v_vals, color = col, linewidth = 2.5, label = lbl)
end

v_css_vals = v_css_norm_conv.(x_hat, τ̂_conv)
lines!(ax22, x_hat, v_css_vals, color = :black, linewidth = 2.5,
       linestyle = :dash, label = "CSS  (N → ∞)")

# Far-field plate velocity reference: v·T/ω → 1 as x̂ → ∞
hlines!(ax22, [1.0], color = (:gray, 0.4), linewidth = 1.5, linestyle = :dot)
text!(ax22, 4.4, 1.03; text = L"v_{\mathrm{plate}}", fontsize = 15, color = :gray)

axislegend(ax22, position = :rb, framevisible = true, labelsize = 18)

save("fig/fig_css_behavior.png", fig2)
println("fig/fig_css_behavior.png  saved")