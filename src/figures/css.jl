using CairoMakie, SpecialFunctions, QuadGK
using LaTeXStrings

include("../models/viscous.jl")
include("../models/maxwell.jl")

# ── Dimensionless units: a = T = b = v_plate = 1 ──────────────────────────────
# x̂ = x/(aT),  τ̂ = τ/T ∈ [0,1),  u normalised by b,  v normalised by b/T
a_d  = 1.0;  T_d = 1.0;  b_d = 1.0;  vp_d = b_d / T_d

x_hat = collect(range(0, 5.0, length = 400))

# The RSS-velocity-vs-phase panel (locked width at aT across tau_hat) was
# dropped: that transition is already shown in Figure~\ref{fig:viscous_halfspace}
# (fig_halfspace_step_cycle.png), so this figure now focuses solely on the
# convergence-with-N-cycles result below.

# Convergence of finite velocity sum to CSS at τ̂ = 0.5.
# Velocity (not displacement) converges: each past dike's velocity decays as 1/t.
τ̂_conv = 0.5

function v_finite_norm(x̂, τ̂, N)
    # v normalised by ω/T (= 1 in dimensionless units)
    sum(get_v_step(x̂, τ̂ + j, b_d, a_d) for j in 0:(N-1))
end

v_css_norm_conv(x̂, τ̂) = get_v_css(x̂, τ̂, b_d, vp_d, a_d)  # already normalised

Ns     = [1, 3, 10, 100]
conv_cols   = [:royalblue, :darkorange, :forestgreen, :red]
conv_labels = ["N = 1", "N = 3", "N = 10", "N = 100"]

mkpath("fig")

# ── Figure layout ─────────────────────────────────────────────────────────────
fig2 = Figure(size = (800, 840), fontsize = 30)

ax22 = Axis(fig2[1, 1],
    xlabel = L"x/(aT)",
    ylabel = L"v \cdot T / b",
    title  = L"Convergence of velocity sum to RSS ($\tau/T = 0.5$)",
    limits = (0, 5, 0, nothing))

for (N, col, lbl) in zip(Ns, conv_cols, conv_labels)
    v_vals = v_finite_norm.(x_hat, τ̂_conv, N)
    lines!(ax22, x_hat, v_vals, color = col, linewidth = 2.5, label = lbl)
end

v_css_vals = v_css_norm_conv.(x_hat, τ̂_conv)
lines!(ax22, x_hat, v_css_vals, color = :black, linewidth = 2.5,
       linestyle = :dash, label = "RSS  (N → ∞)")

# Far-field plate velocity reference: v·T/ω → 1/2 as x̂ → ∞ (each flank
# carries half the plate rate; see vault note "Factor-of-2 error in the EV
# RSS velocity formula").
hlines!(ax22, [0.5], color = (:gray, 0.4), linewidth = 1.5, linestyle = :dot)
text!(ax22, 4.4, 0.53; text = L"V/2", fontsize = 22.5, color = :gray)

axislegend(ax22, position = :rb, framevisible = true, labelsize = 27)

save("fig/fig_css_behavior.png", fig2)
println("fig/fig_css_behavior.png  saved")