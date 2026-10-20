using CairoMakie, LaTeXStrings, SpecialFunctions

mkpath("fig")

# Raw (not peak-normalised) counterpart of step_selfsimilar.jl: the same
# similarity variables, but each curve keeps its own prefactor, so absolute
# amplitudes are compared as well as the tail decay.
# EV:      u^STEP/b        = (1/pi) atan(1/xi),                 xi   := x/(at)
#          t*v^STEP/b      = (1/pi) xi/(1+xi^2)
# Foulger: u_Foulger/b     = (1/2) erfc(xi_F/2),                xi_F := x/sqrt(kappa*t)
#          t*v_Foulger/b   = xi_F/(4 sqrt(pi)) exp(-xi_F^2/4)
# (Foulger: same initial condition b*H(-x) under Fickian diffusion, so
#  u(0,t) = b/2 as in EV; velocity is t * d/dt of the displacement.)
# Panel (a): displacement (top); panel (b): velocity (bottom).
xi_tail = 10 .^ range(log10(0.05), log10(20.0), length = 400)
u_ev(x)      = atan(1 / x) / pi
u_foulger(x) = erfc(x / 2) / 2
v_ev(x)      = x / (1 + x^2) / pi
v_foulger(x) = x / (4 * sqrt(pi)) * exp(-x^2 / 4)

println("u_EV(0+)      = 1/2       = ", 0.5)
println("u_Foulger(0+) = 1/2       = ", 0.5)
println("max v_EV      = 1/(2pi)   = ", 1 / (2pi), "  at xi = 1")
println("max v_Foulger = e^-1/2/(2 sqrt(2pi)) = ", exp(-0.5) / (2 * sqrt(2pi)), "  at xi_F = sqrt(2)")

fig = Figure(size = (900, 1200), fontsize = 33.75)
guide_x = [2.0, 20.0]

# ── Panel (a): displacement u/b ──────────────────────────────────────────────
ax_u = Axis(fig[1, 1],
    xlabel = L"\xi = x/(at)", ylabel = L"u/b",
    xscale = log10, yscale = log10,
    limits = (0.05, 20, 1e-4, 2))
lines!(ax_u, xi_tail, u_ev.(xi_tail), color = :forestgreen, linewidth = 3.5,
    label = L"EV:  $\frac{1}{\pi}\arctan(1/\xi)$")
lines!(ax_u, xi_tail, u_foulger.(xi_tail), color = :crimson, linewidth = 3.5,
    label = L"Foulger:  $\frac{1}{2}\mathrm{erfc}(\xi_F/2)$")
lines!(ax_u, guide_x, u_ev(2.0) .* (2.0 ./ guide_x),
    color = (:gray, 0.6), linestyle = :dot, linewidth = 2.0)
text!(ax_u, 6.0, 0.04; text = L"\propto \xi^{-1}", fontsize = 24.55, color = :gray)
axislegend(ax_u, position = :lb, framevisible = false, labelsize = 27.61)

# ── Panel (b): velocity t*v/b ────────────────────────────────────────────────
ax_v = Axis(fig[2, 1],
    xlabel = L"\xi = x/(at)", ylabel = L"t\,|v|/b",
    xscale = log10, yscale = log10,
    limits = (0.05, 20, 1e-4, 1))
lines!(ax_v, xi_tail, v_ev.(xi_tail), color = :forestgreen, linewidth = 3.5,
    label = L"EV:  $\frac{1}{\pi}\,\xi/(1+\xi^2)$")
lines!(ax_v, xi_tail, v_foulger.(xi_tail), color = :crimson, linewidth = 3.5,
    label = L"Foulger:  $\frac{\xi_F}{4\sqrt{\pi}}\,e^{-\xi_F^2/4}$")
lines!(ax_v, guide_x, v_ev(2.0) .* (2.0 ./ guide_x),
    color = (:gray, 0.6), linestyle = :dot, linewidth = 2.0)
text!(ax_v, 6.0, 0.03; text = L"\propto \xi^{-1}", fontsize = 24.55, color = :gray)
axislegend(ax_v, position = :lb, framevisible = false, labelsize = 27.61)

for (k, lab) in enumerate(["(a)", "(b)"])
    Label(fig[k, 1, TopLeft()], lab, fontsize = 36, font = :bold,
        padding = (0, 10, 25, 0), halign = :right)
end

save("fig/fig_step_selfsimilar_raw.png", fig, px_per_unit = 2)
println("Saved: fig/fig_step_selfsimilar_raw.png")
