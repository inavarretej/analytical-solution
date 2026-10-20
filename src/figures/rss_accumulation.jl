using CairoMakie, SpecialFunctions, LaTeXStrings

include("../models/viscous.jl")

mkpath("fig")

# Cumulative EV RSS displacement over three cycles (after a figure by J.-A. Olive).
# Dimensionless units a = T = b = 1 (v = b/T = 1 in get_u_css). get_u_css(x, tau)
# is the displacement accumulated since the last intrusion (tau in [0, T)); each
# completed cycle adds b/2 per flank (Eq. conservation), so
#   u_cum(x, n T + tau) = n (b/2) sign(x) + u_RSS(x, tau).
a_d = 1.0; T_d = 1.0; b_d = 1.0; v_d = b_d / T_d
n_cycles = 3

function u_cum(x, t)
    n   = floor(Int, t / T_d + 1e-9)
    tau = t - n * T_d
    return n * (b_d / 2) * sign(x) + get_u_css(x, tau, b_d, v_d, a_d)
end

colors = cgrad(:viridis, 4 * n_cycles, categorical = true)

fig = Figure(size = (1200, 1400), fontsize = 26)

# (a) Profiles every T/4 over three cycles; thick lines at t = nT.
ax_a = Axis(fig[1, 1], xlabel = L"x/(aT)", ylabel = L"\mathrm{cumulative\ displacement}/b",
    limits = (-3, 3.45, nothing, nothing))
x_pos = collect(range(1e-4, 3.0, length = 600))
x_neg = -reverse(x_pos)
for k in 1:(4 * n_cycles)
    t = k * T_d / 4
    thick = k % 4 == 0
    col = colors[k]
    # t = nT is taken just before the next intrusion (tau -> T), where the
    # profile is uniform at n b/2 per flank.
    tt = thick ? t - 1e-6 : t
    lines!(ax_a, x_pos, u_cum.(x_pos, tt), color = col, linewidth = thick ? 4 : 2)
    lines!(ax_a, x_neg, u_cum.(x_neg, tt), color = col, linewidth = thick ? 4 : 2)
    lab = thick ? "$(k ÷ 4)T" : "$(k / 4)T"
    text!(ax_a, 3.03, u_cum(3.0, tt); text = lab, color = col, fontsize = 18, align = (:left, :center))
end
hlines!(ax_a, [0.0], color = (:gray, 0.5), linewidth = 1)
vlines!(ax_a, [0.0], color = (:gray, 0.5), linewidth = 1)

# (b) Displacement histories at fixed distances (x > 0 side).
ax_b = Axis(fig[2, 1], xlabel = L"t/T", ylabel = L"\mathrm{cumulative\ displacement}/b")
t_grid = collect(range(1e-4, n_cycles * T_d - 1e-6, length = 3000))
xs      = [1e-4, 0.1, 0.5, 1.0, 3.0]
x_labs  = [L"x = 0^+", L"x = 0.1\,aT", L"x = 0.5\,aT", L"x = aT", L"x = 3\,aT"]
x_cols  = [:black, :darkorange, :royalblue, :forestgreen, :mediumpurple]
for (x, lab, col) in zip(xs, x_labs, x_cols)
    lines!(ax_b, [0.0; t_grid], [0.0; u_cum.(x, t_grid)], color = col, linewidth = 3, label = lab)
end
lines!(ax_b, [0.0, n_cycles * T_d], [0.0, n_cycles * b_d / 2], color = :black,
    linestyle = :dot, linewidth = 2, label = L"\mathrm{mean\ rate}\ (b/2)\,t/T")
axislegend(ax_b, position = :lt, framevisible = true, labelsize = 22)

for (i, lab) in enumerate(("(a)", "(b)"))
    Label(fig[i, 1, TopLeft()], lab, fontsize = 30, padding = (0, 10, 10, 0), halign = :right)
end

save("fig/fig_rss_accumulation.png", fig, px_per_unit = 2)
println("Saved: fig/fig_rss_accumulation.png")
