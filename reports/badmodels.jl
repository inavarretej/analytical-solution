using CairoMakie

# --- 1. Define Parameters ---
V = 10.0  # Extension velocity
H = 15.0  # Depth H
x = range(-50, 50, length=500)
x_safe = [val == 0.0 ? 1e-6 : val for val in x]

# --- 2. Calculate Velocities ---
v_dislocation = (V / π) .* atan.(x_safe ./ H)
v_dike        = (V / π) .* atan.(H ./ x_safe)

# --- 3. Blue plot: dislocation ---
fig1 = Figure(size = (800, 500), font = "Arial")
ax1 = Axis(fig1[1, 1],
    xlabel = "x / H",
    ylabel = "v / V",
    xgridvisible = true,
    ygridvisible = true,
    xgridstyle = :dash,
    ygridstyle = :dash,
    xlabelsize = 22,
    ylabelsize = 22

)
hlines!(ax1, [0], color = :black, linestyle = :dot)
lines!(ax1, x_safe, v_dislocation, color = :blue, linewidth = 3,
       label = "Dislocation (Locked to H)")
save("fig/svg/velocity_profile_dislocation.svg", fig1)

# --- 4. Red plot: continuous dike ---
fig2 = Figure(size = (800, 500), font = "Arial")
ax2 = Axis(fig2[1, 1],
    xlabel = "x / H",
    ylabel = "v / V",
    xgridvisible = true,
    ygridvisible = true,
    xgridstyle = :dash,
    ygridstyle = :dash
)
hlines!(ax2, [0], color = :black, linestyle = :dot)
lines!(ax2, x_safe, v_dike, color = :red, linewidth = 3,
       linestyle = :dash, label = "Continuous Dike (Depth H)")
save("fig/svg/velocity_profile_lafemina.svg", fig2)




