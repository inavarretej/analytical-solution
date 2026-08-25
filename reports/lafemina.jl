using CairoMakie

# --- 1. Define EVZ Parameters ---
# These parameters represent the EVZ dislocation model from the paper
u_evz = 13.0  # Opening rate (mm/yr)
d_evz = 5.0   # Locking depth (km)

# Symmetric x-axis from -50 to 50 km
x = range(-50, 50, length=500)
# Standard safety check for zero-division, though atan handle 0 gracefully
x_safe = [val == 0.0 ? 1e-6 : val for val in x]

# --- 2. Calculate EVZ Velocity ---
# The surface velocity is calculated using the atan dislocation formula
# Formula: v(x) = (u / π) * atan(x / d)
v_evz = (u_evz / π) .* atan.(x_safe ./ d_evz)

# --- 3. Plotting with CairoMakie ---
fig = Figure(size = (800, 500), font = "Arial")
ax = Axis(fig[1, 1],
    xlabel = "x / H",
    ylabel = "v / V",
    xgridstyle = :dash,
    ygridstyle = :dash,

    xlabelsize = 22,
    ylabelsize = 22
)

# Reference lines for clarity
hlines!(ax, [0], color = :black, linestyle = :dot)

# Plotting the EVZ contribution
lines!(ax, x, v_evz, color = :red, linewidth = 3, 
       )

save("fig/svg/velocity_profile_lafemina.svg", fig) # Uncomment to save as vector PDF


# Display the plot
fig