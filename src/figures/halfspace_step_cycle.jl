using CairoMakie, QuadGK, SpecialFunctions

include("../models/viscous.jl")
include("../models/maxwell.jl")
include("../parameters.jl")
include("../halfspace.jl")

mkpath("fig")

a_nom = 1.18e-4     # m/s
b_nom = 0.662       # m


taus_yr = [0.5, 2.0, 5.0, 10.0, 20.0, 23.0]
taus    = taus_yr .* yr
labels  = ["$(t) yr" for t in taus_yr]

x_span_m = 600e3    # +/- 600 km (widened for the tau=100yr column: L=a*tau ~ 372km there)
z_top_m  = 20e3     # 20 km: shallow margin below the surface (also sets the room available for arrow length)
z_bot_m  = 300e3    # 300 km depth
x_prof = collect(range(-x_span_m, x_span_m, length = 400))
x_prof_H = x_prof ./ H

# Interior quiver keeps its own, narrower x-range (independent of x_span_m
# above): widening it in step with the surface-profile row would have halved
# the on-screen arrow length (DataAspect ties x and z to the same physical
# scale, and z_top_m/nz -- unrelated to x_span_m -- already caps arrow length
# to avoid overlapping the row above/below it).
x_quiver_span_m = 300e3

nx, nz = 15, 9
x_grid = collect(range(-x_quiver_span_m, x_quiver_span_m, length = nx))
z_grid = collect(range(z_top_m, z_bot_m, length = nz))
xs = Float64[]; zs = Float64[]
for zi in z_grid, xi in x_grid
    push!(xs, xi); push!(zs, zi)
end
x_grid_H = x_grid ./ H
z_grid_H = z_grid ./ H

# Exclude the grid point sitting right above the dike (x=0, shallowest z):
# the near-source field there is huge compared to the rest of the domain and
# was dominating both the color scale and the arrow-length scale.
origin_mask = isapprox.(xs, 0.0; atol = 1.0) .& isapprox.(zs, z_top_m; atol = 1.0)
keep  = .!origin_mask
xs_H  = xs[keep] ./ H
zs_H  = zs[keep] ./ H

dx_H      = x_grid_H[2] - x_grid_H[1]
margin_H  = z_grid_H[1]
max_len_H = min(2.2 * dx_H, 1.7 * margin_H)

T_css       = 23.0 * yr
v_plate_css = b_nom / T_css

function compute_speed_max(taus, make_V, nondim)
    speed_max = 0.0
    for tau_i in taus
        V = make_V(tau_i)
        us = halfspace_vx.(xs, zs, Ref(V)) .* nondim
        ws = halfspace_vz.(xs, zs, Ref(V)) .* nondim
        speed_max = max(speed_max, maximum(sqrt.(us .^ 2 .+ ws .^ 2)[keep]))
    end
    return speed_max
end

function plot_model_block!(fig, r1, r2, taus, labels, make_V, make_U, nondim, label_suffix; color_max)
    # Surface velocity profile V(x) for each tau column, nondimensionalized.
    all_v_prof = [make_V(tau_i).(x_prof) .* nondim for tau_i in taus]

    # Surface displacement profile u(x), nondimensionalized by b_nom itself
    # (u is already an O(b) quantity, so u/b_nom is O(1) without any extra
    # tau-dependent rescaling). Fixed symmetric limit (-2, 2), shared by both
    # the STEP (row 1) and CSS (row 3) blocks, so the two are directly
    # comparable rather than each being auto-scaled to its own data.
    all_u_prof = [make_U(tau_i).(x_prof) ./ b_nom for tau_i in taus]
    u_lim      = 2.0

    # Interior velocity field (vx, vz) on the (x_grid, z_grid) mesh, one
    # (us, ws) pair per tau column, via the halfspace convolution (origin
    # point excluded).
    all_us = Vector{Vector{Float64}}(undef, length(taus))
    all_ws = Vector{Vector{Float64}}(undef, length(taus))
    for (col, tau_i) in enumerate(taus)
        V = make_V(tau_i)
        all_us[col] = (halfspace_vx.(xs, zs, Ref(V)) .* nondim)[keep]
        all_ws[col] = (halfspace_vz.(xs, zs, Ref(V)) .* nondim)[keep]
    end
    speeds    = [sqrt.(all_us[col] .^ 2 .+ all_ws[col] .^ 2) for col in eachindex(taus)]
    speed_min = minimum(minimum.(speeds))   # log color scale needs a positive lower bound, not 0
    # Arrow length uses THIS block's own max (post-origin-removal) so STEP and
    # CSS are each well populated visually; color still uses the shared
    # color_max so the color mapping is comparable across both blocks.
    len_scale = max_len_H / maximum(maximum.(speeds))

    for (col, (tau_i, lab)) in enumerate(zip(taus, labels))
        L = a_nom * tau_i
        v_prof_nondim = all_v_prof[col]
        u_prof_nondim = all_u_prof[col]

        # Row 1: surface velocity profile v_x(x) for this column's tau, with
        # the characteristic length L = a*tau marked as a vertical guide.
        ax1 = Axis(fig[r1, col],
            xlabel = "x / H", ylabel = "v_x $label_suffix  (dimensionless)",
            title  = "τ = $lab   (L/H = a·τ/H = $(round(L / H, sigdigits = 3)))",
            limits = ((-x_span_m / H, x_span_m / H), (-5.0, 5.0)))
        hlines!(ax1, [0.0], color = (:gray, 0.4), linestyle = :dot)
        hlines!(ax1, [1.0], color = (:gray, 0.4), linestyle = :dot)
        hlines!(ax1, [-1.0], color = (:gray, 0.4), linestyle = :dot)
        vlines!(ax1, [0.0], color = (:gray, 0.4), linestyle = :dot)
        vlines!(ax1, [L / H], color = (:firebrick, 0.6), linestyle = :dash, linewidth = 1.5)
        l_v = lines!(ax1, x_prof_H, v_prof_nondim, color = :forestgreen, linewidth = 2.5)

        # Displacement u(x) shares the panel via a twin axis on the right,
        # linked in x but with its own (symmetric) y-scale since u and v have
        # unrelated magnitudes and units.
        ax1b = Axis(fig[r1, col],
            ylabel = "u_x / b  (dimensionless)",
            yaxisposition = :right,
            limits = ((-x_span_m / H, x_span_m / H), (-u_lim, u_lim)))
        hidespines!(ax1b)
        hidexdecorations!(ax1b)
        linkxaxes!(ax1, ax1b)
        l_u = lines!(ax1b, x_prof_H, u_prof_nondim, color = :darkorange, linewidth = 2.5, linestyle = :dash)
        if col == 1
            axislegend(ax1, [l_v, l_u], ["v_x $label_suffix", "u_x / b"],
                position = :lt, framevisible = false, labelsize = 16)
        end

        us, ws = all_us[col], all_ws[col]
        speed  = speeds[col]
        println("τ = $lab:  |v| $label_suffix range = [$(minimum(speed)), $(maximum(speed))] (dimensionless)")

        # Row 2: (x, depth) panel for the interior velocity quiver.
        ax2 = Axis(fig[r2, col],
            xlabel = "x / H", ylabel = "z / H",
            yreversed = true, aspect = DataAspect(),
            limits = ((x_grid_H[1], x_grid_H[end]), (0.0, z_grid_H[end])))

        arrows!(ax2, xs_H, zs_H, us .* len_scale, ws .* len_scale,
            arrowsize   = 7,
            linewidth   = 2.5,
            lengthscale = 1.0,
            color       = speed, colorrange = (0.0, color_max), colormap = :viridis)
    end
    return speed_min, color_max
end

make_V_step(tau_i) = xp -> get_v_step(xp, tau_i, b_nom, a_nom)
make_V_css(tau_i)  = xp -> get_v_css(xp, tau_i, b_nom, v_plate_css, a_nom)
make_U_step(tau_i) = xp -> get_u_step(xp, tau_i, b_nom, a_nom)
make_U_css(tau_i)  = xp -> get_u_css(xp, tau_i, b_nom, v_plate_css, a_nom)

# Shared color_max across STEP and CSS, computed after excluding the origin
# point, so the color mapping is comparable between blocks. Arrow length is
# scaled per block (see plot_model_block!) so both blocks stay well populated.
speed_max_step_only = compute_speed_max(taus, make_V_step, T_css / b_nom)
speed_max_css_only  = compute_speed_max(taus, make_V_css, T_css / b_nom)
color_max = max(speed_max_step_only, speed_max_css_only)

fig = Figure(size = (5000, 1800), fontsize = 22)

# ── STEP model (rows 1-2)
println("STEP model:")
plot_model_block!(fig, 1, 2, taus, labels, make_V_step, make_U_step, T_css / b_nom, "/ V"; color_max = color_max)
Colorbar(fig[2, length(taus) + 1], limits = (0.0, color_max), colormap = :viridis, label = "|v| / V  (dimensionless)")
Label(fig[1, 1:length(taus), Top()], "STEP model", fontsize = 26, padding = (0, 0, 30, 0))

# ── CSS model (rows 3-4)
v_plate_css_nondim = v_plate_css * (T_css / b_nom)
println("\nCSS model (T = 23 yr):  v_plate (nondim) = $v_plate_css_nondim")
plot_model_block!(fig, 3, 4, taus, labels, make_V_css, make_U_css, T_css / b_nom, "/ V"; color_max = color_max)
Colorbar(fig[4, length(taus) + 1], limits = (0.0, color_max), colormap = :viridis, label = "|v| / V  (dimensionless)")
Label(fig[3, 1:length(taus), Top()], "CSS model (T = 23 yr)", fontsize = 26, padding = (0, 0, 30, 0))

save("fig/fig_halfspace_step_cycle.png", fig, px_per_unit = 2)
println("Saved: fig/fig_halfspace_step_cycle.png")
