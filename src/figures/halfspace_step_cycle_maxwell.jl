using CairoMakie, QuadGK, SpecialFunctions

include("../models/viscous.jl")
include("../models/maxwell.jl")
include("../parameters.jl")
include("../halfspace.jl")

mkpath("fig")

# Same illustrative dike parameters as fig_halfspace_step_cycle.jl (Natron
# best-fit a, b), for a directly comparable companion figure.
a_nom = 1.18e-4     # m/s
b_nom = 0.662       # m

# tau_R is NOT redefined here: it comes straight from src/parameters.jl as
# tau_R = eta_a / mu_A, i.e. viscosity over the mantle's own elastic shear
# modulus -- the Maxwell relaxation time itself, independent of a_nom and
# not the shared-modulus closure of Section 2.6.3 (Eq. tau_R_shared).
println("tau_R (= eta_a / mu_A) = $(tau_R / yr) yr")

# Maxwell length (paper Eq. 51): ell_M := a * tau_R. Fixed (not per-column):
# the purely elastic length scale that governs the Maxwell regimes, marked
# as a vertical guide in every panel instead of the viscous L = a*tau.
ell_M = a_nom * tau_R
println("ell_M = a*tau_R = $(ell_M / 1e3) km = $(ell_M / H) H")

taus_yr = [0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 100.0]
taus    = taus_yr .* yr
labels  = ["$(t) yr" for t in taus_yr]

x_span_m = 600e3    # +/- 600 km (widened for the tau=100yr column, kept consistent with the viscous companion figure)
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

# Local, capped-maxevals copies of the halfspace angular convolution (NOT
# src/halfspace.jl -- that shared, uncapped version is still what the viscous
# figures use, where it is fine because get_v_step/get_v_css are closed-form
# and cheap). Wrapping the Maxwell velocities (which themselves need an
# internal k-integral) inside an UNCAPPED adaptive angular quadrature causes
# runaway refinement -- confirmed empirically: a single halfspace_vx call
# hung past several minutes with the default (unbounded) maxevals, while
# maxevals=500 converges every tested grid corner in <3s with an error
# estimate many orders below what a qualitative figure needs.
function halfspace_vx_capped(x, z, V::Function)
    integrand(t) = (2 / pi) * sin(t)^2 * V(x - z * tan(t))
    val, _ = quadgk(integrand, -pi/2, pi/2; maxevals = 500)
    return val
end

function halfspace_vz_capped(x, z, V::Function)
    integrand(t) = (2 / pi) * sin(t) * cos(t) * V(x - z * tan(t))
    val, _ = quadgk(integrand, -pi/2, pi/2; maxevals = 500)
    return val
end

# Compute the surface profile and interior (us, ws) field ONCE per block
# (previously compute_speed_max and plot_model_block! each recomputed the
# same grid independently -- 2x redundant work). Each (grid point, tau) pair
# is an independent quadrature, so flatten and thread across all of them.
function compute_block_fields(taus, make_V, make_U, nondim)
    n_tau, n_pts, n_prof = length(taus), length(xs), length(x_prof)
    Vs = [make_V(tau_i) for tau_i in taus]
    Us = [make_U(tau_i) for tau_i in taus]

    # Displacement u nondimensionalized by b_nom itself (u/b, same convention
    # as fig_maxwell_regimes.jl / fig_maxwell_vs_viscous.jl and the viscous
    # companion figure), NOT by nondim_maxwell -- that factor scales velocity
    # by tau_R/b, an unrelated combination for a quantity already O(b).
    v_prof_mat = Matrix{Float64}(undef, n_prof, n_tau)
    u_prof_mat = Matrix{Float64}(undef, n_prof, n_tau)
    Threads.@threads for idx in 1:(n_tau * n_prof)
        j = (idx - 1) ÷ n_prof + 1
        i = (idx - 1) % n_prof + 1
        v_prof_mat[i, j] = Vs[j](x_prof[i]) * nondim
        u_prof_mat[i, j] = Us[j](x_prof[i]) / b_nom
    end

    us_mat = Matrix{Float64}(undef, n_pts, n_tau)
    ws_mat = Matrix{Float64}(undef, n_pts, n_tau)
    Threads.@threads for idx in 1:(n_tau * n_pts)
        j = (idx - 1) ÷ n_pts + 1
        i = (idx - 1) % n_pts + 1
        us_mat[i, j] = halfspace_vx_capped(xs[i], zs[i], Vs[j]) * nondim
        ws_mat[i, j] = halfspace_vz_capped(xs[i], zs[i], Vs[j]) * nondim
    end

    all_v_prof = [v_prof_mat[:, j] for j in 1:n_tau]
    all_u_prof = [u_prof_mat[:, j] for j in 1:n_tau]
    all_us     = [us_mat[keep, j] for j in 1:n_tau]
    all_ws     = [ws_mat[keep, j] for j in 1:n_tau]
    return all_v_prof, all_u_prof, all_us, all_ws
end

function plot_model_block!(fig, r1, r2, labels, all_v_prof, all_u_prof, all_us, all_ws, label_suffix; color_max)
    speeds    = [sqrt.(all_us[col] .^ 2 .+ all_ws[col] .^ 2) for col in eachindex(labels)]
    speed_min = minimum(minimum.(speeds))   # log color scale needs a positive lower bound, not 0
    # Arrow length uses THIS block's own max (post-origin-removal) so STEP and
    # CSS are each well populated visually; color still uses the shared
    # color_max so the color mapping is comparable across both blocks.
    len_scale = max_len_H / maximum(maximum.(speeds))
    # Symmetric right-axis limit for u/b, fixed across this block's columns
    # (same idea as the viscous companion figure).
    u_lim = 1.05 * maximum(maximum.(abs.(u) for u in all_u_prof))

    for (col, lab) in enumerate(labels)
        v_prof_nondim = all_v_prof[col]
        u_prof_nondim = all_u_prof[col]

        # Row 1: surface velocity profile v_x(x) for this column's tau, with
        # the fixed Maxwell length ell_M marked as a vertical guide (unlike
        # the viscous L = a*tau, ell_M does not depend on tau). Displacement
        # u_x(x) shares the panel via a twin axis on the right.
        ax1 = Axis(fig[r1, col],
            xlabel = "x / H", ylabel = "v_x $label_suffix  (dimensionless)",
            title  = "τ = $lab",
            limits = ((-x_span_m / H, x_span_m / H), (-1.0, 1.0)))
        hlines!(ax1, [0.0], color = (:gray, 0.4), linestyle = :dot)
        vlines!(ax1, [0.0], color = (:gray, 0.4), linestyle = :dot)
        vlines!(ax1, [ell_M / H], color = (:firebrick, 0.6), linestyle = :dash, linewidth = 1.5)
        l_v = lines!(ax1, x_prof_H, v_prof_nondim, color = :forestgreen, linewidth = 2.5)

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

make_V_step_maxwell(tau_i) = xp -> get_v_relax(xp, tau_i, b_nom, a_nom, tau_R)
make_V_css_maxwell(tau_i)  = xp -> get_v_css_maxwell_closed(xp, tau_i, b_nom, v_plate_css, a_nom, tau_R; k_limit_scale = 1.0)
make_U_step_maxwell(tau_i) = xp -> get_u_viscoelastic(xp, tau_i, b_nom, a_nom, tau_R)
make_U_css_maxwell(tau_i)  = xp -> get_u_css_maxwell_closed(xp, tau_i, b_nom, v_plate_css, a_nom, tau_R; k_limit_scale = 1.0)

# Both blocks nondimensionalized by tau_R/b (v * tau_R / b), per the Maxwell
# relaxation time -- NOT by 1/V as in the viscous companion figure.
nondim_maxwell = tau_R / b_nom

println("Computing Maxwell STEP fields...")
@time step_v_prof, step_u_prof, step_us, step_ws = compute_block_fields(taus, make_V_step_maxwell, make_U_step_maxwell, nondim_maxwell)

println("Computing Maxwell CSS fields...")
@time css_v_prof, css_u_prof, css_us, css_ws = compute_block_fields(taus, make_V_css_maxwell, make_U_css_maxwell, nondim_maxwell)

# Shared color_max across STEP and CSS, so the color mapping is comparable
# between blocks. Arrow length is scaled per block (see plot_model_block!)
# so both blocks stay well populated.
speed_max_step_only = maximum(maximum.(sqrt.(u .^ 2 .+ w .^ 2) for (u, w) in zip(step_us, step_ws)))
speed_max_css_only  = maximum(maximum.(sqrt.(u .^ 2 .+ w .^ 2) for (u, w) in zip(css_us, css_ws)))
color_max = max(speed_max_step_only, speed_max_css_only)

fig = Figure(size = (5000, 1800), fontsize = 22)

# ── Maxwell STEP model (rows 1-2)
println("Maxwell STEP model:")
plot_model_block!(fig, 1, 2, labels, step_v_prof, step_u_prof, step_us, step_ws, "· τ_R / b"; color_max = color_max)
Colorbar(fig[2, length(taus) + 1], limits = (0.0, color_max), colormap = :viridis, label = "|v| · τ_R / b  (dimensionless)")
Label(fig[1, 1:length(taus), Top()], "Maxwell STEP model", fontsize = 26, padding = (0, 0, 30, 0))

# ── Maxwell CSS model (rows 3-4)
println("\nMaxwell CSS model (T = 23 yr):")
plot_model_block!(fig, 3, 4, labels, css_v_prof, css_u_prof, css_us, css_ws, "· τ_R / b"; color_max = color_max)
Colorbar(fig[4, length(taus) + 1], limits = (0.0, color_max), colormap = :viridis, label = "|v| · τ_R / b  (dimensionless)")
Label(fig[3, 1:length(taus), Top()], "Maxwell CSS model (T = 23 yr)", fontsize = 26, padding = (0, 0, 30, 0))

save("fig/fig_halfspace_step_cycle_maxwell.png", fig, px_per_unit = 2)
println("Saved: fig/fig_halfspace_step_cycle_maxwell.png")
