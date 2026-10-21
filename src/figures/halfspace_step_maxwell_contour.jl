include("halfspace_common.jl")

# EVE (Maxwell) STEP, the counterpart of halfspace_step_contour.jl (Fig. 3):
# same layout, phases, dense grid (nx_dense x nz_dense), 1H zoom, color range
# and plotting routine (plot_model_rows!). Only the model changes, and the
# interior kernels are the maxevals-capped versions: the uncapped
# halfspace_vx/vz (src/halfspace.jl) hang on the nested Maxwell quadrature.
#
# Expensive: ~0.2 s per kernel evaluation, i.e. ~0.4 s per grid point, over
# 2 x 400 x 240 points x 6 phases. Run with threads, e.g.
#   julia -t auto --project=. src/figures/halfspace_step_maxwell_contour.jl
# (~11 h with 12 threads).

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

# First row: tau = 0+ via the closed form get_v_relax_0plus (maxwell.jl),
# because get_v_relax at tau = 0.001 tau_R costs ~7 s per interior-kernel
# evaluation (~100 h for that row alone). At 0+ the displacement is zero for
# x != 0 (the whole opening sits at the dike). Documentation: Obsidian note
# 30-models/analytical/EVE contour figures - tau = 0+ limit.md
v_relax_0plus(x) = get_v_relax_0plus(x, b_nom, a_nom, tau_R)

function generate_step_maxwell_contour_figure(; orientation::Symbol = :vertical, fontsize::Real = 30,
        zoom_extent = (10e3, 10e3), tau_R_mult_step = [0.0, 0.25, 0.5, 1.0, 2.0, 3.0])
    # Same phases as the EV STEP figure (halfspace_step_contour.jl), except the
    # first row: tau = 0+ (closed form above) instead of 0.001 tau_R.
    taus   = tau_R_mult_step .* tau_R
    labels = [m == 0 ? "0⁺" : "$(m)·τ_R" for m in tau_R_mult_step]

    vel_nondim_step  = tau_R / b_nom
    disp_nondim_step = 1.0 / b_nom

    make_V_step(tau_i) = tau_i == 0 ? v_relax_0plus :
        (xp -> get_v_relax(xp, tau_i, b_nom, a_nom, tau_R))
    make_U_step(tau_i) = tau_i == 0 ? (xp -> 0.0) :
        (xp -> get_u_viscoelastic(xp, tau_i, b_nom, a_nom, tau_R))

    color_min_step = 0.01
    color_max_step = 1.0

    n_panels = zoom_extent === nothing ? 2 : 3
    panel_scale = n_panels / 2
    fig_size = orientation == :vertical ?
        (2000 * panel_scale, 430 * n_rows) : (430 * n_rows, 2000 * panel_scale)
    fig_step = Figure(size = fig_size, fontsize = fontsize)
    println("EVE STEP model (contour), $(Threads.nthreads()) threads:")
    plot_model_rows!(fig_step, taus, labels, make_V_step, make_U_step;
        orientation = orientation, fontsize = fontsize, interior_panel = :contour,
        zoom_extent = zoom_extent,
        title = "EVE STEP model", colorbar_label = "|v| / (b/τ_R)  (dimensionless)",
        vel_nondim = vel_nondim_step, vel_symbol = "(b/τ_R)",
        disp_nondim = disp_nondim_step, disp_symbol = "b", u_lim = 2.0, vel_lim = 2.0,
        length_scale = H, length_symbol = "H", guide_formula = nothing,
        time_symbol = "t",
        color_min = color_min_step, color_max = color_max_step, ref_speed = 0.01, colorscale = log10,
        hvx = halfspace_vx_capped, hvz = halfspace_vz_capped,  # dashed guide at a*tau, as in Fig. 3
        cache_prefix = "step_maxwell_contour")  # resume from cache/ after an interruption
    out_path = "fig/fig_halfspace_step_maxwell_contour$(orientation_suffix(orientation)).png"
    save(out_path, fig_step, px_per_unit = 2)
    println("Saved: $out_path")

    return fig_step
end

@time generate_step_maxwell_contour_figure()
