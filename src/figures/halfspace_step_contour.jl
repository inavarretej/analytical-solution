include("halfspace_common.jl")

# STEP, interior panel as a filled |v| contour map instead of a quiver: same
# model, same profile column, but the interior field is evaluated on a dense
# structured grid (nx_dense x nz_dense in halfspace_common.jl) rather than
# the sparse arrow grid. `orientation` lays phases out top-to-bottom
# (:vertical) or left-to-right (:horizontal); `fontsize` is the Figure's base
# text size (title/legend sizes scale with it, see plot_model_rows! in
# halfspace_common.jl). `zoom_extent = (x_max_m, z_max_m)` adds a third panel
# zoomed into the near-dike peak (x in [-x_max_m,x_max_m], z in [0,z_max_m]),
# same color scale as the main contour panel; `nothing` omits it.
function generate_step_contour_figure(; orientation::Symbol = :vertical, fontsize::Real = 30,
        zoom_extent = (10e3, 10e3), tau_R_mult_step = [0.001, 0.25, 0.5, 1.0, 2.0, 3.0])
    # Own phase list (last row 3 tau_R instead of the shared 4.9 tau_R, where
    # the near-dike zoom falls entirely below color_min_step).
    taus   = tau_R_mult_step .* tau_R
    labels = ["$(m)·τ_R" for m in tau_R_mult_step]

    # tau_R = eta_a/mu_A is the Maxwell time set by the half-space viscosity and
    # the elastic layer's modulus (shared-modulus closure) -- independent of the
    # CSS recurrence interval T, so STEP no longer borrows CSS's v_plate scale.
    vel_nondim_step  = tau_R / b_nom
    disp_nondim_step = 1.0 / b_nom

    make_V_step(tau_i) = xp -> get_v_step(xp, tau_i, b_nom, a_nom)
    make_U_step(tau_i) = xp -> get_u_step(xp, tau_i, b_nom, a_nom)

    # Fixed range, narrowed to the two decades where the bulk of the data
    # actually sits (see the per-phase ranges this script prints).
    color_min_step = 0.01
    color_max_step = 1.0

    # Width/height scale with the panel count (profile + interior [+ zoom]):
    # the base 2000/430 figures were tuned for 2 panel columns/rows.
    n_panels = zoom_extent === nothing ? 2 : 3
    panel_scale = n_panels / 2
    fig_size = orientation == :vertical ?
        (2000 * panel_scale, 430 * n_rows) : (430 * n_rows, 2000 * panel_scale)
    fig_step = Figure(size = fig_size, fontsize = fontsize)
    println("STEP model (contour):")
    plot_model_rows!(fig_step, taus, labels, make_V_step, make_U_step;
        orientation = orientation, fontsize = fontsize, interior_panel = :contour,
        zoom_extent = zoom_extent,
        title = "STEP model", colorbar_label = "|v| / (b/τ_R)  (dimensionless)",
        vel_nondim = vel_nondim_step, vel_symbol = "(b/τ_R)",
        disp_nondim = disp_nondim_step, disp_symbol = "b", u_lim = 2.0, vel_lim = 2.0,
        length_scale = H, length_symbol = "H", guide_formula = nothing,
        color_min = color_min_step, color_max = color_max_step, ref_speed = 0.01, colorscale = log10)
    out_path = "fig/fig_halfspace_step_contour$(orientation_suffix(orientation)).png"
    save(out_path, fig_step, px_per_unit = 2)
    println("Saved: $out_path")

    return fig_step
end

generate_step_contour_figure()
