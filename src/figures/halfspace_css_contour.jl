include("halfspace_common.jl")

# CSS/RSS, interior panel as a filled |v| contour map instead of a quiver:
# same model, same profile column, but the interior field is evaluated on a
# dense structured grid (nx_dense x nz_dense in halfspace_common.jl) rather
# than the sparse arrow grid. `orientation` lays phases out top-to-bottom
# (:vertical) or left-to-right (:horizontal); `fontsize` is the Figure's base
# text size (title/legend sizes scale with it, see plot_model_rows! in
# halfspace_common.jl). `zoom_extent = (x_max_m, z_max_m)` adds a third panel
# zoomed into the near-dike peak (x in [-x_max_m,x_max_m], z in [0,z_max_m]),
# same color scale as the main contour panel; `nothing` omits it.
function generate_css_contour_figure(; orientation::Symbol = :vertical, fontsize::Real = 30,
        zoom_extent = (10e3, 10e3), T_mult_css = [0.001, 0.2, 0.4, 0.6, 0.8, 0.99])
    T_css       = 5.0 * tau_R  # recurrence interval tied to tau_R, not an independent literature/historical value
    v_plate_css = b_nom / T_css  # b_nom (Natron best-fit) held fixed; v_plate derived

    L_aT       = a_nom * T_css
    vel_nondim_css  = T_css / b_nom     # = 1/V
    disp_nondim_css = 1.0 / b_nom
    # Own phase list as fractions of T (evenly spread over one cycle), instead
    # of the shared tau_R multiples in halfspace_common.jl.
    taus       = T_mult_css .* T_css
    labels_css = ["$(m)·T" for m in T_mult_css]

    # get_v_css/get_u_css's closed (Digamma) form is valid only for tau in
    # [0, T): periodicity was proven by construction of the infinite sum it
    # closes, not by the Digamma function itself (which is not periodic), so any
    # tau_i outside one cycle must be reduced mod T_css first. This matters
    # because snapshot phases are multiples of tau_R (up to 5*tau_R), which can
    # reach or exceed T_css = 5*tau_R.
    make_V_css(tau_i) = xp -> get_v_css(xp, mod(tau_i, T_css), b_nom, v_plate_css, a_nom)
    make_U_css(tau_i) = xp -> get_u_css(xp, mod(tau_i, T_css), b_nom, v_plate_css, a_nom)

    # Fixed range, narrowed to the two decades where the bulk of the data
    # actually sits (see the per-phase ranges this script prints).
    color_min_css = 0.01
    color_max_css = 1.0

    T_css_yr_str = round(T_css / yr, digits = 1)
    println("\nRSS model (contour) (T = $T_css_yr_str yr = 5·τ_R):  v_plate (nondim) = $(v_plate_css * vel_nondim_css)")

    # Width/height scale with the panel count (profile + interior [+ zoom]):
    # the base 2000/430 figures were tuned for 2 panel columns/rows.
    n_panels = zoom_extent === nothing ? 2 : 3
    panel_scale = n_panels / 2
    fig_size = orientation == :vertical ?
        (2000 * panel_scale, 430 * n_rows) : (430 * n_rows, 2000 * panel_scale)
    fig_css = Figure(size = fig_size, fontsize = fontsize)
    plot_model_rows!(fig_css, taus, labels_css, make_V_css, make_U_css;
        orientation = orientation, fontsize = fontsize, interior_panel = :contour,
        zoom_extent = zoom_extent,
        title = "RSS model (T = $T_css_yr_str yr = 5·τ_R)",
        colorbar_label = "|v| / V  (dimensionless)",
        vel_nondim = vel_nondim_css, vel_symbol = "V",
        disp_nondim = disp_nondim_css, disp_symbol = "b", u_lim = 2.0, vel_lim = 2.0,
        length_scale = L_aT, length_symbol = "aT", guide_formula = nothing,
        color_min = color_min_css, color_max = color_max_css, ref_speed = 0.5, colorscale = log10)
    out_path = "fig/fig_halfspace_css_contour$(orientation_suffix(orientation)).png"
    save(out_path, fig_css, px_per_unit = 2)
    println("Saved: $out_path")

    return fig_css
end

generate_css_contour_figure()
