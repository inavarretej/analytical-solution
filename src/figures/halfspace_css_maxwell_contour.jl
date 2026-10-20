include("halfspace_common.jl")

# EVE (Maxwell) RSS, the counterpart of halfspace_css_contour.jl (Fig. 5):
# same layout, phases (fractions of T), T = 5 tau_R, dense grid, zoom, color
# range and plotting routine (plot_model_rows!). Only the model changes: the
# closed periodic EVE RSS forms (get_v/u_css_maxwell_closed, infinite periodic
# history, converging to get_v/u_css as tau_R -> 0), and the maxevals-capped
# interior kernels (the uncapped halfspace_vx/vz hang on nested Maxwell
# quadrature).
#
# Expensive, like the EVE STEP figure. Run with threads, e.g.
#   julia -t auto --project=. src/figures/halfspace_css_maxwell_contour.jl

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

# First row: tau = 0+ (label "0⁺"), as in the EVE STEP figure. The numerical
# closed form at tau = 0.001 T costs ~4 s per interior-kernel evaluation
# (~40 h for that row). get_v_css_maxwell_0plus (maxwell.jl) gives the exact
# 0+ limit but is itself slow (~0.3-4 s per call), so it is tabulated once on
# a log-spaced |x| grid (1 m .. 1e10 m, 3001 points; v is odd in x) and
# linearly interpolated in log|x|. Outside the table: clamped to the end
# values (at 1e8 m v is already V/2 within 2e-4). The interpolation error is
# checked against direct calls at off-grid points every run (logged; the run
# stops if it exceeds 1e-4). Documentation: Obsidian note
# 30-models/analytical/EVE contour figures - tau = 0+ limit.md
function make_v_css_0plus_interp(v_plate_css; n = 3001, logx_min = 0.0, logx_max = 10.0)
    logxs = collect(range(logx_min, logx_max, length = n))
    vals = cached_grid("css_maxwell_contour", "v0plus_table") do
        tab = Vector{Float64}(undef, n)
        n_done = Threads.Atomic{Int}(0); t0 = time()
        Threads.@threads for i in 1:n
            tab[i] = get_v_css_maxwell_0plus(10.0^logxs[i], b_nom, v_plate_css, a_nom, tau_R)
            c = Threads.atomic_add!(n_done, 1) + 1
            c % 150 == 0 && log_progress("v(0⁺) table  $(round(100c/n, digits = 1))%  elapsed $(fmt_hms(time() - t0))")
        end
        tab
    end
    Δ = logxs[2] - logxs[1]
    interp(x) = begin
        x == 0 && return 0.0
        s = clamp(log10(abs(x)), logx_min, logx_max)
        i = min(floor(Int, (s - logx_min) / Δ) + 1, n - 1)
        w = (s - logxs[i]) / Δ
        sign(x) * ((1 - w) * vals[i] + w * vals[i + 1])
    end
    # Interpolation check at off-grid points (midpoints between nodes).
    check_logx = logxs[1:150:end-1] .+ Δ / 2
    errs = [abs(interp(10.0^s) / get_v_css_maxwell_0plus(10.0^s, b_nom, v_plate_css, a_nom, tau_R) - 1) for s in check_logx]
    log_progress("v(0⁺) table interpolation check: max relative error $(maximum(errs)) over $(length(errs)) off-grid points")
    maximum(errs) < 1e-4 || error("v(0⁺) interpolation error too large: $(maximum(errs))")
    return interp
end

function generate_css_maxwell_contour_figure(; orientation::Symbol = :vertical, fontsize::Real = 30,
        zoom_extent = (10e3, 10e3), T_mult_css = [0.0, 0.2, 0.4, 0.6, 0.8, 0.99])
    T_css       = 5.0 * tau_R  # same illustrative T as the EV RSS figure
    v_plate_css = b_nom / T_css

    L_aT            = a_nom * T_css
    vel_nondim_css  = T_css / b_nom     # = 1/V
    disp_nondim_css = 1.0 / b_nom
    taus       = T_mult_css .* T_css
    labels_css = [m == 0 ? "0⁺" : "$(m)·T" for m in T_mult_css]

    v_css_0plus = make_v_css_0plus_interp(v_plate_css)

    # Velocity feeds the interior kernels thousands of times: default
    # k_limit_scale. The surface displacement profile (400 points per phase)
    # uses k_limit_scale = 30, which removes the small spurious ripple of the
    # closed u form at scale 1 (see get_u_css_maxwell_closed in maxwell.jl).
    # At 0+ the RSS displacement since the last intrusion is zero for x != 0.
    make_V_css(tau_i) = tau_i == 0 ? v_css_0plus :
        (xp -> get_v_css_maxwell_closed(xp, mod(tau_i, T_css), b_nom, v_plate_css, a_nom, tau_R))
    make_U_css(tau_i) = tau_i == 0 ? (xp -> 0.0) :
        (xp -> get_u_css_maxwell_closed(xp, mod(tau_i, T_css), b_nom, v_plate_css, a_nom, tau_R;
            k_limit_scale = 30.0))

    color_min_css = 0.01
    color_max_css = 1.0

    T_css_yr_str = round(T_css / yr, digits = 1)
    println("\nEVE RSS model (contour) (T = $T_css_yr_str yr = 5·τ_R), $(Threads.nthreads()) threads")

    n_panels = zoom_extent === nothing ? 2 : 3
    panel_scale = n_panels / 2
    fig_size = orientation == :vertical ?
        (2000 * panel_scale, 430 * n_rows) : (430 * n_rows, 2000 * panel_scale)
    fig_css = Figure(size = fig_size, fontsize = fontsize)
    plot_model_rows!(fig_css, taus, labels_css, make_V_css, make_U_css;
        orientation = orientation, fontsize = fontsize, interior_panel = :contour,
        zoom_extent = zoom_extent,
        title = "EVE RSS model (T = $T_css_yr_str yr = 5·τ_R)",
        colorbar_label = "|v| / V  (dimensionless)",
        vel_nondim = vel_nondim_css, vel_symbol = "V",
        disp_nondim = disp_nondim_css, disp_symbol = "b", u_lim = 2.0, vel_lim = 2.0,
        length_scale = L_aT, length_symbol = "aT", guide_formula = nothing,
        color_min = color_min_css, color_max = color_max_css, ref_speed = 0.5, colorscale = log10,
        hvx = halfspace_vx_capped, hvz = halfspace_vz_capped,  # dashed guide at a*tau, as in Fig. 5
        cache_prefix = "css_maxwell_contour")  # resume from cache/ after an interruption
    out_path = "fig/fig_halfspace_css_maxwell_contour$(orientation_suffix(orientation)).png"
    save(out_path, fig_css, px_per_unit = 2)
    println("Saved: $out_path")

    return fig_css
end

@time generate_css_maxwell_contour_figure()
