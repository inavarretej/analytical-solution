function plot_asal_velocity_profile_css(; a_css_free, b_css_free, T_css_free,
        a_css_fixed, b_css_fixed, T_fixed, a_maxcss_free, b_maxcss_free, T_maxcss_free,
        tauR_maxcss_free, a_maxcss_fixed, b_maxcss_fixed, tauR_maxcss_fixed, yr)
    mkpath("fig")
    t_epoch    = 19.0 * yr
    x_center_m = 1000.0

    df       = CSV.read("data/asal_profile.csv", DataFrame)
    obs_x_km = df.X
    obs_v    = df.Y
    obs_sig  = df.SD
    obs_x_m  = obs_x_km .* 1e3
    N        = length(obs_x_km)

    x_rel(x_m) = x_m .- x_center_m
    pred_css_free(x_m)  = get_v_css.(x_rel(x_m), t_epoch, b_css_free,  b_css_free  / T_css_free, a_css_free)  .* 1000 .* yr
    pred_css_fixed(x_m) = get_v_css.(x_rel(x_m), t_epoch, b_css_fixed, b_css_fixed / T_fixed,    a_css_fixed) .* 1000 .* yr
    pred_maxcss_free(x_m; k_limit_scale=1.0)  = get_v_css_maxwell_closed.(x_rel(x_m), t_epoch, b_maxcss_free,  b_maxcss_free  / T_maxcss_free, a_maxcss_free,  tauR_maxcss_free; k_limit_scale=k_limit_scale)  .* 1000 .* yr
    pred_maxcss_fixed(x_m; k_limit_scale=1.0) = get_v_css_maxwell_closed.(x_rel(x_m), t_epoch, b_maxcss_fixed, b_maxcss_fixed / T_fixed,       a_maxcss_fixed, tauR_maxcss_fixed; k_limit_scale=k_limit_scale) .* 1000 .* yr

    function fit_vplate(pred_fn)
        v_mod = pred_fn(obs_x_m)
        w     = 1.0 ./ obs_sig.^2
        sum(w .* (obs_v .- v_mod)) / sum(w)
    end

    Vp_css_free     = fit_vplate(pred_css_free)
    Vp_css_fixed    = fit_vplate(pred_css_fixed)
    Vp_maxcss_free  = fit_vplate(pred_maxcss_free)
    Vp_maxcss_fixed = fit_vplate(pred_maxcss_fixed)

    @printf("Background V_plate fits:\n")
    @printf("  CSS (T free)      V_plate = %.2f mm/yr\n", Vp_css_free)
    @printf("  CSS (T=500yr)     V_plate = %.2f mm/yr\n", Vp_css_fixed)
    @printf("  MaxCSS (T free)   V_plate = %.2f mm/yr\n", Vp_maxcss_free)
    @printf("  MaxCSS (T=500yr)  V_plate = %.2f mm/yr\n", Vp_maxcss_fixed)

    chi2(pred_fn, Vp) = sum(((obs_v .- (pred_fn(obs_x_m) .+ Vp)) ./ obs_sig).^2)
    @printf("\nChi² / (N-1)  [k_eff = 1 for V_plate]:\n")
    @printf("  CSS (T free)      χ²_r = %.2f\n", chi2(pred_css_free,     Vp_css_free)     / (N-1))
    @printf("  CSS (T=500yr)     χ²_r = %.2f\n", chi2(pred_css_fixed,    Vp_css_fixed)    / (N-1))
    @printf("  MaxCSS (T free)   χ²_r = %.2f\n", chi2(pred_maxcss_free,  Vp_maxcss_free)  / (N-1))
    @printf("  MaxCSS (T=500yr)  χ²_r = %.2f\n", chi2(pred_maxcss_fixed, Vp_maxcss_fixed) / (N-1))

    x_grid_km = collect(range(-12, 12, length = 800))
    x_grid_m  = x_grid_km .* 1e3

    fig = Figure(size = (1600, 900), fontsize = 28)
    ax  = Axis(fig[1, 1],
        xlabel = "Distance from rift axis (km)",
        ylabel = "Rift-perpendicular velocity, Somalia frame (mm yr⁻¹)",
        title  = "Asal velocity profile — CSS viscous vs Maxwell — t = 19 yr post-rifting",
        limits = (-12, 12, -2, 22))

    vspan!(ax, -3, 3, color = (:lightgray, 0.25))
    vspan!(ax, -1, 1, color = (:lightgray, 0.40))
    vlines!(ax, [0.0], color = (:gray, 0.4), linestyle = :dot,  linewidth = 1.0)
    vlines!(ax, [1.0], color = (:black, 0.3), linestyle = :dash, linewidth = 1.0)

    lines!(ax, x_grid_km, pred_css_free(x_grid_m) .+ Vp_css_free,
        color = :purple, linewidth = 2.5,
        label = "CSS (T free)  a=$(@sprintf("%.2e", a_css_free)) m/s  T=$(round(T_css_free/yr,digits=1)) yr  V_pl=$(round(Vp_css_free,digits=1)) mm/yr")
    lines!(ax, x_grid_km, pred_css_fixed(x_grid_m) .+ Vp_css_fixed,
        color = :orange, linewidth = 2.5, linestyle = :dash,
        label = "CSS (T=500yr)  a=$(@sprintf("%.2e", a_css_fixed)) m/s  V_pl=$(round(Vp_css_fixed,digits=1)) mm/yr")
    lines!(ax, x_grid_km, pred_maxcss_free(x_grid_m; k_limit_scale=30.0) .+ Vp_maxcss_free,
        color = :navy, linewidth = 2.5,
        label = "MaxCSS (T free)  a=$(@sprintf("%.2e", a_maxcss_free)) m/s  T=$(round(T_maxcss_free/yr,digits=1)) yr  V_pl=$(round(Vp_maxcss_free,digits=1)) mm/yr")
    lines!(ax, x_grid_km, pred_maxcss_fixed(x_grid_m; k_limit_scale=30.0) .+ Vp_maxcss_fixed,
        color = :brown, linewidth = 2.5, linestyle = :dash,
        label = "MaxCSS (T=500yr)  a=$(@sprintf("%.2e", a_maxcss_fixed)) m/s  V_pl=$(round(Vp_maxcss_fixed,digits=1)) mm/yr")

    errorbars!(ax, obs_x_km, obs_v, obs_sig, color = :black, whiskerwidth = 5)
    scatter!(ax, obs_x_km, obs_v, color = :black, markersize = 9, label = "Vigny et al. (2007)")

    axislegend(ax, position = :rt, framevisible = true, labelsize = 14)
    text!(ax, -11.5, 0.8, text = "shaded: 1978 activation zone (light) & Petit Rift (dark)",
        fontsize = 16, color = (:gray, 0.65))

    save("fig/fig_asal_velocity_profile_css.png", fig, px_per_unit = 2)
    println("Saved: fig/fig_asal_velocity_profile_css.png")
    return fig
end
