function plot_asal_timeseries_css(; obs_year, obs_dL, obs_sig, x0_m, t_ref_s, t0_year, yr,
        a_css_free, b_css_free, T_css_free, a_css_fixed, b_css_fixed, T_fixed,
        a_maxcss_free, b_maxcss_free, T_maxcss_free, tauR_maxcss_free,
        a_maxcss_fixed, b_maxcss_fixed, tauR_maxcss_fixed)
    predict_css(t, p)      = 2 .* get_u_css.(x0_m, t, p[2], p[2] / p[3], p[1]) .* 1000
    predict_css_tref(t, p) = predict_css(t, p) .- predict_css([t_ref_s], p)[1]

    mkpath("fig")
    t_yr_grid = range(0.05, stop = 26.0, length = 600)
    t_s_grid  = t_yr_grid .* yr

    fig = Figure(size = (1600, 900), fontsize = 28)
    ax  = Axis(fig[1, 1],
        xlabel = "Time since rifting (yr)",
        ylabel = "EP00–DF00 distance change (mm)",
        limits = (t0_year - 1, 2005, -80, 580))

    lines!(ax, t_yr_grid .+ t0_year, predict_css_tref(t_s_grid, [a_css_free, b_css_free, T_css_free]),
        color = :purple, linewidth = 2.5,
        label = "CSS (T free)  a=$(sci(a_css_free)) m/s  b=$(round(b_css_free; digits=3)) m  T=$(round(T_css_free/yr; digits=1)) yr")
    lines!(ax, t_yr_grid .+ t0_year, predict_css_tref(t_s_grid, [a_css_fixed, b_css_fixed, T_fixed]),
        color = :orange, linewidth = 2.5, linestyle = :dash,
        label = "CSS (T=500yr)  a=$(sci(a_css_fixed)) m/s  b=$(round(b_css_fixed; digits=3)) m")
    lines!(ax, t_yr_grid .+ t0_year,
        2 .* get_u_css_maxwell_closed.(x0_m, t_s_grid, b_maxcss_free, b_maxcss_free / T_maxcss_free, a_maxcss_free, tauR_maxcss_free; k_limit_scale=30.0) .* 1000 .-
        (2 .* get_u_css_maxwell_closed(x0_m, t_ref_s, b_maxcss_free, b_maxcss_free / T_maxcss_free, a_maxcss_free, tauR_maxcss_free; k_limit_scale=30.0) .* 1000),
        color = :navy, linewidth = 2.5,
        label = "MaxCSS (T free)  a=$(sci(a_maxcss_free)) m/s  b=$(round(b_maxcss_free; digits=3)) m  T=$(round(T_maxcss_free/yr; digits=1)) yr  τ_R=$(round(tauR_maxcss_free/yr; digits=2)) yr")
    lines!(ax, t_yr_grid .+ t0_year,
        2 .* get_u_css_maxwell_closed.(x0_m, t_s_grid, b_maxcss_fixed, b_maxcss_fixed / T_fixed, a_maxcss_fixed, tauR_maxcss_fixed; k_limit_scale=30.0) .* 1000 .-
        (2 .* get_u_css_maxwell_closed(x0_m, t_ref_s, b_maxcss_fixed, b_maxcss_fixed / T_fixed, a_maxcss_fixed, tauR_maxcss_fixed; k_limit_scale=30.0) .* 1000),
        color = :brown, linewidth = 2.5, linestyle = :dash,
        label = "MaxCSS (T=500yr)  a=$(sci(a_maxcss_fixed)) m/s  b=$(round(b_maxcss_fixed; digits=3)) m  τ_R=$(round(tauR_maxcss_fixed/yr; digits=2)) yr")

    scatter!(ax, obs_year, obs_dL,
        color = :black, marker = :circle, markersize = 9,
        label = "Vigny et al. (2007)")
    errorbars!(ax, obs_year, obs_dL, obs_sig, color = :black)

    axislegend(ax, position = :rb, framevisible = true, labelsize = 16)
    save("fig/fig_asal_timeseries_css_comparison.png", fig)
    println("fig/fig_asal_timeseries_css_comparison.png  saved")
    return fig
end
