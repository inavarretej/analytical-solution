function plot_krafla_fits_css(; obs_x_km, obs_u, obs_sig, x_km, x_m,
        a_css_free, b_css_free, T_css_free, a_css_fixed, b_css_fixed, T_fixed,
        a_maxcss_free, b_maxcss_free, T_maxcss_free, tauR_maxcss_free,
        a_maxcss_fixed, b_maxcss_fixed, tauR_maxcss_fixed, tau_8, tau, yr)
    pred_css(x, a, b, T) = diff_css(x, tau_8, tau, b, b / T, a)
    pred_css_maxwell(x, a, b, T, tau_R; k_limit_scale=1.0) = begin
        get_u_css_maxwell_closed.(x, tau,   b, b / T, a, tau_R; k_limit_scale=k_limit_scale) .-
        get_u_css_maxwell_closed.(x, tau_8, b, b / T, a, tau_R; k_limit_scale=k_limit_scale)
    end

    fig = Figure(size=(1600, 900), fontsize=28)
    ax  = Axis(fig[1, 1], xlabel = "Distance (km)", ylabel = "Displacement (cm)", limits = (-150, 150, -13, 13))

    lines!(ax, x_km, pred_css(x_m, a_css_free, b_css_free, T_css_free) .* 100,
        color=:purple, linewidth=2,
        label="CSS 3p (T free) | a=$(sci(a_css_free)) m/s  b=$(sci(b_css_free)) m  T=$(round(T_css_free/yr,digits=1)) yr")
    lines!(ax, x_km, pred_css(x_m, a_css_fixed, b_css_fixed, T_fixed) .* 100,
        color=:orange, linewidth=2, linestyle=:dash,
        label="CSS 2p (T=200yr fixed) | a=$(sci(a_css_fixed)) m/s  b=$(sci(b_css_fixed)) m")
    lines!(ax, x_km, pred_css_maxwell(x_m, a_maxcss_free, b_maxcss_free, T_maxcss_free, tauR_maxcss_free; k_limit_scale=30.0) .* 100,
        color=:navy, linewidth=2,
        label="MaxCSS 3p (T free) | a=$(sci(a_maxcss_free)) m/s  b=$(sci(b_maxcss_free)) m  T=$(round(T_maxcss_free/yr,digits=1)) yr  τ_R=$(round(tauR_maxcss_free/yr,digits=2)) yr")
    lines!(ax, x_km, pred_css_maxwell(x_m, a_maxcss_fixed, b_maxcss_fixed, T_fixed, tauR_maxcss_fixed; k_limit_scale=30.0) .* 100,
        color=:brown, linewidth=2, linestyle=:dash,
        label="MaxCSS 2p (T=200yr fixed) | a=$(sci(a_maxcss_fixed)) m/s  b=$(sci(b_maxcss_fixed)) m  τ_R=$(round(tauR_maxcss_fixed/yr,digits=2)) yr")

    scatter!(ax, obs_x_km, obs_u .* 100, color=:black, markersize=10, label="Foulger et al. (1992)")
    errorbars!(ax, obs_x_km, obs_u .* 100, obs_sig .* 100, color=:black)

    axislegend(ax, position=:rb, framevisible=true, labelsize=18)
    save("fig/fig08_fits_css_comparison.png", fig)
    println("\nfig/fig08_fits_css_comparison.png saved")
    return fig
end
