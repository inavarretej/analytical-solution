function plot_asal_chi2_landscape_css(; obs_year, obs_dL, obs_sig, x0_m, t_ref_s, t0_year, yr,
        a_css_fixed, b_css_fixed, T_fixed)
    mkpath("fig")
    obs_t = (obs_year .- t0_year) .* yr
    N     = length(obs_t)

    weighted_chi2(pred, obs, sig) = sum(((obs .- pred) ./ sig) .^ 2)
    predict_css(t, a, b)      = 2 .* get_u_css.(x0_m, t, b, b / T_fixed, a) .* 1000
    predict_css_tref(t, a, b) = predict_css(t, a, b) .- predict_css([t_ref_s], a, b)[1]

    function argmin_grid(a_vec, b_vec, chi2N)
        imin = argmin(chi2N)
        return a_vec[imin[1]], b_vec[imin[2]]
    end

    # NelderMead can "converge" on Asal's shallow, ill-conditioned (a,b) ridge
    # well short of the true minimum (checked: off by Δχ² ~ 12, i.e. outside the
    # 1σ contour it should sit at the centre of). Rather than trust the optimizer
    # alone, re-centre on the grid's OWN argmin at each step -- self-consistent
    # by construction, and the final high-resolution grid is precise enough that
    # its argmin is an accurate stand-in for the true continuous optimum.
    function zoomed_grid(chi2_fn, a_opt, b_opt; n_search = 45, n_final = 300, max_iter = 8)
        a_half = 0.2
        b_half = max(0.15 * b_opt, 0.01)
        local a_vec, b_vec, chi2N, chi2N_min, lv
        for _ in 1:max_iter
            a_vec = 10 .^ range(log10(a_opt) - a_half, log10(a_opt) + a_half, length = n_search)
            b_vec = range(max(1e-4, b_opt - b_half), b_opt + b_half, length = n_search)
            chi2N = [chi2_fn(a, b) / N for a in a_vec, b in b_vec]
            a_opt, b_opt = argmin_grid(a_vec, b_vec, chi2N)
            chi2N_min = minimum(chi2N)
            lv = chi2N_min .+ [2.30, 6.17, 11.83] ./ N
            mask3 = chi2N .≤ lv[3]
            touches_edge = any(mask3[1, :]) || any(mask3[end, :]) || any(mask3[:, 1]) || any(mask3[:, end])
            touches_edge || break
            a_half *= 1.8
            b_half *= 1.8
        end
        a_vec = 10 .^ range(log10(a_opt) - a_half, log10(a_opt) + a_half, length = n_final)
        b_vec = range(max(1e-4, b_opt - b_half), b_opt + b_half, length = n_final)
        chi2N = [chi2_fn(a, b) / N for a in a_vec, b in b_vec]
        a_opt, b_opt = argmin_grid(a_vec, b_vec, chi2N)
        chi2N_min = minimum(chi2N)
        lv = chi2N_min .+ [2.30, 6.17, 11.83] ./ N
        return a_vec, b_vec, chi2N, chi2N_min, lv, a_opt, b_opt
    end

    println("Building viscous CSS landscape …")
    a_v, b_v, chi2N_v, min_v, lv_v, a_opt, b_opt = zoomed_grid(
        (a, b) -> weighted_chi2(predict_css_tref(obs_t, a, b), obs_dL, obs_sig), a_css_fixed, b_css_fixed)

    println("N = $N data points")
    println("χ²/N minimum = $(round(min_v, sigdigits=4)) at a=$a_opt, b=$b_opt  (grid argmin, refined from optimizer)")

    conf_colors = [:white, :gold, :tomato]
    conf_labels = ["68.3%  (1σ)", "95.4%  (2σ)", "99.7%  (3σ)"]

    fig = Figure(size = (1100, 740), fontsize = 20)
    ax  = Axis(fig[1, 1],
        xscale = log10,
        xlabel = L"a\ \ (\mathrm{m\,s^{-1}})",
        ylabel = L"b\ \ (\mathrm{m})",
        title  = L"\chi^2/N\ \mathrm{landscape\ —\ viscous\ CSS\ (Asal,\ }T=500\ \mathrm{yr)}")

    n_levels  = 20
    clip_hi   = lv_v[3] + (lv_v[3] - min_v) * 0.4
    cmap_disc = cgrad(ColorSchemes.bamako100, n_levels, categorical = true)

    hm = heatmap!(ax, collect(a_v), collect(b_v), chi2N_v, colormap = cmap_disc, colorrange = (min_v, clip_hi))

    cb_ticks_val = round.([min_v; lv_v], sigdigits = 3)
    cb_ticks_lbl = [L"\chi^2_{\min}/N", "1σ", "2σ", "3σ"]
    Colorbar(fig[1, 2], hm, label = L"\chi^2/N", ticks = (cb_ticks_val, cb_ticks_lbl))

    for (lev, col) in zip(lv_v, conf_colors)
        contour!(ax, collect(a_v), collect(b_v), chi2N_v, levels = [lev], color = col, linewidth = 2.5)
    end

    scatter!(ax, [a_opt], [b_opt], color = :cyan, markersize = 16, marker = :star5, strokewidth = 0)

    Legend(fig[1, 3],
        [LineElement(color = conf_colors[1], linewidth = 2.5),
         LineElement(color = conf_colors[2], linewidth = 2.5),
         LineElement(color = conf_colors[3], linewidth = 2.5),
         MarkerElement(color = :cyan, marker = :star5, markersize = 16)],
        [conf_labels..., L"(a^*,\;b^*)"],
        framevisible = true, labelsize = 16)

    save("fig/fig_chi2_landscape_asal_css.png", fig)
    println("Saved  fig/fig_chi2_landscape_asal_css.png")
    return fig
end
