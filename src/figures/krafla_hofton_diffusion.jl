function plot_krafla_hofton_diffusion()
    df       = CSV.read("data/hofton.csv", DataFrame)
    obs_x_m  = df.X .* 1000.0
    obs_u    = df.Y  ./ 100.0
    obs_sig  = df.SD ./ 100.0
    obs_x_km = obs_x_m ./ 1e3
    N        = length(obs_x_m)

    yr    = 3600.0 * 24 * 365
    tau_i = (1987 - 1975) * yr
    tau_f = (1992 - 1975) * yr

    kappa_hofton = 3.11
    b_hofton     = 1.27

    chi2_c = c -> begin
        u_pred = diff_foulger(obs_x_m, tau_i, tau_f, b_hofton, kappa_hofton) .+ c[1]
        sum(((obs_u .- u_pred) ./ obs_sig).^2)
    end
    result = optimize(chi2_c, [-0.2], [0.2], [0.0], Fminbox(NelderMead()))
    c_opt  = Optim.minimizer(result)[1]
    chi2r  = Optim.minimum(result) / (N - 1)

    @printf("Foulger diffusion model at Hofton & Foulger (1996) published parameters:\n")
    @printf("  kappa = %.3e m^2/s  (published)\n", kappa_hofton)
    @printf("  b     = %.3f m  (published)\n", b_hofton)
    @printf("  c     = %.3f m  (reference-frame offset)\n", c_opt)
    @printf("  chi2/nu = %.3f  (N=%d)\n", chi2r, N)

    x_km = range(-180, stop=180, length=2000)
    x_m  = x_km .* 1e3

    fig = Figure(size=(1600, 900), fontsize=28)
    ax  = Axis(fig[1, 1], xlabel = "Distance from Krafla volcanic system (km)",
        ylabel = "Displacement (cm)", limits = (-180, 180, -12, 20))

    lines!(ax, x_km, (diff_foulger(x_m, tau_i, tau_f, b_hofton, kappa_hofton) .+ c_opt) .* 100,
        color=:red, linewidth=2.5,
        label="Foulger diffusion | κ=$(sci(kappa_hofton)) m²/s  b=$(round(b_hofton, digits=2)) m")

    scatter!(ax, obs_x_km, obs_u .* 100, color=:black, markersize=10, label="Hofton & Foulger (1996), Fig. 16")
    errorbars!(ax, obs_x_km, obs_u .* 100, obs_sig .* 100, color=:black)

    axislegend(ax, position=:rb, framevisible=true)
    save("fig/fig_hofton_krafla_diffusion.png", fig, px_per_unit=2)
    println("Saved: fig/fig_hofton_krafla_diffusion.png")
    return fig
end
