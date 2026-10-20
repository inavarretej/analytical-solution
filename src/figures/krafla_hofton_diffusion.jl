using CairoMakie, CSV, DataFrames, Printf, Optim

include("../models/viscous.jl")
sci(x) = @sprintf("%.2e", x)

# Joint refit of the Foulger et al. (1992) diffusion model, under the same
# conditions as the RSS fits (Discussion, non-local versus Fickian
# relaxation): kappa = 2.70 m^2/s, b = 2.81 m, chi2_r = 4.90. Writes
# fig/fig_krafla_diffusion_joint.png (not used in the manuscript).
# Figure 1 of the paper: Foulger diffusion model against both Krafla
# rifting-episode datasets in one panel. A single (kappa, b) is fit by
# minimizing the summed chi-squared across both datasets, each evaluated at
# its own (tau_i, tau_f) -- same joint-fit convention as fit_krafla_css_joint
# in the Krafla fits.
#
# Foulger et al. (1992) ties its network to a set of near-rift-axis reference
# points the original authors judged undeformed, so no offset is applied to
# that dataset. Hofton & Foulger (1996) state explicitly, in the caption of
# their Figure 5 (the same profile recreated here), that "the zero reference
# point for displacements is arbitrary" -- the two datasets' zero levels are
# not tied to a common, recoverable frame. The offset between them is
# therefore a genuine free (nuisance) parameter c, fit jointly with (kappa, b)
# rather than precomputed from the data's mean.
function plot_krafla_hofton_diffusion()
    yr = 3600.0 * 24 * 365

    foulger          = CSV.read("data/extracted_data.csv", DataFrame)
    obs_x_m_foulger  = foulger.X .* 1000.0
    obs_u_foulger    = foulger.Y  ./ 100.0
    obs_sig_foulger  = foulger.SD ./ 100.0
    obs_x_km_foulger = obs_x_m_foulger ./ 1e3
    N_foulger        = length(obs_x_m_foulger)
    tau_i_foulger    = (1987 - 1975) * yr   # common time origin t0 = 1975, as in
    tau_f_foulger    = (1990 - 1975) * yr   # krafla_css_v_constrained.jl

    hofton           = CSV.read("data/hofton.csv", DataFrame)
    obs_x_m_hofton   = hofton.X .* 1000.0
    obs_u_hofton     = hofton.Y  ./ 100.0
    obs_sig_hofton   = hofton.SD ./ 100.0
    obs_x_km_hofton  = obs_x_m_hofton ./ 1e3
    N_hofton         = length(obs_x_m_hofton)
    tau_i_hofton     = (1987 - 1975) * yr
    tau_f_hofton     = (1992 - 1975) * yr

    chi2_joint = p -> begin
        kappa, b, c = p
        (kappa <= 0 || b <= 0) && return 1e10
        chi2_f = sum(((obs_u_foulger .- diff_foulger(obs_x_m_foulger, tau_i_foulger, tau_f_foulger, b, kappa)) ./ obs_sig_foulger) .^ 2)
        chi2_h = sum((((obs_u_hofton .- c) .- diff_foulger(obs_x_m_hofton, tau_i_hofton, tau_f_hofton, b, kappa)) ./ obs_sig_hofton) .^ 2)
        chi2_f + chi2_h
    end

    lb     = [0.05, 0.2, -0.20]
    ub     = [10.0, 10.0, 0.20]
    p0     = [1.86, 2.08, 0.0545]   # kappa, b: midpoint of the datasets' previous independent estimates; c: prior mean-offset estimate
    result = optimize(chi2_joint, lb, ub, p0, Fminbox(NelderMead()))
    kappa_opt, b_opt, c_opt = Optim.minimizer(result)

    obs_u_hofton_centered = obs_u_hofton .- c_opt

    chi2_foulger_val = sum(((obs_u_foulger .- diff_foulger(obs_x_m_foulger, tau_i_foulger, tau_f_foulger, b_opt, kappa_opt)) ./ obs_sig_foulger) .^ 2)
    chi2_hofton_val  = sum(((obs_u_hofton_centered .- diff_foulger(obs_x_m_hofton, tau_i_hofton, tau_f_hofton, b_opt, kappa_opt)) ./ obs_sig_hofton) .^ 2)
    N_total     = N_foulger + N_hofton
    nu          = N_total - 3   # kappa, b, c
    chi2r_total = (chi2_foulger_val + chi2_hofton_val) / nu

    @printf("Foulger diffusion model, joint fit to Foulger (1992) + Hofton & Foulger (1996):\n")
    @printf("  kappa = %.3e m^2/s  (joint fit)\n", kappa_opt)
    @printf("  b     = %.3f m  (joint fit)\n", b_opt)
    @printf("  c     = %.3f m  (Hofton & Foulger reference-frame offset, jointly fit)\n", c_opt)
    @printf("  chi2 (Foulger 1992)          = %.3f  (N=%d)\n", chi2_foulger_val, N_foulger)
    @printf("  chi2 (Hofton & Foulger 1996) = %.3f  (N=%d)\n", chi2_hofton_val, N_hofton)
    @printf("  chi2/nu (joint)              = %.3f  (N=%d, nu=%d)\n", chi2r_total, N_total, nu)

    x_km = range(-180, stop=180, length=2000)
    x_m  = x_km .* 1e3

    fig = Figure(size=(1600, 1050), fontsize=60)
    ax  = Axis(fig[1, 1], xlabel = "Distance from Krafla volcanic system (km)",
        ylabel = "Displacement (cm)", limits = (-180, 180, -15, 20))

    # Campaign-consistent colors (validated categorical pair, CVD-safe): each
    # dataset's data points and its own curve share one color; marker shape
    # (square/triangle) gives a second, color-independent way to tell the
    # two campaigns apart.
    color_foulger = RGBf(0x2a/255, 0x78/255, 0xd6/255)   # blue
    color_hofton  = RGBf(0xeb/255, 0x68/255, 0x34/255)   # orange

    lines!(ax, x_km, diff_foulger(x_m, tau_i_foulger, tau_f_foulger, b_opt, kappa_opt) .* 100,
        color=color_foulger, linewidth=2.5,
        label="Diffusion model, τ = 12–15 yr")
    lines!(ax, x_km, diff_foulger(x_m, tau_i_hofton, tau_f_hofton, b_opt, kappa_opt) .* 100,
        color=color_hofton, linewidth=2.5, linestyle=:dash,
        label="Diffusion model, τ = 12–17 yr")

    scatter!(ax, obs_x_km_foulger, obs_u_foulger .* 100, color=color_foulger, marker=:rect, markersize=12,
        label="Foulger et al. (1992)")
    errorbars!(ax, obs_x_km_foulger, obs_u_foulger .* 100, obs_sig_foulger .* 100, color=color_foulger)

    scatter!(ax, obs_x_km_hofton, obs_u_hofton_centered .* 100, color=color_hofton, marker=:utriangle, markersize=13,
        label="Hofton & Foulger (1996), Fig. 5")
    errorbars!(ax, obs_x_km_hofton, obs_u_hofton_centered .* 100, obs_sig_hofton .* 100, color=color_hofton)

    Legend(fig[2, 1], ax, nbanks=2, framevisible=false, labelsize=40, tellwidth=false, tellheight=true)
    save("fig/fig_krafla_diffusion_joint.png", fig, px_per_unit=2)
    println("Saved: fig/fig_krafla_diffusion_joint.png")
    return fig
end

mkpath("fig")
plot_krafla_hofton_diffusion()
