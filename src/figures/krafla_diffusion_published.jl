using CairoMakie, CSV, DataFrames, Printf

mkpath("fig")

include("../models/viscous.jl")

# Fig. 1 of the paper: Krafla post-rifting profiles and the Fickian diffusion
# curves exactly as each original paper reports them. No refitting and no
# reference-frame correction (that is handled in our own fits, Section 4).

# ── Data ─────────────────────────────────────────────────────────────────────
yr = 3600.0 * 24 * 365

foulger = CSV.read("data/extracted_data.csv", DataFrame)   # cm
hofton  = CSV.read("data/hofton.csv", DataFrame)           # cm, their Fig. 5 (same points as their Fig. 16)

# ── Published models ─────────────────────────────────────────────────────────
# Both papers use the method of Foulger et al. (1992, eq. 4): the displacement
# over a survey interval dt is treated as the instantaneous velocity of
# u = U0 erfc(x / 2 sqrt(kappa t)) at time t since diking, times dt.
# Time origin 1979 (1990 = 11 yr, as stated by Foulger et al., 1992).
#
# Foulger et al. (1992): 1987-1990, t = 11 yr, dt = 3 yr,
#   U0 = 1.0 +/- 0.1 m, kappa = 1.1 +/- 0.3 m^2/s. Network tied to near-rift
#   reference points judged undeformed, so no offset.
U0_foulger = 1.0
kappa_foulger = 1.1
t_foulger = 11.0 * yr
dt_foulger = 3.0 * yr

# Hofton & Foulger (1996), Fig. 16 (solid curve): 1987-1992, t = 13 yr,
#   dt = 5 yr, kappa = 3.11 +/- 0.68 m^2/s, half dike thickness 1.27 +/- 0.14 m.
#   With this convention the curve peaks at x ~ 50 km with half-amplitude
#   ~11.8 cm, as in their figure. Their zero reference is arbitrary and their
#   plotted curve is shifted by a constant; c_hofton is read from their
#   Fig. 16 (midpoint between the curve's maximum and minimum), not fitted.
U0_hofton = 1.27
kappa_hofton = 3.11
t_hofton = 13.0 * yr
dt_hofton = 5.0 * yr
c_hofton = 6.7   # cm

published(x_m, U0, kappa, t, dt) = get_v_foulger.(x_m, t, kappa, 2 * U0) .* dt .* 100   # cm

# ── Figure ───────────────────────────────────────────────────────────────────
x_km = range(-180, stop=180, length=2000)
x_m  = x_km .* 1e3

u_foulger = published(x_m, U0_foulger, kappa_foulger, t_foulger, dt_foulger)
u_hofton  = published(x_m, U0_hofton, kappa_hofton, t_hofton, dt_hofton) .+ c_hofton

i_max = argmax(u_hofton); i_min = argmin(u_hofton)
@printf("Hofton & Foulger (1996) curve: max %.1f cm at %.0f km, min %.1f cm at %.0f km\n",
    u_hofton[i_max], x_km[i_max], u_hofton[i_min], x_km[i_min])
@printf("Foulger et al. (1992) curve:   max %.1f cm at %.0f km\n",
    maximum(u_foulger), x_km[argmax(u_foulger)])

# Full-width (figure*) layout, same size and font scale as one row of
# fig_krafla_v_constrained_combined.png so the two Krafla figures match.
# Published half-widths U0 are shown as the full opening b = 2 U0.
fig = Figure(size = (2600, 1000), fontsize = 40)

function published_panel!(pos, title, d, u_model, U0, kappa; ylabel = "Displacement (cm)")
    ax = Axis(pos, xlabel = "Distance from Krafla volcanic system (km)", ylabel = ylabel,
              title = title, limits = (-180, 180, -12, 22))
    lines!(ax, x_km, u_model, color = :red, linewidth = 3,
           label = "Published diffusion model: b = $(round(2U0, digits = 2)) m, κ = $(kappa) m²/s")
    errorbars!(ax, d.X, d.Y, d.SD, color = :black)
    scatter!(ax, d.X, d.Y, color = :black, markersize = 12, label = "GPS data")
    axislegend(ax, position = :rb, framevisible = true, labelsize = 30)
    ax
end

ax_a = published_panel!(fig[1, 1], "Foulger et al. (1992), 1987–1990", foulger,
                        u_foulger, U0_foulger, kappa_foulger)
ax_b = published_panel!(fig[1, 2], "Hofton & Foulger (1996), 1987–1992", hofton,
                        u_hofton, U0_hofton, kappa_hofton; ylabel = "")
linkyaxes!(ax_a, ax_b)

for (ax, lab) in zip((ax_a, ax_b), ("(a)", "(b)"))
    text!(ax, 0.02, 0.98, text = lab, space = :relative, align = (:left, :top), fontsize = 44)
end

save("fig/fig_krafla_diffusion_published.png", fig, px_per_unit = 2)
println("fig/fig_krafla_diffusion_published.png saved")
