using CairoMakie, LaTeXStrings

mkpath("fig")

# ── Best-fit characteristic velocities (from krafla_fit.jl / asal_fit.jl) ────
a_krafla = 5.17e-5   # m/s (viscous CSS, T=200yr fixed)
a_asal   = 9.47e-6   # m/s (viscous CSS, T=500yr fixed)
yr       = 3600.0 * 24 * 365.25

# ── Adopted recurrence intervals and network apertures (Section 4.2 text) ───
T_krafla_yr    = 200.0
T_asal_yr      = 500.0
network_krafla_km = 230.0   # Foulger et al. (1992) station span
network_asal_km   = 17.0    # Vigny et al. (2007) GPS profile aperture

aT_km(a, T_yr) = a * T_yr * yr / 1e3

aT_krafla = aT_km(a_krafla, T_krafla_yr)
aT_asal   = aT_km(a_asal,   T_asal_yr)

println("Krafla:  T = $(T_krafla_yr) yr  →  aT = $(round(aT_krafla, digits=1)) km  (network = $(network_krafla_km) km)")
println("Asal:    T = $(T_asal_yr) yr  →  aT = $(round(aT_asal, digits=1)) km  (network = $(network_asal_km) km)")

# ── Figure ────────────────────────────────────────────────────────────────────
T_yr_grid = 10 .^ range(log10(10.0), log10(2000.0), length = 400)

fig = Figure(size = (1000, 700), fontsize = 20)
ax  = Axis(fig[1, 1],
    xlabel = "Recurrence interval  T  (yr)",
    ylabel = L"CSS\ deformation\ half-width\ \ aT\ \ (\mathrm{km})",
    xscale = log10, yscale = log10,
    xticks = ([10, 30, 100, 300, 1000], ["10", "30", "100", "300", "1000"]),
    yticks = ([1, 3, 10, 30, 100, 300, 1000], ["1", "3", "10", "30", "100", "300", "1000"]),
    limits = (10, 2000, 1, 2000))

lines!(ax, T_yr_grid, aT_km.(a_krafla, T_yr_grid),
    color = :firebrick, linewidth = 2.5,
    label = L"Krafla,\ a = 5.17\times10^{-5}\ \mathrm{m\,s^{-1}}")
lines!(ax, T_yr_grid, aT_km.(a_asal, T_yr_grid),
    color = :royalblue, linewidth = 2.5,
    label = L"Asal,\ a = 9.47\times10^{-6}\ \mathrm{m\,s^{-1}}")

# Adopted recurrence intervals
vlines!(ax, [T_krafla_yr], color = (:firebrick, 0.5), linestyle = :dash, linewidth = 1.5)
vlines!(ax, [T_asal_yr],   color = (:royalblue, 0.5), linestyle = :dash, linewidth = 1.5)

# Network apertures
hlines!(ax, [network_krafla_km], color = (:firebrick, 0.6), linestyle = :dashdot, linewidth = 1.5)
hlines!(ax, [network_asal_km],   color = (:royalblue, 0.6), linestyle = :dashdot, linewidth = 1.5)

# Implied aT points
scatter!(ax, [T_krafla_yr], [aT_krafla], color = :firebrick, markersize = 14, marker = :star5)
scatter!(ax, [T_asal_yr],   [aT_asal],   color = :royalblue, markersize = 14, marker = :star5)

text!(ax, T_krafla_yr * 1.15, aT_krafla * 1.15;
    text = "aT ≈ $(round(Int, aT_krafla)) km", color = :firebrick, fontsize = 15)
text!(ax, T_asal_yr * 1.15, aT_asal * 1.15;
    text = "aT ≈ $(round(Int, aT_asal)) km", color = :royalblue, fontsize = 15)

text!(ax, 11, network_krafla_km * 1.12;
    text = "Krafla network span (230 km)", color = (:firebrick, 0.8), fontsize = 13)
text!(ax, 11, network_asal_km * 1.12;
    text = "Asal network aperture (17 km)", color = (:royalblue, 0.8), fontsize = 13)

axislegend(ax, position = :lt, framevisible = true, labelsize = 15)

save("fig/fig_css_width.png", fig)
println("Saved  fig/fig_css_width.png")
