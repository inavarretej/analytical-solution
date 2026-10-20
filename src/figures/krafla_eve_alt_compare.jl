using CairoMakie, CSV, DataFrames, Printf, Optim, Dates, ColorSchemes, LaTeXStrings

# Krafla: EV RSS and paper EVE RSS (best fits read from
# reports/krafla_v_constrained_combined_params.csv, not refit) against the
# alternative EVE RSS of src/models/maxwell_alt.jl, fitted here with the same
# data, t0 = 1975, V = 1.8 cm/yr, T = b/V and analytic Hofton offset c as
# krafla_v_constrained_combined.jl. Comparison only: no existing script,
# figure or report is modified.
#   fig/fig_krafla_eve_alt_compare.png
#   reports/krafla_eve_alt_report.txt, reports/krafla_eve_alt_params.csv
mkpath("fig"); mkpath("reports")

include("../models/viscous.jl")
include("../models/maxwell.jl")
include("../models/maxwell_alt.jl")
include("../parameters.jl")   # H, E_L, nu_L
sci(x) = @sprintf("%.2e", x)

yr    = 3600.0 * 24 * 365     # same year length as krafla_v_constrained_combined.jl
V_lit = 0.018 / yr
T0    = 1975
tau_R_from_a(a) = H / (a * (1 - nu_L))

foulger = CSV.read("data/extracted_data.csv", DataFrame)
hofton  = CSV.read("data/hofton.csv", DataFrame)
xf, uf, sf = foulger.X .* 1e3, foulger.Y ./ 100, foulger.SD ./ 100
xh, uh, sh = hofton.X  .* 1e3, hofton.Y  ./ 100, hofton.SD  ./ 100
tf = ((1987 - T0) * yr, (1990 - T0) * yr)
th = ((1987 - T0) * yr, (1992 - T0) * yr)
N  = length(xf) + length(xh)

pred_ev(x, ti, tf_, a, b) = diff_css(x, ti, tf_, b, V_lit, a)
pred_eve(x, ti, tf_, a, b) = (tR = tau_R_from_a(a);
    get_u_css_maxwell_split.(x, tf_, b, V_lit, a, tR) .- get_u_css_maxwell_split.(x, ti, b, V_lit, a, tR))
pred_alt(x, ti, tf_, a, b) = (tR = tau_R_from_a(a);
    get_u_css_maxwell_alt.(x, tf_, b, V_lit, a, tR) .- get_u_css_maxwell_alt.(x, ti, b, V_lit, a, tR))

const wh = 1 ./ sh .^ 2
function joint_chi2(pred, a, b)
    (a <= 0 || b <= 0) && return (1e12, 0.0, 1e12, 1e12)
    pf = pred(xf, tf..., a, b)
    ph = pred(xh, th..., a, b)
    c  = sum(wh .* (uh .- ph)) / sum(wh)
    cf = sum(((uf .- pf) ./ sf) .^ 2)
    ch = sum(((uh .- ph .- c) ./ sh) .^ 2)
    (cf + ch, c, cf, ch)
end

P = CSV.read("reports/krafla_v_constrained_combined_params.csv", DataFrame)
row(m) = only(eachrow(filter(r -> r.model == m, P)))
ev, eve = row("EV"), row("EVE")

# ── Fit the alternative EVE: coarse grid in (log a, b), then Nelder-Mead ─────
println("Fitting alternative EVE RSS (T = b/V) jointly to Foulger + Hofton …")
const DCHI2 = [2.30, 6.17, 11.83]
function chi2_grid(a_vec, b_vec)
    G = Matrix{Float64}(undef, length(a_vec), length(b_vec))
    Threads.@threads for j in eachindex(b_vec)
        for i in eachindex(a_vec)
            G[i, j] = joint_chi2(pred_alt, a_vec[i], b_vec[j])[1]
        end
    end
    G
end
a_c = 10 .^ range(-5.5, -3.0, length = 26); b_c = collect(range(0.8, 3.0, length = 23))
Gc  = chi2_grid(a_c, b_c); ic = argmin(Gc)
res = optimize(q -> joint_chi2(pred_alt, 10^q[1], q[2])[1], [log10(a_c[ic[1]]), b_c[ic[2]]], NelderMead())
a_alt, b_alt = 10^Optim.minimizer(res)[1], Optim.minimizer(res)[2]
@printf("  a = %.4e m/s, b = %.4f m, χ² = %.3f\n", a_alt, b_alt, Optim.minimum(res))

# ── Landscape around the alternative best fit ────────────────────────────────
function zoomed_grid(a0, b0; n_search = 31, n_final = 121, max_iter = 8)
    a_half, b_half = 0.2, 0.1 * b0
    grid(n) = (10 .^ range(log10(a0) - a_half, log10(a0) + a_half, length = n),
               collect(range(max(0.05, b0 - b_half), b0 + b_half, length = n)))
    for _ in 1:max_iter
        a_vec, b_vec = grid(n_search)
        G = chi2_grid(a_vec, b_vec)
        m3 = G .≤ minimum(G) + DCHI2[3]
        (any(m3[1, :]) || any(m3[end, :]) || any(m3[:, 1]) || any(m3[:, end])) || break
        a_half *= 1.6; b_half *= 1.6
    end
    a_vec, b_vec = grid(n_final)
    (a_vec, b_vec, chi2_grid(a_vec, b_vec))
end
println("Building alternative-EVE landscape …")
t_land = @elapsed a_vec, b_vec, G = zoomed_grid(a_alt, b_alt)
cmin = min(minimum(G), Optim.minimum(res))
@printf("  done [%.0f s]\n", t_land)
in1 = G .≤ cmin + DCHI2[1]
a1 = extrema(a_vec[[i for i in axes(G, 1) if any(in1[i, :])]])
b1 = extrema(b_vec[[j for j in axes(G, 2) if any(in1[:, j])]])

models = [(name = "EV RSS",                a = ev.a,  b = ev.b,  pred = pred_ev,  color = :purple,     style = :solid),
          (name = "EVE RSS (paper)",       a = eve.a, b = eve.b, pred = pred_eve, color = :navy,       style = :dash),
          (name = "EVE RSS (alternative)", a = a_alt, b = b_alt, pred = pred_alt, color = :darkorange, style = :solid)]
stats = map(models) do m
    chi2, c, cf, ch = joint_chi2(m.pred, m.a, m.b)
    (; m.name, m.a, m.b, T = m.b / V_lit, tauR = m.name == "EV RSS" ? NaN : tau_R_from_a(m.a), chi2, c, cf, ch)
end

# ── Report ───────────────────────────────────────────────────────────────────
k = 3; aic(c) = c + 2k
open("reports/krafla_eve_alt_report.txt", "w") do io
    println(io, "Krafla — EV RSS vs paper EVE RSS vs alternative EVE RSS (T = b/V, V = 1.8 cm/yr, t0 = 1975, N = $N)")
    println(io, "Generated: $(Dates.now())")
    println(io, "EV and paper-EVE best fits read from reports/krafla_v_constrained_combined_params.csv (not refit); alternative EVE fitted here.")
    println(io, "Free: a, b; Hofton offset c solved analytically (k = 3).\n")
    @printf(io, "  %-22s %10s %8s %7s %7s %7s %9s %9s %9s %7s %8s\n",
            "Model", "a (m/s)", "b (m)", "T (yr)", "τR (yr)", "c (cm)", "χ²", "χ²_F92", "χ²_HF96", "χ²_r", "ΔAIC")
    for s in stats
        @printf(io, "  %-22s %10.3e %8.4f %7.1f %7.2f %7.2f %9.2f %9.2f %9.2f %7.3f %+8.2f\n",
                s.name, s.a, s.b, s.T / yr, s.tauR / yr, 100s.c, s.chi2, s.cf, s.ch, s.chi2 / (N - k),
                aic(s.chi2) - aic(stats[1].chi2))
    end
    @printf(io, "\n  Alternative EVE 1σ (Δχ² ≤ 2.30 projections): a ∈ [%.3e, %.3e] m/s, b ∈ [%.3f, %.3f] m, T ∈ [%.1f, %.1f] yr\n",
            a1..., b1..., (b1 ./ V_lit ./ yr)...)
    @printf(io, "  η (H = 10 km, E = 40 GPa, ν = 0.25): EV %.2e, EVE paper %.2e, EVE alternative %.2e Pa s\n",
            [E_L * H / (2 * s.a * (1 - nu_L^2)) for s in stats]...)
end
CSV.write("reports/krafla_eve_alt_params.csv",
          DataFrame(model = [s.name for s in stats], a = [s.a for s in stats], b = [s.b for s in stats],
                    c_hofton = [s.c for s in stats], T_yr = [s.T / yr for s in stats],
                    tau_R_yr = [s.tauR / yr for s in stats], chi2 = [s.chi2 for s in stats]))
print(read("reports/krafla_eve_alt_report.txt", String))

# ── Figure ───────────────────────────────────────────────────────────────────
fig = Figure(size = (2600, 2000), fontsize = 36)
x_km = range(-180, 180, length = 721); x_m = collect(x_km .* 1e3)
panels = [(fig[1, 1], "Foulger et al. (1992), τ = 12–15 yr", xf, uf, sf, tf, false),
          (fig[1, 2], "Hofton & Foulger (1996), τ = 12–17 yr", xh, uh, sh, th, true)]
axs = map(panels) do (pos, title, xo, uo, so, tt, offset)
    ax = Axis(pos, title = title, xlabel = "Distance from Krafla volcanic system (km)",
              ylabel = "Displacement (cm)", limits = (-180, 180, -15, 22))
    for (m, s) in zip(models, stats)
        lines!(ax, x_km, (m.pred(x_m, tt..., m.a, m.b) .+ (offset ? s.c : 0.0)) .* 100,
               color = m.color, linestyle = m.style, linewidth = 4,
               label = @sprintf("%s  (χ² = %.1f)", m.name, s.chi2))
    end
    errorbars!(ax, xo ./ 1e3, uo .* 100, so .* 100, color = :black)
    scatter!(ax, xo ./ 1e3, uo .* 100, color = :black, markersize = 12)
    ax
end
axislegend(axs[1], position = :rb, labelsize = 28)

# (c) zoom near the axis, Hofton profile, where the three forms differ most
ax_c = Axis(fig[2, 1], title = "Hofton & Foulger (1996), near the axis", xlabel = "Distance (km)",
            ylabel = "Displacement (cm)", limits = (-40, 40, -8, 22))
xz = collect(range(-40e3, 40e3, length = 801))
for (m, s) in zip(models, stats)
    lines!(ax_c, xz ./ 1e3, (m.pred(xz, th..., m.a, m.b) .+ s.c) .* 100, color = m.color, linestyle = m.style, linewidth = 4)
end
errorbars!(ax_c, xh ./ 1e3, uh .* 100, sh .* 100, color = :black)
scatter!(ax_c, xh ./ 1e3, uh .* 100, color = :black, markersize = 12)

ax_d = Axis(fig[2, 2], xscale = log10, xtickformat = vs -> [@sprintf("%.1e", v) for v in vs],
            xlabel = L"a\ \ (\mathrm{m\,s^{-1}})", ylabel = L"b\ \ (\mathrm{m})", title = "EVE RSS (alternative), Δχ²")
cmap = cgrad(ColorSchemes.bamako100, 20, categorical = true)
hm = heatmap!(ax_d, a_vec, b_vec, G .- cmin, colormap = cmap, colorrange = (0, DCHI2[3] * 1.4), highclip = cmap[end])
for (lev, c) in zip(DCHI2, [:white, :gold, :tomato])
    contour!(ax_d, a_vec, b_vec, G .- cmin, levels = [lev], color = c, linewidth = 3)
end
scatter!(ax_d, [a_alt], [b_alt], color = :magenta, marker = :xcross, markersize = 26)
Colorbar(fig[2, 3], hm, label = L"\Delta\chi^2", ticks = ([0; DCHI2], ["0", "1σ", "2σ", "3σ"]))
for (ax, l) in zip((axs..., ax_c, ax_d), ("(a)", "(b)", "(c)", "(d)"))
    text!(ax, 0.02, 0.98, text = l, space = :relative, align = (:left, :top), fontsize = 44)
end
save("fig/fig_krafla_eve_alt_compare.png", fig, px_per_unit = 2)
println("Saved fig/fig_krafla_eve_alt_compare.png")
