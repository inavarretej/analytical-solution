using CairoMakie, LaTeXStrings, Statistics, Random, SpecialFunctions, QuadGK

include("../models/viscous.jl")
include("../models/maxwell.jl")

# EV RSS velocity profile under stochastic recurrence (fig:stochastic_css).
# Panels (2026-10-21; T = b/V holds on average in every case):
#   (a) contrast: b_j random, T fixed (opening not matched to the interval)
#   (b) coupled:  T_j random, b_j = V T_j (opening = divergence since the
#       previous intrusion)
#   (c) dike position: x_j ~ Normal(0, sigma_x), T and b periodic
#   (d) coupled + dike position
# A uniform x_j distribution with the same standard deviation gives bands
# within 5-10 per cent of the normal case (checked 2026-10-21 with a preview
# script, since removed).

# ── Editable parameters ───────────────────────────────────────────────────────
cv_T     = 0.30   # coefficient of variation for T (σ_T / T_mean)
cv_b     = 0.30   # coefficient of variation for b (σ_b / b_mean)
sigma_x  = 0.10   # standard deviation of the dike position, in units of aT
N_events = 300    # past events per realization — large N reduces finite-sum far-field residual
N_real   = 500    # ensemble size
p_lo     = 10     # lower percentile bound for band
p_hi     = 90     # upper percentile bound for band

# ── Dimensionless units: a = T_mean = b_mean = 1 ─────────────────────────────
a_d    = 1.0
T_d    = 1.0
b_d    = 1.0
vp_d   = b_d / T_d
τ_phase = 0.5 * T_d   # fixed phase: midpoint of mean cycle

# Both sides of the rift; avoid x = 0 (velocity is 0 there by symmetry anyway)
x_hat  = [collect(range(-5.0, -0.02, length = 300)); collect(range(0.02, 5.0, length = 300))]
x_phys = x_hat .* (a_d * T_d)

# ── Deterministic RSS (closed-form Digamma) ───────────────────────────────────
v_det = [get_v_css(xp, τ_phase, b_d, vp_d, a_d) * T_d / b_d for xp in x_phys]

# ── One stochastic realization ────────────────────────────────────────────────
# Events in the past at times t_past[j] = τ_phase + Σ_{k<j} T_k, opened at x_j.
# T_k ~ Normal(T_d, cv_T·T_d), b_k ~ Normal(b_d, cv_b·b_d), both clamped to
# ≥ 5% of their mean; coupled = true sets b_k = (b_d/T_d) T_k; x_k ~
# Normal(0, sx·aT). The raw sum of STEP velocities gives the RSS far-field
# value b/(2T) directly (vault note "Factor-of-2 error in the EV RSS
# velocity formula").
function stoch_profile(cvT, cvb; coupled = false, sx = 0.0)
    T_j    = max.(T_d .+ cvT * T_d .* randn(N_events), 0.05 * T_d)
    b_j    = coupled ? (b_d / T_d) .* T_j : max.(b_d .+ cvb * b_d .* randn(N_events), 0.05 * b_d)
    x_j    = sx > 0 ? sx * a_d * T_d .* randn(N_events) : zeros(N_events)
    t_past = τ_phase .+ [0.0; cumsum(T_j[1:N_events-1])]
    v = zeros(length(x_hat))
    for j in 1:N_events
        @. v += b_j[j] * a_d / π * (x_phys - x_j[j]) / ((x_phys - x_j[j])^2 + (a_d * t_past[j])^2)
    end
    return v .* (T_d / b_d)
end

# ── Ensemble runner ───────────────────────────────────────────────────────────
function run_ensemble(f; seed::Int)
    Random.seed!(seed)
    mat = Matrix{Float64}(undef, N_real, length(x_hat))
    for r in 1:N_real
        mat[r, :] = f()
    end
    lo  = [quantile(view(mat, :, i), p_lo / 100) for i in axes(mat, 2)]
    hi  = [quantile(view(mat, :, i), p_hi / 100) for i in axes(mat, 2)]
    med = [quantile(view(mat, :, i), 0.50)        for i in axes(mat, 2)]
    return lo, hi, med
end

println("Ensemble 1/4: variable b, T fixed (contrast) …")
lo_b, hi_b, med_b = run_ensemble(() -> stoch_profile(0.0, cv_b); seed = 2)
println("Ensemble 2/4: coupled b_j = V T_j …")
lo_C, hi_C, med_C = run_ensemble(() -> stoch_profile(cv_T, 0.0; coupled = true); seed = 4)
println("Ensemble 3/4: dike position …")
lo_X, hi_X, med_X = run_ensemble(() -> stoch_profile(0.0, 0.0; sx = sigma_x); seed = 5)
println("Ensemble 4/4: coupled + dike position …")
lo_CX, hi_CX, med_CX = run_ensemble(() -> stoch_profile(cv_T, 0.0; coupled = true, sx = sigma_x); seed = 6)
println("Done.")

# Band widths in units of the far-field velocity V/2 (= 0.5 here), and the
# largest departure of the median from the deterministic RSS
for (name, lo, hi, med) in [("(a) variable b", lo_b, hi_b, med_b), ("(b) coupled", lo_C, hi_C, med_C),
                            ("(c) position", lo_X, hi_X, med_X), ("(d) coupled + position", lo_CX, hi_CX, med_CX)]
    print(rpad(name, 26), " band/(V/2):")
    for ξ in (0.02, 0.1, 0.3, 1.0, 2.0, 5.0)
        i = argmin(abs.(x_hat .- ξ))
        print("  ξ=", round(x_hat[i], digits = 2), ": ", round((hi[i] - lo[i]) / 0.5, digits = 3))
    end
    println("   max|median − RSS|/(V/2) = ", round(maximum(abs.(med .- v_det)) / 0.5, digits = 3))
end

# ── Figure (single column, 4 rows) ───────────────────────────────────────────
fig = Figure(size = (700, 1550), fontsize = 26.25)

cases = [
    (1, lo_b,  hi_b,  med_b,  :darkorange,  "(a)  Variable b, T fixed   (cv = $(cv_b))", false),
    (2, lo_C,  hi_C,  med_C,  :purple,      "(b)  Coupled   (b_j = V T_j,  cv_T = $(cv_T))", false),
    (3, lo_X,  hi_X,  med_X,  :firebrick,   "(c)  Dike position   (σ_x/aT = $(sigma_x))", true),
    (4, lo_CX, hi_CX, med_CX, :teal,        "(d)  Coupled + dike position", true),
]

for (row, lo, hi, med, color, title, inset) in cases
    ax = Axis(fig[row, 1],
        xlabel = row == 4 ? L"\xi = x/(aT)" : "",
        ylabel = L"v \cdot T / b",
        title  = title,
        limits = (-5, 5, -0.95, 1.1))

    vlines!(ax, [0.0],      color = (:black, 0.20), linewidth = 1.0, linestyle = :dot)
    vlines!(ax, [-1.0, 1.0], color = (:black, 0.15), linewidth = 1.0, linestyle = :dash)
    hlines!(ax, [0.0],      color = (:black, 0.15), linewidth = 0.8)

    band!(ax, x_hat, lo, hi, color = (color, 0.25))
    lines!(ax, x_hat, med,   color = color,  linewidth = 2.0,
           label = "Stochastic median  ($(p_lo)–$(p_hi) pct)")
    lines!(ax, x_hat, v_det, color = :black, linewidth = 2.5,
           linestyle = :dash, label = "Deterministic RSS")

    axislegend(ax, position = :lt, framevisible = true, labelsize = 23.16)

    if inset   # zoom on |ξ| ≤ 0.5, where dike-position scatter acts
        ax_in = Axis(fig[row, 1], width = Relative(0.36), height = Relative(0.42),
                     halign = 0.97, valign = 0.12, limits = (-0.5, 0.5, -0.6, 0.6),
                     xticklabelsize = 15, yticklabelsize = 15, xticks = [-0.5, 0, 0.5], yticks = [-0.5, 0, 0.5],
                     backgroundcolor = :white)
        translate!(ax_in.scene, 0, 0, 10)
        sel = abs.(x_hat) .≤ 0.5
        band!(ax_in, x_hat[sel], lo[sel], hi[sel], color = (color, 0.25))
        lines!(ax_in, x_hat[sel], med[sel], color = color, linewidth = 1.5)
        lines!(ax_in, x_hat[sel], v_det[sel], color = :black, linewidth = 1.8, linestyle = :dash)
    end
end

mkpath("fig")
save("fig/fig_stochastic_css.png", fig)
println("fig/fig_stochastic_css.png  saved")
