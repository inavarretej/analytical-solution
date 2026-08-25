using CairoMakie, LaTeXStrings, Statistics, Random, SpecialFunctions, QuadGK

include("../models/viscous.jl")
include("../models/maxwell.jl")

# ── Editable parameters ───────────────────────────────────────────────────────
cv_T     = 0.30   # coefficient of variation for T (σ_T / T_mean)
cv_ω     = 0.30   # coefficient of variation for ω (σ_ω / ω_mean)
N_events = 300    # past events per realization — large N reduces finite-sum far-field residual
N_real   = 500    # ensemble size
p_lo     = 10     # lower percentile bound for band
p_hi     = 90     # upper percentile bound for band

# ── Dimensionless units: a = T_mean = ω_mean = 1 ─────────────────────────────
a_d    = 1.0
T_d    = 1.0
ω_d    = 1.0
vp_d   = ω_d / T_d
τ_phase = 0.5 * T_d   # fixed phase: midpoint of mean cycle

# Both sides of the rift; avoid x = 0 (velocity is 0 there by symmetry anyway)
x_hat  = [collect(range(-5.0, -0.02, length = 300)); collect(range(0.02, 5.0, length = 300))]
x_phys = x_hat .* (a_d * T_d)

# ── Deterministic CSS (closed-form Digamma) ───────────────────────────────────
v_det = [get_v_css(xp, τ_phase, ω_d, vp_d, a_d) * T_d / ω_d for xp in x_phys]

# ── One stochastic realization ────────────────────────────────────────────────
# Events in the past at times t_past[j] = τ_phase + Σ_{k<j} T_k,
# where T_k ~ Normal(T_d, cv_T·T_d) and ω_k ~ Normal(ω_d, cv_ω·ω_d),
# both clamped to ≥ 5% of their mean to stay positive.
function stoch_profile(cv_T_val::Float64, cv_ω_val::Float64)
    T_j    = max.(T_d .+ cv_T_val * T_d .* randn(N_events), 0.05 * T_d)
    ω_j    = max.(ω_d .+ cv_ω_val * ω_d .* randn(N_events), 0.05 * ω_d)
    t_past = τ_phase .+ [0.0; cumsum(T_j[1:N_events-1])]

    v = zeros(length(x_hat))
    for j in 1:N_events
        t_j = t_past[j]
        # viscous step velocity: ω·a/π · x / (x² + (a·t)²)
        @. v += ω_j[j] * a_d / π * x_phys / (x_phys^2 + (a_d * t_j)^2)
    end
    # factor 2: v_CSS = 2 × Σ v_step (CSS uses 2b/πT prefactor vs b/πT for one step)
    return 2 .* v .* (T_d / ω_d)
end

# ── Coupled profile: ω_j = v_plate · T_j for every event ─────────────────────
# Far-field plate velocity is preserved exactly in every realisation.
function stoch_coupled(cv_T_val::Float64)
    T_j    = max.(T_d .+ cv_T_val * T_d .* randn(N_events), 0.05 * T_d)
    ω_j    = (ω_d / T_d) .* T_j          # coupling: ω ∝ T for every event
    t_past = τ_phase .+ [0.0; cumsum(T_j[1:N_events-1])]

    v = zeros(length(x_hat))
    for j in 1:N_events
        t_j = t_past[j]
        @. v += ω_j[j] * a_d / π * x_phys / (x_phys^2 + (a_d * t_j)^2)
    end
    return 2 .* v .* (T_d / ω_d)
end

# ── Ensemble runner ───────────────────────────────────────────────────────────
function run_ensemble(cv_T_val::Float64, cv_ω_val::Float64; seed::Int = 1)
    Random.seed!(seed)
    mat = Matrix{Float64}(undef, N_real, length(x_hat))
    for r in 1:N_real
        mat[r, :] = stoch_profile(cv_T_val, cv_ω_val)
    end
    lo  = [quantile(view(mat, :, i), p_lo / 100) for i in axes(mat, 2)]
    hi  = [quantile(view(mat, :, i), p_hi / 100) for i in axes(mat, 2)]
    med = [quantile(view(mat, :, i), 0.50)        for i in axes(mat, 2)]
    return lo, hi, med
end

println("Ensemble 1/4: variable T …")
lo_T, hi_T, med_T = run_ensemble(cv_T, 0.0;  seed = 1)
println("Ensemble 2/4: variable ω …")
lo_ω, hi_ω, med_ω = run_ensemble(0.0,  cv_ω; seed = 2)
println("Ensemble 3/4: both variable …")
lo_B, hi_B, med_B = run_ensemble(cv_T, cv_ω; seed = 3)
println("Ensemble 4/4: coupled ω = v_plate·T …")
Random.seed!(4)
mat_C = Matrix{Float64}(undef, N_real, length(x_hat))
for r in 1:N_real; mat_C[r, :] = stoch_coupled(cv_T); end
lo_C = [quantile(view(mat_C, :, i), p_lo / 100) for i in axes(mat_C, 2)]
hi_C = [quantile(view(mat_C, :, i), p_hi / 100) for i in axes(mat_C, 2)]
med_C = [quantile(view(mat_C, :, i), 0.50)       for i in axes(mat_C, 2)]
println("Done.")

# ── Figure (single column, 4 rows) — fits one column of the two-column
# GJI layout without squeezing any panel below legibility ───────────────────
fig = Figure(size = (700, 1550), fontsize = 17)

cases = [
    (1, lo_T, hi_T, med_T, :royalblue,
     "(a)  Variable T   (cv = $(cv_T),  ω fixed)"),
    (2, lo_ω, hi_ω, med_ω, :darkorange,
     "(b)  Variable ω   (T fixed,  cv = $(cv_ω))"),
    (3, lo_B, hi_B, med_B, :forestgreen,
     "(c)  Both variable   (cv_T = $(cv_T),  cv_ω = $(cv_ω))"),
    (4, lo_C, hi_C, med_C, :purple,
     "(d)  Coupled   (ω_j = v_plate · T_j,  cv_T = $(cv_T))"),
]

for (row, lo, hi, med, color, title) in cases
    ax = Axis(fig[row, 1],
        xlabel = row == 4 ? L"\hat{x} = x/(aT)" : "",
        ylabel = L"v \cdot T / \omega",
        title  = title,
        limits = (-5, 5, nothing, nothing))

    vlines!(ax, [0.0],      color = (:black, 0.20), linewidth = 1.0, linestyle = :dot)
    vlines!(ax, [-1.0, 1.0], color = (:black, 0.15), linewidth = 1.0, linestyle = :dash)
    hlines!(ax, [0.0],      color = (:black, 0.15), linewidth = 0.8)
    text!(ax,  1.08, -0.88; text = L"\hat{x}=1",  fontsize = 13, color = (:black, 0.45))
    text!(ax, -1.08, -0.88; text = L"\hat{x}=-1", fontsize = 13, color = (:black, 0.45),
          align = (:right, :center))

    band!(ax, x_hat, lo, hi, color = (color, 0.25))
    lines!(ax, x_hat, med,   color = color,  linewidth = 2.0,
           label = "Stochastic median  ($(p_lo)–$(p_hi) pct)")
    lines!(ax, x_hat, v_det, color = :black, linewidth = 2.5,
           linestyle = :dash, label = "Deterministic CSS")

    axislegend(ax, position = :rb, framevisible = true, labelsize = 15)
end

mkpath("fig")
save("fig/fig_stochastic_css.png", fig)
println("fig/fig_stochastic_css.png  saved")
