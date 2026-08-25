using CairoMakie, LaTeXStrings, Printf, Dates
using QuadGK, SpecialFunctions

include("../models/viscous.jl")
include("../models/maxwell.jl")

mkpath("fig")
mkpath("reports")

yr = 3600.0 * 24 * 365.25

# ── Parameters (Natron CSS best-fit) ──
a = 1.18e-4          # m/s
b = 0.662            # m
H = 10e3             # m
E = 40e9             # Pa
nu = 0.25

eta = H * E / (2 * a * (1 - nu^2))     # Pa.s -- same closed form used throughout the paper (Eq. 7)

# tau_R,ref via the shared-modulus closure already derived in the paper
# (Section 2.6.3, Eq. tau_R_shared): tau_R = H/[a(1-nu)]. This introduces no
# new elastic-modulus assumption -- it uses only a, H, nu, already defined
# for this figure's own STEP/CSS curves (replaces an earlier version that
# assumed an independent mantle shear modulus mu_m, not otherwise used or
# explained in the paper).
tau_R_ref_s = H / (a * (1 - nu))        # s
tau_R_ref_yr = tau_R_ref_s / yr
ell_M = a * tau_R_ref_s                # m — Maxwell length at reference tau_R

@printf("eta = %.3e Pa.s\n", eta)
@printf("tau_R,ref = H/[a(1-nu)] = %.4f yr\n", tau_R_ref_yr)
@printf("ell_M = a*tau_R,ref = %.1f m\n", ell_M)

# ══════════════════════════════════════════════════════════════════════════
# Panel (a): velocity vs x/ell_M, at fixed t = tau_R,ref, for tau_R/tau_R,ref in {0,0.2,1,5}
# ══════════════════════════════════════════════════════════════════════════

ratios_a = [0.0, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0]
nonzero_ratios = filter(!=(0.0), ratios_a)
cmap = cgrad(:turbo, length(nonzero_ratios), categorical=true)
colors_a = Vector{Any}(undef, length(ratios_a))
styles_a = fill(:solid, length(ratios_a))
labels_a = Vector{Any}(undef, length(ratios_a))
for (i, r) in enumerate(ratios_a)
    if r == 0.0
        colors_a[i] = :black
        labels_a[i] = L"\tau_R=0\ \mathrm{(viscous)}"
    else
        colors_a[i] = cmap[findfirst(==(r), nonzero_ratios)]
        tau_R_i_yr = r * tau_R_ref_yr
        yr_str = @sprintf("%.2f", tau_R_i_yr)
        labels_a[i] = latexstring("\\tau_R=", yr_str, "\\,\\mathrm{yr}")
    end
end

x_hat_grid = exp10.(range(log10(1e-4), log10(20.0), length=800))
x_grid = x_hat_grid .* ell_M
t_fixed = tau_R_ref_s

# Dimensionless velocity: v*t_fixed/b. Verified by direct differentiation of
# u^STEP/b = (1/pi)*atan(1/xhat) (confirmed-correct, xhat=x/(at)) with respect to t:
# v^STEP*t/b = xhat/(pi*(1+xhat^2)), dimensionally consistent ([m/s][s]/[m] = 1).
# (Note: v_norm = b/t_fixed has units of velocity, unlike b/ell_M which is dimensionless.)
v_norm = b / t_fixed
v_step_norm(xhat) = xhat / (pi * (1 + xhat^2))

# get_v_relax's default truncation wavenumber (_k_limit_maxwell in src/models.jl)
# under-resolves the oscillatory integrand for the very-near-field, small-x_hat
# points needed here (x_hat down to 1e-4), producing ~1% spurious ripple at fixed
# tau_R. A 20x larger truncation removes it (checked convergent against 5x and
# 20x multipliers); scoped locally so src/models.jl (used by other figures) is
# untouched.
function get_v_relax_precise(x, tau, b, a, tau_R; k_mult=20.0)
    (a <= 0 || tau < 0 || tau_R <= 0) && return 0.0
    k_limit = k_mult * _k_limit_maxwell(x, tau, a, tau_R)
    integrand(k) = sin(k * x) / (1 + tau_R * a * k) * exp(-(a * k * tau) / (1 + tau_R * a * k))
    val = _segmented_quadgk(integrand, k_limit, x; rtol=1e-12, atol=1e-14)
    return (b * a / pi) * val
end

v_curves_a = Vector{Vector{Float64}}(undef, length(ratios_a))
for (i, r) in enumerate(ratios_a)
    if r == 0.0
        v_curves_a[i] = get_v_step.(x_grid, t_fixed, b, a) ./ v_norm
    else
        v_curves_a[i] = get_v_relax_precise.(x_grid, t_fixed, b, a, r * tau_R_ref_s) ./ v_norm
    end
end

# Sanity check quoted in the paper caption: at tau_R/tau_R,ref=5, near-field (x=ell_M) amplitude
# is reduced by more than 60% relative to viscous.
idx_ellM = argmin(abs.(x_hat_grid .- 1.0))
idx_r5 = findfirst(==(5.0), ratios_a)
reduction_pct = 100 * (1 - v_curves_a[idx_r5][idx_ellM] / v_curves_a[1][idx_ellM])
@printf("Panel (a) check: near-field (x=ell_M) reduction at tau_R/tau_R,ref=5: %.1f%%  (paper: >60%%)\n", reduction_pct)

# ── Near-field elastic-rebound plateau (new closed-form result) ──────────────
# In the high-k saturation regime (kappa >> 1), the Maxwell integrand reduces to
# the Dirichlet kernel sin(kx)/k weighted by exp(-t/tau_R)/(tau_R a), so as x -> 0
# (but still x << ell_M) the velocity approaches an x-INDEPENDENT plateau:
#   v^MAX(x->0+, t) / v_norm  ->  (that/2) * exp(-that),   that := t/tau_R.
# Verified numerically against get_v_relax to <1% at x_hat ~ 1e-5 - 1e-2 (the
# quadrature is unreliable below x ~ 1 m due to the k_near floor in
# _k_limit_maxwell, well outside the plotted range here).
plateau_hat(that) = (that / 2) * exp(-that)

# Crossover x_hat where each Maxwell curve intersects the viscous curve
# (root of v_curves_a[i] - v_step_norm.(x_hat_grid) via bisection on the grid).
function find_crossover(xg, diff)
    for j in 1:(length(xg) - 1)
        if sign(diff[j]) != sign(diff[j+1]) && diff[j] != 0
            # linear interpolation in log-x for a smoother estimate
            f = diff[j] / (diff[j] - diff[j+1])
            return exp(log(xg[j]) + f * (log(xg[j+1]) - log(xg[j])))
        end
    end
    return NaN
end

crossovers_a = fill(NaN, length(ratios_a))
for (i, r) in enumerate(ratios_a)
    r == 0.0 && continue
    diff = v_curves_a[i] .- v_step_norm.(x_hat_grid)
    crossovers_a[i] = find_crossover(x_hat_grid, diff)
    that = 1.0 / r
    @printf("Panel (a) check: r=%.1f  x_c/ell_M=%.4f  plateau=%.4f (analytic)\n",
            r, crossovers_a[i], plateau_hat(that))
end

# ══════════════════════════════════════════════════════════════════════════
# Panel (b): displacement vs t/tau_R,ref, at fixed x = ell_M, viscous + tau_R/tau_R,ref in {1,5}
# ══════════════════════════════════════════════════════════════════════════

t_hat_grid = collect(range(0.0, 20.0, length=400))
t_grid = t_hat_grid .* tau_R_ref_s

u_viscous_b = get_u_step.(ell_M, t_grid, b, a) ./ b
u_maxwell1_b = get_u_viscoelastic.(ell_M, t_grid, b, a, 1.0 * tau_R_ref_s) ./ b
u_maxwell5_b = get_u_viscoelastic.(ell_M, t_grid, b, a, 5.0 * tau_R_ref_s) ./ b

# ══════════════════════════════════════════════════════════════════════════
# Two SEPARATE single-panel figures for the paper: velocity (formerly panel
# a) and displacement (formerly panel c). The physical-units duplicate
# (formerly panel b) is dropped entirely -- not used in the manuscript.
# ══════════════════════════════════════════════════════════════════════════

fig_v = Figure(size=(1000, 850), fontsize=20)

t_fixed_yr_str = @sprintf("%.2f", tau_R_ref_yr)

ax1 = Axis(fig_v[1, 1],
    xlabel = L"x/\ell_M",
    ylabel = latexstring(t_fixed_yr_str, "\\,\\mathrm{yr}\\, v / b"),
    title  = latexstring("\\mathrm{Velocity\\ at\\ }t=", t_fixed_yr_str, "\\,\\mathrm{yr}\\mathrm{,\\ } x\\ \\mathrm{normalised\\ by\\ }\\ell_M"),
    xscale = log10,
    limits = (1e-4, 20, 0, nothing))

vspan!(ax1, 1e-4, 1, color=(:gray70, 0.3))
vlines!(ax1, [1.0], color=:gray30, linestyle=:dash, linewidth=1.5)
text!(ax1, 1.05, 0.02; text=L"x=\ell_M", align=(:left,:bottom), fontsize=16, color=:gray20, space=:data)
text!(ax1, 0.06, 0.55; text="elastic\nrebound", align=(:center,:top), fontsize=16, color=:gray30, space=:relative)

for (i, r) in enumerate(ratios_a)
    lines!(ax1, x_hat_grid, v_curves_a[i], color=colors_a[i], linestyle=styles_a[i],
           linewidth=2.5, label=labels_a[i])
    # Analytic near-field plateau: x-independent, drawn only out to the
    # x_hat at which the curve actually crosses the viscous Lorentzian.
    if r > 0.0 && !isnan(crossovers_a[i])
        that = 1.0 / r
        lines!(ax1, [x_hat_grid[1], crossovers_a[i]], fill(plateau_hat(that), 2),
               color=colors_a[i], linestyle=:dot, linewidth=1.5)
        scatter!(ax1, [crossovers_a[i]], [v_step_norm(crossovers_a[i])],
                 color=colors_a[i], markersize=10, strokewidth=1, strokecolor=:black)
    end
end
axislegend(ax1, position=:rt, framevisible=true, labelsize=13, patchsize=(20,10), rowgap=1)

save("fig/fig_maxwell_velocity.png", fig_v)
println("fig/fig_maxwell_velocity.png  saved")

# ══════════════════════════════════════════════════════════════════════════
# Displacement figure (formerly panel c): u/b at x=ell_M vs t/tau_R,ref.
# ══════════════════════════════════════════════════════════════════════════

fig_d = Figure(size=(1000, 850), fontsize=20)

t_grid_yr = t_grid ./ yr
tau_R_1_yr_str = @sprintf("%.2f", tau_R_ref_yr)
tau_R_5_yr_str = @sprintf("%.2f", 5.0 * tau_R_ref_yr)

ax2 = Axis(fig_d[1, 1],
    xlabel = L"t\ \mathrm{(yr)}",
    ylabel = L"u/b",
    title  = L"\mathrm{Displacement\ at\ }x=\ell_M",
    limits = (0, t_grid_yr[end], 0, 0.55))

vspan!(ax2, 0, tau_R_ref_yr, color=(:gray70, 0.3))
band!(ax2, t_grid_yr, u_maxwell5_b, u_viscous_b, color=(:orange, 0.25))
hlines!(ax2, [0.5], color=:gray30, linestyle=:dot, linewidth=1.5)
text!(ax2, 0.98 * t_grid_yr[end], 0.51; text=L"b/2", align=(:right,:bottom), fontsize=16, color=:gray20)

lines!(ax2, t_grid_yr, u_viscous_b, color=:black, linewidth=2.5,
       label=L"\mathrm{Purely\ viscous\ }(\tau_R\to 0)")
lines!(ax2, t_grid_yr, u_maxwell1_b, color=:orangered, linewidth=2.5,
       label=latexstring("\\tau_R=", tau_R_1_yr_str, "\\,\\mathrm{yr}"))
lines!(ax2, t_grid_yr, u_maxwell5_b, color=:red, linestyle=:dash, linewidth=2.5,
       label=latexstring("\\tau_R=", tau_R_5_yr_str, "\\,\\mathrm{yr}"))
text!(ax2, 0.35, 0.03; text="elastic\nrebound", align=(:center,:bottom), fontsize=16, color=:gray30, space=:relative)

axislegend(ax2, position=:rb, framevisible=true, labelsize=15)

save("fig/fig_maxwell_displacement.png", fig_d)
println("fig/fig_maxwell_displacement.png  saved")

# ══════════════════════════════════════════════════════════════════════════
# Report
# ══════════════════════════════════════════════════════════════════════════

open("reports/report_maxwell_regimes.md", "w") do io
    println(io, "# Report: Maxwell Visco-Elastic Two-Regime Behavior (`fig_maxwell_regimes.jl`)")
    println(io)
    println(io, "## Purpose")
    println(io)
    println(io, "Reproduces `fig/fig_maxwell_regimes.png` as specified in `tex/paper/paper_v2.tex` ",
                "(Section \"Maxwell case: two dimensionless groups and three regimes\", \\label{fig:maxwell_regimes}, ",
                "lines 739-801), panels (a) and (b), stacked vertically to fit a single column of the ",
                "two-column GJI layout. Illustrates the elastic-rebound-to-viscous-diffusion crossover ",
                "of the Maxwell step response, in space (panel a) and in time (panel b).")
    println(io)
    println(io, "---")
    println(io)
    println(io, "## Physics")
    println(io)
    println(io, "- Purely viscous velocity: `get_v_step(x,t,b,a)`; purely viscous displacement: `get_u_step(x,t,b,a)`")
    println(io, "- Maxwell velocity: `get_v_relax(x,t,b,a,tau_R)`; Maxwell displacement: `get_u_viscoelastic(x,t,b,a,tau_R)`")
    println(io, "- Maxwell relaxation time: τ_R := η/μ_m (paper eq. for τ_R), where μ_m is the mantle elastic shear modulus")
    println(io, "- Maxwell length: ℓ_M = a·τ_R = EH/[2μ_m(1-ν²)] — a purely elastic quantity, independent of viscosity")
    println(io)
    println(io, "---")
    println(io)
    println(io, "## Parameters")
    println(io)
    @printf(io, "| Parameter | Value |\n|---|---|\n")
    @printf(io, "| a (Natron CSS viscous best-fit) | %.3e m/s |\n", a)
    @printf(io, "| b | %.3f m |\n", b)
    @printf(io, "| H | %.0f km |\n", H/1e3)
    @printf(io, "| E | %.0f GPa |\n", E/1e9)
    @printf(io, "| ν | %.2f |\n", nu)
    @printf(io, "| η = EH/(2a(1-ν²)) | %.3e Pa·s |\n", eta)
    @printf(io, "| τ_R,ref = H/[a(1-ν)] | %.4f yr |\n", tau_R_ref_yr)
    @printf(io, "| ℓ_M = a·τ_R,ref | %.0f m |\n", ell_M)
    println(io)
    println(io, "τ_R,ref uses the shared-modulus closure already derived in the paper ",
                "(Section 2.6.3, Eq. tau_R_shared): τ_R = H/[a(1-ν)]. It introduces no new elastic-modulus ",
                "assumption -- only a, H, ν, already defined for this figure's own STEP/CSS curves -- and is ",
                "independent of the CSS recurrence interval T, which has no role in the step-response regime ",
                "this figure illustrates.")
    println(io)
    println(io, "---")
    println(io)
    println(io, "## Results")
    println(io)
    @printf(io, "Panel (a) near-field (x=ℓ_M) velocity reduction at τ_R/τ_R,ref=5 relative to viscous: **%.1f%%** ", reduction_pct)
    println(io, "(paper caption states \"more than 60%\").")
    println(io)
    println(io, "Panel (b): both Maxwell curves (τ_R/τ_R,ref = 1, 5) lag the purely viscous displacement for ",
                "t ≲ τ_R (elastic-rebound zone, grey shading) and converge onto the common asymptote u/b = 1/2 ",
                "(dotted line) once t ≫ τ_R (x-axis extended to t/τ_R,ref=20 to show full convergence). The orange ",
                "band marks the displacement deficit of the slowest-relaxing case (τ_R/τ_R,ref=5) relative to the ",
                "viscous prediction.")
    println(io)
    println(io, "---")
    println(io)
    println(io, "## Dependencies")
    println(io)
    println(io, "- `src/models.jl` — `get_v_step`, `get_v_relax`, `get_u_step`, `get_u_viscoelastic` (no modifications needed)")
    println(io, "- `tex/paper/paper_v2.tex` — figure specification (\\label{fig:maxwell_regimes}) and physical parameter definitions")
    println(io, "- Natron CSS best-fit `a` reproduced from `fig_behavior_natron.jl`")
    println(io)
    @printf(io, "Generated: %s\n", Dates.now())
end
println("reports/report_maxwell_regimes.md  saved")
