# An analytical recurrent-steady-state model for inter-rifting deformation at magmatic rifts

Iván Navarrete and Jean-Arthur Olive — Laboratoire de Géologie, École Normale Supérieure, PSL Research University, Paris.

**Manuscript:** [`tex/paper/paper.pdf`](tex/paper/paper.pdf) (LaTeX source: [`tex/paper/paper.tex`](tex/paper/paper.tex), Geophysical Journal International format).

## Summary

Inter-rifting deformation, the strain that accumulates between discrete magmatic episodes, is the least-constrained phase of the rifting cycle. This repository contains the analytical model, the fits and the figures of a paper that treats it as the superposed relaxation of a long history of discrete dike intrusions in an elastic plate overlying a viscous half-space.

- **EV model (elastic plate over a Newtonian viscous half-space).** The non-local Stokes coupling between the plate and the half-space turns post-diking relaxation into a fractional diffusion governed by the half-Laplacian and a single characteristic velocity `a = EH/[2η(1−ν²)]`. A single intrusion gives a closed-form step response (`u = (b/π) arctan(at/x)`), which spreads at the constant speed `a` with a `1/x` tail. Periodic intrusions converge to a **Recurrent Steady State (RSS)** with a closed form in the Digamma and log-Gamma functions. The RSS conserves the opening of each intrusion over one cycle, so far-field plate divergence emerges without boundary conditions, and its only length scale is the deformation half-width `aT`.
- **EVE model (elastic plate over a Maxwell visco-elastic half-space).** Derived through the elastic–visco-elastic correspondence principle, with the Maxwell time tied to `a` by a shared-modulus closure, `τ_R = H/[a(1−ν)]`. It departs from the EV model only near the axis and during the first few `τ_R`.
- **Applications.** Both RSS models are fitted to post-rifting geodetic records from **Krafla, Iceland** (two displacement profiles, Foulger et al. 1992; Hofton & Foulger 1996) and **Asal, Djibouti** (EP00–DF00 baseline time series and velocity profile, Vigny et al. 2007). The recurrence interval is never a free parameter: it is derived as `T = b/V`, from the fitted dike opening `b` and an independent long-term rift opening rate `V`.
- **Irregular recurrence.** Stochastic ensembles test the RSS against scatter in recurrence interval, dike opening and dike position.

## Repository layout

```text
Project.toml, Manifest.toml   Julia environment
src/
  models/                     forward models (EV, EVE, alternative EVE)
  parameters.jl               reference elastic and viscous parameters
  halfspace.jl                interior velocity field of the half-space (Poisson kernel)
  figures/                    scripts for the manuscript figures and the fits
  checks/                     numerical checks of the alternative EVE initial condition
  annex/                      numerical checks of the fractional-Laplacian identities
data/                         digitized geodetic data (not tracked by git)
fig/                          figures written by the scripts (not tracked by git)
reports/                      fit reports and best-fit parameters (not tracked by git)
cache/                        checkpoints of the expensive EVE contour figures (not tracked by git)
notebooks/                    Pluto notebook on the EVE RSS velocity
tex/paper/                    manuscript (paper.tex, refs.bib, GJI class files)
matlab/                       earlier MATLAB version of the model (not used)
```

`data/`, `fig/`, `reports/` and `cache/` are listed in `.gitignore`. A fresh clone therefore needs the data files below before the fit scripts can run.

## Getting started

Julia 1.12 is used. From the repository root:

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'   # install the environment once
julia --project=. src/figures/<script>.jl             # run a script
julia -t auto --project=. src/figures/<script>.jl     # multithreaded (landscapes, contour figures)
```

All scripts are run from the repository root and write to `fig/` and `reports/`.

To compile the manuscript:

```bash
cd tex/paper && latexmk -pdf paper.tex
```

### Data files (`data/`)

| File | Content | Used by |
| --- | --- | --- |
| `extracted_data.csv` | Krafla, Foulger et al. (1992) Fig. 2b, 20 stations, 1987–1990 | Krafla scripts |
| `hofton.csv` | Krafla, Hofton & Foulger (1996) Fig. 5, 17 stations, 1987–1992 | Krafla scripts |
| `asal_timeseries_merged.csv` | Asal, EP00–DF00 baseline change, Vigny et al. (2007) Fig. 3, 33 points | Asal scripts |
| `asal_profile.csv` | Asal, rift-perpendicular velocity profile, Vigny et al. (2007), 15 points | Asal scripts |

The other files in `data/` (Natron velocities and faults, earlier digitizations) are not used by the manuscript.

## Models (`src/models/`)

| File | Content |
| --- | --- |
| `viscous.jl` | EV model: step response `get_u_step`, `get_v_step`; RSS displacement and velocity in closed form, `get_u_css`, `get_v_css` (Digamma / log-Gamma); Foulger et al. (1992) diffusion model, `get_u_foulger`, `get_v_foulger`, written with the full dike opening `b` (`U₀ = b/2`); differenced displacements between two epochs. |
| `maxwell.jl` | EVE model: step response by numerical quadrature (`get_u_viscoelastic`, `get_v_relax`) and its `τ → 0⁺` limit; EVE RSS in closed periodic form (`get_*_css_maxwell_closed`) and in the "split" form used by the fits (`get_*_css_maxwell_split`: the large-`k` limit is integrated in closed form and the remainder numerically). |
| `maxwell_alt.jl` | Alternative EVE RSS in which the state just after each intrusion is obtained from the visco-elastic kernel integrated across the intrusion. Used only by the two `*_eve_alt_compare.jl` scripts; the formulation is under discussion and is not used in the manuscript. |

`src/parameters.jl` sets the reference plate parameters (`H = 10 km`, `E = 40 GPa`, `ν = 0.25`, shear modulus 16 GPa) and a reference viscosity used by the model-behaviour figures. `src/halfspace.jl` computes the interior velocity field from the surface velocity (Poisson-kernel convolution).

## Scripts (`src/figures/`)

### Manuscript figures

| Script | Output (`fig/`) | Manuscript | Notes |
| --- | --- | --- | --- |
| `krafla_diffusion_published.jl` | `fig_krafla_diffusion_published.png` | Fig. 1 | Krafla profiles with the diffusion curves published by Foulger et al. (1992) and Hofton & Foulger (1996), reproduced without refitting. |
| — | `EsquemaModelo.png` | Fig. 2 | Model schematic, drawn by hand (no script). |
| `halfspace_step_contour.jl` | `fig_halfspace_step_contour.png` | Fig. 3 | EV step response: surface profiles and interior speed field at six times. Uses `halfspace_common.jl`. |
| `step_selfsimilar_raw.jl` | `fig_step_selfsimilar_raw.png` | Fig. 4 | Self-similar collapse of the EV and Fickian step responses. |
| `rss_accumulation.jl` | `fig_rss_accumulation.png` | Fig. 5 | Cumulative EV RSS displacement over three cycles. |
| `halfspace_css_contour.jl` | `fig_halfspace_css_contour.png` | Fig. 6 | EV RSS: surface profiles and interior field at six phases. Uses `halfspace_common.jl`. |
| `css.jl` | `fig_css_behavior.png` | Fig. 7 | Convergence of the finite sum of past intrusions to the RSS. |
| `halfspace_step_maxwell_contour.jl` | `fig_halfspace_step_maxwell_contour.png` | Fig. 8 | EVE step response, same layout as Fig. 3. **Expensive** (about 11 h with 12 threads); checkpoints in `cache/`, progress in a log file. |
| `halfspace_css_maxwell_contour.jl` | `fig_halfspace_css_maxwell_contour.png` | Fig. 9 | EVE RSS, same layout as Fig. 6. **Expensive** (several hours with threads); checkpoints in `cache/`. |
| `krafla_v_constrained_combined.jl` | `fig_krafla_v_constrained_combined.png` | Fig. 10 | Krafla joint EV and EVE RSS fits and their χ² landscapes. Re-derives the best fits and writes `reports/krafla_v_constrained_combined_params.csv` and `reports/krafla_chi2_landscape_v_constrained.txt`. Needs `reports/krafla_css_v_constrained_params.csv` as seed. |
| `asal_v_constrained_combined.jl` | `fig_asal_v_constrained_combined.png` | Fig. 11 | Asal joint EV and EVE RSS fits (time series and profile) and their χ² landscapes. Reads `reports/asal_css_v_constrained_params.csv`. Uses `asal_v_constrained_common.jl`. |
| `rheology_1sigma_krafla.jl` | `fig_rheology_1sigma.png` | Fig. 12 | Krafla iso-`a` trade-off in (H, η), with literature viscosities. Best-fit `a` and 1σ range written in the script. |
| `rheology_1sigma_asal.jl` | `fig_rheology_1sigma_asal.png` | Fig. 13 | Same for Asal. |
| `css_width.jl` | `fig_css_width.png` | Fig. 14 | Deformation half-width `aT` against `T` for both sites, with network apertures. |
| `stochastic_css.jl` | `fig_stochastic_css.png` | Fig. 15 | RSS under stochastic recurrence: random openings (uncoupled), coupled `b_j = V T_j`, random dike position, and coupled + position. |

Shared include files: `halfspace_common.jl` (grid, phases, plotting routine and cached interior fields for Figs. 3, 6, 8, 9) and `asal_v_constrained_common.jl` (Asal data, forward models, χ² and the `T ≥ 37.2 yr` floor).

### Fits and numbers cited in the text

These scripts produce numbers quoted in the manuscript or files read by the figure scripts. Some also write diagnostic figures that the manuscript does not use.

| Script | Writes | Used for |
| --- | --- | --- |
| `krafla_css_v_constrained.jl` | `reports/krafla_css_v_constrained_params.csv` | Seed for `krafla_v_constrained_combined.jl`. |
| `asal_css_v_constrained.jl` | `reports/asal_css_v_constrained_params.csv`, `reports/asal_css_v_constrained_report.txt` | Asal best fits (Table 1), read by the Asal figure and landscape scripts. |
| `asal_chi2_landscape_v_constrained.jl` | `reports/asal_chi2_landscape_v_constrained.txt` | Asal 1σ ranges (Table 1). |
| `krafla_hofton_diffusion.jl` | console output | Joint refit of the diffusion model at Krafla (Section 4.1 and Table 1: κ = 2.70 m² s⁻¹, b = 2.81 m, c = 6.4 cm, χ²_r = 4.90). |

### Run order for the fits

```bash
# Krafla
julia --project=. src/figures/krafla_css_v_constrained.jl
julia -t auto --project=. src/figures/krafla_v_constrained_combined.jl
# Asal
julia --project=. src/figures/asal_css_v_constrained.jl
julia -t auto --project=. src/figures/asal_chi2_landscape_v_constrained.jl
julia -t auto --project=. src/figures/asal_v_constrained_combined.jl
```

The rheology (Figs. 12–13) and half-width (Fig. 14) scripts take the best-fit values written in their headers; update them if the fits change.

### Alternative EVE initial condition (under discussion, not in the manuscript)

| Script | Writes |
| --- | --- |
| `krafla_eve_alt_compare.jl` | `fig/fig_krafla_eve_alt_compare.png`, `reports/krafla_eve_alt_report.txt`, `reports/krafla_eve_alt_params.csv` |
| `asal_eve_alt_compare.jl` | `fig/fig_asal_eve_alt_compare.png`, `reports/asal_eve_alt_report.txt`, `reports/asal_eve_alt_params.csv` |

They refit both sites with `maxwell_alt.jl` and compare with the EV and EVE fits of the manuscript. `src/checks/` contains the numerical checks of the same formulation (Laplace inversion, time stepping, velocity comparison).

## Other numerical checks

`src/annex/` checks three properties of the half-Laplacian formulation numerically:

- `verify_fractional_kernel.jl`: the normalization constant of the fractional Laplacian at `s = 1/2`, and the agreement between its principal-value form and the Hadamard finite part of the Stokes kernel.
- `verify_hilbert_symbol.jl`: the derivative of the Hilbert transform has the Fourier symbol `|k|`, that of `(−Δ)^{1/2}` (checked with FFT).
- `verify_sobolev_divergence.jl`: the Heaviside initial condition has an infinite `H^{1/2}` (Gagliardo) seminorm.

## Reports (`reports/`)

Fit reports and best-fit parameters written by the scripts above (`*_params.csv`, `*_report.txt`, `*_chi2_landscape_*.txt`), and three files tracking the verification of the bibliography (`reference_*.md`).
