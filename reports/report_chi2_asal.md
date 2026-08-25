# Reporte: Paisaje χ² — Asal (`fig_chi2_landscape_asal.jl`)

## Propósito

Mapeo del espacio de parámetros (a, b) para el modelo Step viscoso ajustado a Asal usando la **serie temporal** de desplazamiento. Genera `fig/fig_chi2_landscape_asal.png`.

---

## Datos

- **Fuente**: `data/asal_timeseries_merged.csv`, Vigny et al. (2007)
- **Observable**: ΔL(t) en mm (separación EP00–DF00), x₀ = 2000 m
- **N** = 33 observaciones (1979–2004)
- **t₀**: 1978.85 (onset del rifting)

---

## Método

- **Grid**: 500 × 400 puntos en (a, b)
  - a ∈ [8.0×10⁻⁶, 1.1×10⁻⁵] m/s
  - b ∈ [0.65, 0.78] m
- **χ²/N**: offset temporal (t_ref = primer dato, ~1.1 yr) restado del modelo — sin parámetro libre de offset
- **Predicción**: ΔL(t) = 2 × u(x₀, t; a, b) × 1000 [mm]

### Contornos de confianza

| Contorno | Δχ²/N | Nivel |
|----------|--------|-------|
| 1σ | 2.30/N | 68.3% |
| 2σ | 6.17/N | 95.4% |
| 3σ | 11.83/N | 99.7% |

---

## Resultados

| Parámetro | Valor | Límite 1σ inferior | Límite 1σ superior |
|-----------|-------|-------------------|-------------------|
| a* | 9.503×10⁻⁶ m/s | 9.13×10⁻⁶ m/s | 9.89×10⁻⁶ m/s |
| b* | 0.7103 m | — | — |
| χ²/N mínimo | ver figura | | |

### Implicaciones para η (H = 10 km, E = 40 GPa, ν = 0.25)

| | η (Pa·s) |
|--|---------|
| Best-fit | 2.24×10¹⁹ |
| 1σ inferior (a máx) | ~2.07×10¹⁹ |
| 1σ superior (a mín) | ~2.36×10¹⁹ |

Nota: la incertidumbre estadística en η es pequeña (~5%); la incertidumbre en E domina (ver `report_rheology_asal.md`).

---

## Dependencias

- `src/models.jl` — `get_u_step`
- `data/asal_timeseries_merged.csv`
- Resultados usados por: `fig_rheology_1sigma_asal.jl`
