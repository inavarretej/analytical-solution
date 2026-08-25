# Reporte: `krafla_fit.jl`

## Propósito

Ajuste del perfil espacial de desplazamiento de Krafla (Islandia) usando datos de Foulger et al. (1992). El observable es la diferencia de desplazamiento entre dos épocas: τ₁ = 8 yr y τ₂ = 11 yr post-rifting (campañas 1984 y 1987). Genera `fig/fig08_fits_2p_3p.png`.

---

## Datos

- **Fuente**: Foulger et al. (1992), perfil de desplazamiento, N = 20 estaciones
- **Observable**: Δu(x) = u(x, τ₂) − u(x, τ₁) [m → cm para la figura]
- **Épocas**: τ₁ = 8 yr, τ₂ = 11 yr post-rifting

---

## Modelos y parámetros best-fit

Parámetros cargados de `reports/relaxation_report.txt` (generado 2026-04-21):

| Modelo | Parámetros | χ²/ν | AIC | BIC |
|--------|-----------|------|-----|-----|
| **Foulger 2p** | κ = 6.19×10⁻¹ m²/s, b = 0.810 m | 0.902 | — | — |
| **Step viscoso 2p** | a = 5.13×10⁻⁵ m/s, b = 1.144 m | 0.718 | — | — |
| **Maxwell 3p** | a = 6.61×10⁻⁵ m/s, b = 1.150 m, τ_R = 1.48 yr | 0.775 | — | — |

### Comparación estadística (N = 20)

| Modelo | χ² | χ²_r | AIC | BIC | MSPE |
|--------|----|------|-----|-----|------|
| Foulger 2p | — | 0.902 | — | — | — |
| Step 2p | — | 0.718 | — | — | — |
| Maxwell 3p | — | 0.775 | — | — | — |

**Test razón de verosimilitud (2p vs 3p)**: Δχ² = 2.53 < 3.84 (umbral 95%, 1 DOF) → τ_R **no es significativo**, el modelo Maxwell 3p no mejora el Step 2p.

---

## Implicaciones físicas

Usando a* del modelo Step viscoso (a = 5.13×10⁻⁵ m/s):

| Parámetro | Valor | Supuestos |
|-----------|-------|-----------|
| η (viscosidad) | ~4.2×10¹⁸ Pa·s | H = 10 km, E = 40 GPa, ν = 0.25 |
| Incertidumbre 1σ en η | factor ~4 | dominada por incertidumbre en E |

---

## Dependencias

- `src/models.jl` — `diff_step`, `diff_foulger`, `diff_viscoelastic`
- `data/extracted_data.csv` — perfil de Foulger et al. (1992)
- `reports/relaxation_report.txt` — parámetros best-fit cargados directamente
