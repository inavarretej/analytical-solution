# Reporte: Trade-off reológico — Asal (`fig_rheology_1sigma_asal.jl`)

## Propósito

Construcción del diagrama de trade-off entre viscosidad η y espesor elástico H a partir del parámetro ajustado a* de Asal. Genera `fig/fig_rheology_1sigma_asal.png`.

---

## Relación física

El parámetro de difusividad `a` es una combinación de η, H y E:

```
a = E·H / (2η(1−ν²))   →   H(η) = 2aη(1−ν²) / E
```

Cada valor de a define una línea recta en el plano log(η) − H. La incertidumbre en E [20, 80 GPa] define una banda.

---

## Parámetros usados

| Parámetro | Valor | Origen |
|-----------|-------|--------|
| a* | 9.503×10⁻⁶ m/s | `fig_chi2_landscape_asal.jl` |
| a (1σ inferior) | 9.13×10⁻⁶ m/s | contorno 1σ en a |
| a (1σ superior) | 9.89×10⁻⁶ m/s | contorno 1σ en a |
| ν | 0.25 | asumido |
| E_ref | 40 GPa | referencia |
| E_min / E_max | 20 / 80 GPa | rango de incertidumbre |

---

## Resultados (a H = 10 km)

| | η (Pa·s) |
|--|---------|
| **Best-fit** | **2.24×10¹⁹** |
| 1σ estadístico (en a) | [2.07, 2.36]×10¹⁹ |
| 1σ combinado (a + E) | [1.12, 4.45]×10¹⁹ |

**Conclusión**: la incertidumbre en E domina completamente sobre la estadística (σ_log10(E) > σ_log10(a)).

---

## Puntos de la literatura

| Referencia | η (Pa·s) | H (km) | Método |
|------------|---------|--------|--------|
| Nooner (2009) — Dabbahu, Afar | 5.2×10¹⁸ | 13.2 | relajación post-rifting |
| Wright et al. (2012) — Dabbahu, Afar | 1–3×10¹⁸ | 10 | relajación post-rifting viscoelástica |
| Cattin et al. (2005) — Asal | [0.5–1.0]×10¹⁸ | 3–12 | modelo FEM 2D termo-mecánico |

---

## Dependencias

- `fig_chi2_landscape_asal.jl` → provee a*, a_1σ_lo, a_1σ_hi (hardcodeados)
- No requiere datos externos
