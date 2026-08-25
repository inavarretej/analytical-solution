# Reporte: Perfil de velocidad Asal (`fig_asal_velocity_profile.jl`)

## Propósito

Predicción hacia adelante del perfil de velocidad espacial de Asal usando los parámetros best-fit obtenidos de la serie temporal (`asal_fit.jl`). Genera `fig/fig_asal_velocity_profile.png`.

---

## Datos

- **Fuente**: `data/asal_profile.csv`, digitalizado de Vigny et al. (2007) Fig. 5
- **N**: 15 estaciones
- **Observable**: velocidad de apertura de rift, marco Somalia (mm/yr)
- **Época del perfil**: t = 19 yr post-rifting (media de campañas GPS ~1997, Vigny §3.2)
- **Centro del modelo**: x = +1 km NE del eje topográfico (Petit Rift, Vigny §3.2)

---

## Parámetros (hardcodeados de `reports/asal_report.txt`)

| Modelo | Parámetros |
|--------|-----------|
| Step viscoso | a = 9.503×10⁻⁶ m/s, b = 0.7103 m |
| Foulger (ajustado) | κ = 1.040×10⁻² m²/s, U₀ = 0.4381 m |
| Foulger (Vigny 2007) | κ = 0.015 m²/s, U₀ = 0.4 m |

---

## Velocidad de placa (V_plate)

La predicción transiente no incluye la velocidad de fondo Somalia–Danakil. Se ajusta un único parámetro V_plate por mínimos cuadrados ponderados (marginalización analítica):

```
V_plate = Σ[wᵢ(v_obs,i − v_model,i)] / Σwᵢ     wᵢ = 1/σᵢ²
```

| Modelo | V_plate (mm/yr) |
|--------|----------------|
| Step viscoso | ~8.1 |
| Foulger (ajustado) | ~7.9 |
| Foulger (Vigny ref) | ~8.0 |

Referencia Somalia–Danakil: ~16 mm/yr (Vigny et al. 2007). Los valores ajustados (~8 mm/yr) son menores porque el perfil GPS incluye estaciones próximas al rift que aún tienen señal transiente.

---

## Calidad del ajuste

| Modelo | χ²_r (con V_plate, k_eff=1) |
|--------|---------------------------|
| Step viscoso | ~18 |
| Foulger (ajustado) | ~33 |

χ²_r >> 1 en ambos modelos: el perfil espacial tiene estructura no modelada (movimiento del bloque Danakil, heterogeneidades). El Step viscoso es preferido sobre Foulger por un factor ~2 en χ²_r.

---

## Notas

- Este script es una **predicción forward**, no un ajuste independiente al perfil.
- Los parámetros transientes vienen exclusivamente de la serie temporal.
- Script actualizado respecto a `asal_profile_fit.jl` (May 25): usa `asal_timeseries_merged.csv` (N=33) en lugar de `asal_timeseries.csv` (N incompleto), y marginaliza V_plate analíticamente.

---

## Dependencias

- `src/models.jl` — `get_v_step`, `get_v_foulger`
- `data/asal_profile.csv`
- Parámetros de: `reports/asal_report.txt`
