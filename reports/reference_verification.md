# Reference Verification — Rheology Figures

Verification of literature points plotted in `fig_rheology_1sigma.png` (Krafla)
and `fig_rheology_1sigma_asal.png` (Asal). Values cross-checked against PDFs in Zotero.

---

## Krafla figure (`fig_rheology_1sigma.jl`)

| Reference | Plotted values | PDF source | Status | Notes |
|---|---|---|---|---|
| Sigmundsson (1991) | η < 10¹⁹ Pa·s, H = 10 km (upper bound) | Zotero: FMIJJDPU | ✅ VERIFIED | Paper: "1×10¹⁹ Pa s is a maximum value for asthenosphere viscosity in Iceland." H=10 km used as approximate elastic layer thickness (maximum hypocentral depth). |
| Hofton & Foulger (1996) | η = 1.1×10¹⁸ Pa·s, H = 10 km | Zotero: D7TLWXZ4 | ✅ VERIFIED | Abstract explicitly: "half-space viscosity of 1.1 × 10¹⁸ Pa s…elastic layer thickness for northeast Iceland of 10 km." Table 3 confirms both values. |
| Wright et al. (2012) | η = 1–3×10¹⁸ Pa·s, H = 13 km | Zotero: FUC9GGQT | ⚠️ PARTIAL — H IS WRONG | η range confirmed. **H should be ~10 km, not 13 km** — paper states "~10 km thick" (p. 244). Also: these values cite refs 70–71 (Hofton & Foulger 1996; Heki et al. 1993), not original to Wright. → **Fix: move bar to H = 10–11 km.** |
| Árnadóttir et al. (2009) | η = 5–9×10¹⁸ Pa·s, H = 8 km | Not in Zotero | ❌ NOT VERIFIED | Could not locate PDF. Values cannot be confirmed. |

## Asal figure (`fig_rheology_1sigma_asal.jl`)

| Reference | Plotted values | PDF source | Status | Notes |
|---|---|---|---|---|
| Cattin et al. (2005) | η = 0.5–1×10¹⁸ Pa·s, H = 3–12 km | Zotero: XM8HB84V | ✅ VERIFIED | Post-1978 Asal-Ghoubbet relaxation. Near-axis effective viscosity from 2D thermo-mechanical FEM. Seismogenic thickness 3 km (beneath Fieale) to ~12 km off-axis. |
| Nooner et al. (2009) | η = 5.2×10¹⁸ Pa·s, H = 13.2 km | Zotero: HKMEXL8K | ✅ VERIFIED | "The best fit model has an elastic thickness of 13.2 km and a viscosity of 5.2 × 10¹⁸ Pa·s." **Region: Dabbahu rift, Afar** — correctly labeled in figure. |
| Hamling et al. (2014) | η ≈ 10¹⁸ Pa·s, H = 21 km | Not in Zotero | ❌ NOT VERIFIED | Could not locate PDF. Values cannot be confirmed. Region: Dabbahu, Afar. |

---

## Action items

1. **Fix Wright et al. (2012)** in `fig_rheology_1sigma.jl`: change H from 13 → 10 km (or 11 km to avoid visual overlap with Hofton & Foulger at H=10).
2. **Árnadóttir et al. (2009)** — decide: remove from Krafla figure until PDF is located, or keep with explicit "unverified" caveat.
3. **Hamling et al. (2014)** — decide: remove from Asal figure until PDF is located, or keep with caveat.

---

*Verified: 2026-06-08*
