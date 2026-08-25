# Reference Status Report — `paper_v2.tex`

Last updated: 2026-07-20. Consolidates and supersedes the counts in
`reference_coverage.md` (2026-06-08) and `reference_verification.md`
(2026-06-08) after the citation fixes made this session. `refs.bib` has
70 entries; 52 are now cited in `paper_v2.tex`.

---

## 1. Citations fixed this session

Six citations existed in `refs.bib` but were written as plain text in the
manuscript prose instead of `\citep`/`\citet`, so they would **not** have
appeared in the printed reference list and were not hyperlinked. All six
now use proper natbib commands:

| Line | Before | After | Bib key |
| --- | --- | --- | --- |
| 464 | `(B\"{u}rgmann and Dresen, 2008)` | `\citep{burgmann2008rheology}` | `burgmann2008rheology` |
| 580 | `(Laurie, 1997)` | `\citep{laurie1997calculation}` | `laurie1997calculation` |
| 864 | `(AIC; Akaike, 1974)` | `\citep[AIC;][]{akaike1974new}` | `akaike1974new` |
| 867 | `(Stone, 1974)` | `\citep{stone1974cross}` | `stone1974cross` |
| 874 | `(Wilks, 1938)` | `\citep{wilks1938large}` | `wilks1938large` |
| 1112 | `(Holtzman, 2016)` | `\citep{holtzman2016questions}` | `holtzman2016questions` |

`holtzman2016questions` was already correctly cited once elsewhere
(§2.1); this fixes a second, loose mention of the same source.

Recompiled clean: `natbib`/`bibtex` report zero undefined citations.

---

## 2. Still missing — no bib entry at all

Two in-text mentions are explicitly flagged `[CITE NEEDED]` and have no
corresponding entry in `refs.bib`. These need a full reference (or
removal) from you before submission — I can't fabricate one.

| Line | Text | Context | Prior status |
| --- | --- | --- | --- |
| 1120 | Hamling et al., 2014 | Dabbahu viscosity comparison, §4.2 | Already flagged in `reference_verification.md` (2026-06-08): "Could not locate PDF." Also absent from the Obsidian vault. |
| 1123 | Smittarello et al. (2016) | Asal magma-supply vs. viscoelastic relaxation, §4.2 | Not previously investigated in any report. |

---

## 3. Cited entries with a known content issue (already corrected)

These are not citation-*mechanics* problems (the `\cite` calls are fine)
but factual issues with the cited content, caught in earlier sessions:

| Bib key | Issue | Status |
| --- | --- | --- |
| `heki1993plate` | Wrong title in `refs.bib` (was "Mt. Spurr GRL") | Fixed 2026-05-27 |
| `nooner2009post` | Wrong title in `refs.bib` (was "Axial Seamount G3") | Fixed 2026-05-27 |
| `wright2012geophysical` | Elastic thickness plotted as H=13 km; source states ~10 km | Fixed this session in the new `fig_rheology_1sigma_krafla.jl` (now H=11 km, avoiding overlap with the Hofton & Foulger marker at H=10 km) |

One value used in the old (deleted) Krafla rheology figure, **Árnadóttir et al. (2009)** (η = 5–9×10¹⁸ Pa·s, H = 8 km), was never added to `refs.bib` and has no vault note — it was dropped rather than carried into the new `fig_rheology_1sigma_krafla.jl`. No action needed unless you want it reinstated with a real source.

---

## 4. Cited entries with a mismatched or missing reading note

The Obsidian vault (`10-literature/notes/`) has a note for a
*same-author, same-year* paper that is **not** the one actually cited in
`refs.bib` — easy to mistake as "already read":

| Bib key (cited) | Bib title | Vault note found | Vault note's actual title | 
| --- | --- | --- | --- |
| `hetland2006interseismic` | *Interseismic strain accumulation: Spin-up, cycle invariance...* | `hetlandEffectsRheologicalLayering2006` | *The effects of rheological layering on post-seismic deformation* | 

That vault note is a different Hetland & Hager (2006) paper. The one
actually cited in the CSS section (`\citep{thatcher1984viscoelastic,hetland2006interseismic}`, line 367) has no reading note yet.

Genuinely unread / no vault PDF (carried over from `project_paper_v2.md`'s
2026-05-27 audit, re-checked today — still true):

| Bib key | Reference |
| --- | --- |
| `barbot2010unified` | Barbot & Fialko (2010) |
| `hamling2009geodetic` | Hamling et al. (2009) — note: different paper from the missing Hamling (2014) in §2 above |
| `pollitz1992postseismic` | Pollitz (1992) — vault only has Pollitz (1997), a different paper, itself marked `status: unread` |
| `ruegg1984seismotectonic` | Ruegg et al. (1984) — vault only has Ruegg 1979 and 1987 |
| `sigmundsson2010geodetic` | Sigmundsson et al. (2010) |
| `cattin2005numerical` | Cattin et al. (2005) |
| `ledain1979seismotectonic` | Le Dain et al. (1979) |
| `reiss2022shallow` | Reiss et al. (2022) |
| `ali2014geodetic` | Ali et al. (2014) |
| `hetland2006thesis` | Hetland (2006) PhD thesis |
| `samko1993fractional`, `burnham2002model` | Reference books — low priority |

Resolved since the 2026-05-27 audit: `savage1983dislocation` and
`thatcher1984viscoelastic` both have vault notes now (created from
domain knowledge, no PDF in Zotero) — but note `savage1983dislocation`
is currently **uncited** in `paper_v2.tex` (see §5).

`akaike1974new`, `stone1974cross`, `wilks1938large`, `laurie1997calculation`,
`burgmann2008rheology` (the 5 newly-fixed citations) have no vault notes.
Four of these are standard statistics/numerical-methods references
(AIC, cross-validation, likelihood-ratio test, quadrature) — low priority,
same tier as `burnham2002model`. `burgmann2008rheology` is a genuine
geodynamics review and probably worth a real note eventually.

---

## 5. Uncited entries in `refs.bib` (18 / 70)

Entries present in the bibliography file but never referenced in
`paper_v2.tex`. Three belong to the separate fractional-Laplacian annex
(`tex/annex/annex_fractional.tex`, which has its own standalone
bibliography and is not `\input`-ed into `paper_v2.tex`) — not a gap,
just a different document:

| Bib key | Reference |
| --- | --- |
| `caffarelli2007extension` | Caffarelli & Silvestre (2007) |
| `daoud2022fractional` | Daoud & Laamri (2022) |
| `dinezza2012hitchhikers` | Di Nezza, Palatucci & Valdinoci (2012) |

The remaining 15 are candidate references prepared (mostly during the
2026-05-27 vault audit) but not yet woven into the current prose. Not
errors — just unused inventory. Flagging in case any belong in the
Discussion (e.g. the Afar/Iceland tectonic-context papers) or were meant
to replace a currently-cited item:

| Bib key | Reference |
| --- | --- |
| `abramowitz1964handbook` | Abramowitz & Stegun — special functions handbook |
| `ebinger1989tectonic` | Ebinger (1989) — western branch EAR tectonics |
| `freed2004evidence` | Freed & Bürgmann (2004) — power-law flow, Mojave |
| `hamling2009geodetic` | Hamling et al. (2009) — Dabbahu dyke intrusions |
| `hooper2011increased` | Hooper et al. (2011) — ice-cap retreat / magma capture, Iceland |
| `kenner2000postseismic` | Kenner & Segall (2000) — 1906 San Francisco EQ |
| `lee2016magmatic` | Lee et al. — CO₂ degassing, Ol Doinyo Lengai |
| `legall2008pliocene` | Le Gall et al. (2008) — Turkana depression rifting |
| `mcclusky2010kinematics` | McClusky et al. (2010) — Red Sea–Afar triple junction |
| `muirhead2016evolution` | Muirhead et al. (2016) — upper crustal faulting, magmatic volatiles |
| `pollitz1992postseismic` | Pollitz (1992) — postseismic relaxation, spherical Earth |
| `reiss2023crustal` | Reiss et al. (2023) — Natron rift crustal structure |
| `saemundsson1978fissure` | Sæmundsson (1978) — Iceland fissure swarms |
| `savage1983dislocation` | Savage (1983) — dislocation model, subduction zone |
| `sigmundsson2015segmented` | Sigmundsson et al. (2015) — Bárðarbunga segmented dyke growth |

---

## 6. Summary

| Category | Count |
| --- | --- |
| Total `refs.bib` entries | 70 |
| Cited in `paper_v2.tex` | 52 |
| Uncited — belongs to fractional-Laplacian annex | 3 |
| Uncited — unused candidate reference | 15 |
| Cited but plain-text (fixed this session) | 6 |
| Missing entirely, `[CITE NEEDED]` | 2 |
| Cited, factual issue (all now fixed) | 3 |
| Cited, vault note mismatched to a different paper | 1 (`hetland2006interseismic`) |
