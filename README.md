# Replication package for The Color of Activity Space
Wenhao Jiang and Yongjun Zhang. Manuscript under review at *Demography* (revise & resubmit).

This repository contains all analysis code needed to reproduce the tables and figures of the
paper and its appendix. It contains **code only**: the mobility data are proprietary and the
historical data are redistributed by their original authors, so `data/` documents how to obtain
and stage every input rather than shipping it (see [Data requirements](#data-requirements)).

---

## What the paper does

The paper measures **integration segregation** for ~4 million U.S. establishments (POIs): the
distance between a location's 2019 visitor racial composition (SafeGraph mobility flows, visitor
race inferred from home census block groups) and its metropolitan racial composition. It then
estimates the causal effect of 1930s HOLC "redlining" on this outcome with a **boundary
difference-in-discontinuities** design after Aaronson, Hartley, and Mazumder (2021, *AEJ:
Policy*; "AHM"): the segregation discontinuity at real B–C and B–D grade boundaries is
benchmarked against placebo discontinuities at same-grade counterfactual grid boundaries,
reweighted by **entropy balancing** (Hainmueller 2012) so the placebo pool reproduces the real
boundaries' pre-map (1910–1930) covariate gaps exactly. Placebo boundaries lack a "worse" side,
so each is oriented at random and every estimate is averaged over 100 placebo-orientation draws
(the mediation and EI analyses use draws 1/3/13; each script states its draws).

## Repository layout

```
RedSeg/
├── README.md                  <- this file
├── data/README.md             <- how to obtain and stage every input dataset
└── code/
    ├── 0_data_construction/   <- raw flows -> analysis files (run once, in order)
    ├── 1_main_results/        <- Table 1 (descriptive gradient + boundary DiD, 5 columns)
    ├── 2_mechanism/           <- mediation: composition, groups, amenity mix, income, Place Pulse
    ├── 3_heterogeneity/       <- establishment types and catchment bins
    ├── 4_appendix/            <- Appendices A-I (robustness, EI, Oster, diversity, balanced borders)
    ├── 5_figures/             <- Figures 1-5 (each reads results/ from the stages above)
    └── 6_fractional_measure/  <- the paper's FINAL estimates: every analysis re-run with
                                  fractional (proportional) visitor-race attribution (see below)
```

> **Which numbers are in the paper?** The published estimates use **fractional attribution**:
> each visit contributes its home block group's full ACS racial distribution, exactly as the
> paper's Eq. 1 defines. Stages 0-5 build the pipeline and the largest-group (one-hot) variant
> retained as a robustness row in Appendix B; `6_fractional_measure/` re-runs every estimator
> on the fractional outcome (`ei_bounds_seg_frac.R` output, plus a fractional re-stream of the
> excluding-local-residents and diversity measures in `holc_nonlocal_seg_FRAC.R`) and produces
> the numbers, tables, and figures that appear in the manuscript. The `*_FRAC.R` scripts are
> one-line path variants of the stage-1-to-4 originals (via the drop-in file
> `data/poi_seg_frac_as_rw.rds`); `holc_frac_core.R`, `holc_frac_batch.R`, and
> `holc_frac_figprep.R` carry the headline DiDs, the mechanism/heterogeneity batch, and the
> figure regeneration.

All R scripts assume the **working directory is the project root** (the directory holding
`data/` and `results/`; scripts create `results/` outputs by relative path). The five figure
scripts write the final PDFs by an absolute path constant (`OUT` or the `ggsave()` call at the
bottom) pointing at the paper's LaTeX directory — **edit that one line** to your own output
location before running them.

## Data requirements

See `data/README.md` for full acquisition instructions and the exact directory tree the scripts
expect. In brief:

| Input | Source | Access |
|---|---|---|
| SafeGraph Monthly Patterns, 2019 backfill (visits + `visitor_home_cbgs` by POI-month) | SafeGraph / Dewey Data | **Proprietary license required** |
| `safegraph_2019_m.rdata` (assembled place-level file: POI attributes, published segregation measures, MSA, CBG race) | built from SafeGraph + ACS in the authors' prior project | derived; structure documented in `data/README.md` |
| AHM (2021) replication data (`redlining_dataset_t2010_4th_new.dta`, `redlining_dataset_tgrid2010_4th_new.dta`, ...) | Aaronson–Hartley–Mazumder replication archive (openICPSR) | public |
| Mapping Inequality HOLC polygons (`mappinginequality.json`) | Digital Scholarship Lab, U. Richmond | public |
| ACS 2019 5-year, CBG race (B03002) and population | Census API via `tidycensus` | public (API key) |
| Place Pulse 2.0 street-image perception ratings | MIT Media Lab | public |
| Wharton land-use index, gentrification indicators, POI–HOLC distance (`POI_redline_distance.csv`) | assembled by the authors (WRLURI; Hwang–Ding-style gentrifiability) | derived |

## Software

R ≥ 4.2 with: `tidyverse`, `data.table`, `haven`, `fixest`, `sf`, `maptiles`, `tidyterra`,
`patchwork`, `tidycensus`. Python ≥ 3.9 (standard library only) for the flow extractor.
No seeds beyond those set explicitly in the scripts; entropy balancing is deterministic.
The heavy steps (`compute_seg_reweighted.R`, the `ei_*_seg` scripts, `holc_nonlocal_seg.R`)
stream the ~12 monthly flow files and need ~32 GB RAM; estimator scripts run in minutes.

## Execution order and script map

### Stage 0 — data construction (`code/0_data_construction/`, run in this order)

| # | Script | What it does | Output |
|---|---|---|---|
| 0.1 | `extract_flows_2019.py` | Streams raw SafeGraph monthly patterns; writes one `POI x home-CBG x visitors` file per month (US CBGs only). Run once per month directory. | `data/flows_2019/*.csv.gz` |
| 0.2 | `build_cbg_race_2019.R` | CBG-level race shares, population, and SafeGraph device counts (for post-stratification). | `data/cbg_race_2019.csv` |
| 0.3 | `pull_cbg_race_frac.R` | Fractional CBG racial composition (ACS B03002) for the ecological-inference appendix. | `data/cbg_race_frac.csv` |
| 0.4 | `compute_seg_reweighted.R` | Recomputes integration segregation from the raw flows, unweighted (`seg_int_unw`) and with SafeGraph's recommended CBG post-stratification reweighting (`seg_int_rw`). This is the largest-group outcome variant (Appendix B robustness row; Appendix I evaluates the reweighting). | `data/poi_seg_reweighted.rds` |
| 0.5 | `holc_poi_assign.R` | Point-level assignment of POIs to real B–(C/D) polygon borders within 1/4 mile (Column 5 and the design figure); also saves the full B/C/D point set for the placebo grid. | `data/poi_real_assign.rds`, `data/poi_bcd_pts.rds` |
| 0.6 | `placepulse_scores.R` | Place Pulse 2.0 composite perceived-quality scores by location (Appendix H). | `data/placepulse_location_scores_us.csv` |

### Stage 1 — main results (`code/1_main_results/`) — Table 1

| Script | Paper exhibit | Notes |
|---|---|---|
| `holc_control_check.R` | Col. 1 (0.049 vs 0.038 under largest-group attribution; the published 0.020 vs 0.017 come from `6_fractional_measure/holc_col1_se_FRAC.R`) | Descriptive C/D gradient, FE-only vs full contemporary controls (income, education, majority race, walkability, Wharton land-use regulation, gentrification), county-clustered. Also the source of the Appendix J zoning/gentrification associations. |
| `holc_cols12_rw.R` | Cols. 1–2 | Descriptive gradient and raw boundary discontinuity on the recomputed, reweighted measure. |
| `holc_cdb_ipw_rw.R` | Col. 3 | Probit inverse-propensity weighting alternative; 100 draws. |
| `holc_cdb_ebal_rw.R` | Col. 4 design (+0.0050 largest-group; the published **+0.0024** comes from `6_fractional_measure/holc_frac_core.R`) | Entropy-balanced boundary DiD, tract-buffer assignment; 100 draws. |
| `holc_cdb_poi.R`, `holc_col5_poi_pub.R` | Col. 5 design (+0.0079 largest-group; the published **+0.0035** comes from `6_fractional_measure/holc_frac_core.R`) | POI-level assignment (each POI by its own distance to the border): `holc_cdb_poi.R` on the recomputed measures, `holc_col5_poi_pub.R` on the published measure used in the table. |

### Stage 2 — mechanism (`code/2_mechanism/`) — Section "Mechanism", Figure 4

| Script | Paper exhibit |
|---|---|
| `holc_cdb_mediation_rw.R` | Compositional mediator: neighborhood distance from MSA composition (a = 0.012, b = 0.195, ~43% mediated). |
| `holc_cdb_mediation_groups.R` | Racial-group decomposition of the channel (Black–white axis); writes `results/holc_mediation_groups.rds`, the input to Figure 4. |
| `holc_cdb_amenity_med.R` | Amenity-mix rival channel: decomposes segregation into category mean + within-category remainder; the mix contributes 0–2% (precise null). |
| `holc_income_intersect.R` | Income-segregation placebo outcome (+0.0005, p = 0.70) and income-control robustness (+0.0044). |
| `holc_income_interaction.R` | Effect does not vary by neighborhood income (interaction p > 0.4). |
| `holc_cdb_med_placepulse_rw.R` | Perceived-quality (Place Pulse) mediation, composite + six dimensions: 0–2% (Appendix H). |

### Stage 3 — heterogeneity (`code/3_heterogeneity/`) — Section "Where the Effect Concentrates", Figure 5

| Script | Paper exhibit |
|---|---|
| `holc_category7.R` | Within-establishment-type DiD (religious +0.009; grocery/pharmacy/retail +0.008; full mapping = Appendix G, Table G1). Writes `results/holc_category7.rds` (Figure 5A input). |
| `holc_catch_multi.R` | DiD within six catchment-size bins (mean same-MSA visitor travel distance); most-local +0.008 to widest ≈ 0. Writes `results/holc_catch_multi.rds` (Figure 5B input). |

### Stage 4 — appendix (`code/4_appendix/`)

| Script | Appendix |
|---|---|
| `holc_appendix_robust.R` | A (Table A1): key-4 covariate balancing, 1/8-mile buffer, daily-routine vs other categories. |
| `holc_nonlocal_seg.R` → `holc_nonlocal_did.R` | A: recompute segregation excluding the POI's own CBG/tract residents; DiD survives. |
| `ei_bounds_seg_frac.R` → `ei_bounds_did_frac.R` | B (Table B1, row 2): fractional CBG racial attribution. |
| `ei_goodman_seg.R` → `ei_goodman_did.R` | B (Table B1, rows 3–4): Goodman (1953) ecological regression per establishment. |
| `ei_bounds_seg.R` → `ei_bounds_did.R` | B: homogeneity-restriction (tau) bounds variant discussed in the text. |
| `holc_oster_broad.R` | D (Table D1): Oster delta* on the descriptive gradients (integration 1.67; diversity 0.27). |
| `holc_diversity_did.R` | E: boundary DiD on the diversity benchmark (precise null). |
| `holc_lowprop_poi.R` | F (Table F2): balanced-boundary second identification (+0.009, p = 0.04; less-balanced +0.020). |
| `holc_lowprop_pretrend.R` | F (Table F1): flat pre-map trajectories on balanced boundaries. |
| `holc_pretrend.R` | Methods claim: reweighting balances pre-map gaps in every census (1910/1920/1930). |
| `holc_nonlocal_hetero.R` | A (Table A2): heterogeneity with local residents excluded, 12 types + 6 catchment bins; also validates the 12-type mapping against Table G1. **Reconstructed** (see below). |
| `holc_bench.R` | C (Table C1): boundary DiD on 2010 income, rent, and homeownership vs the segregation effect, in % of within-boundary SD. **Reconstructed** (see below). |

### Stage 5 — figures (`code/5_figures/`)

| Script | Figure | Reads |
|---|---|---|
| `holc_measure_fig.R` | Fig. 1 (measure construction) | self-contained schematic |
| `holc_design_map2.R` | Fig. 2 (design map, Los Angeles) | `mappinginequality.json`, `data/poi_bcd_pts.rds` |
| `holc_balance_fig.R` | Fig. 3 (covariate balance) | AHM staging data |
| `holc_mediation_groups_fig.R` | Fig. 4 (mechanism) | `results/holc_mediation_groups.rds` |
| `holc_hetero_fig.R` | Fig. 5 (heterogeneity) | `results/holc_category7.rds`, `results/holc_catch_multi.rds` |

## Conventions worth knowing

- **Outcome scaling.** Every segregation index is rescaled by `* 5/8` on load so the summed
  five-group absolute deviation lies in [0, 1], matching the paper's descriptive baseline.
- **Estimator.** `feols(seg ~ lgs + lgs:real | boundary_id + sub_category)`, weights =
  (balancing weight) x log(visits + 1), SEs clustered by boundary. `coef(lgs)` is the placebo
  discontinuity; `coef(lgs:real)` is the reported difference-in-discontinuities.
- **Barriers.** Boundaries split by a railroad or river (`BUFFin*` flags) are excluded
  throughout, as stated in the paper.
- **Tables** are typeset by hand in the paper's LaTeX from each script's console output and
  saved `results/*.rds`; scripts print the exact numbers they contribute.

## Reconstructed scripts

The code that originally produced **Table A2** (heterogeneity excluding local residents) and
**Table C1** (economic benchmarking) was not preserved; `holc_nonlocal_hetero.R` and
`holc_bench.R` reconstruct both from the package's own machinery and were verified against the
published tables. **Table A2** reproduces exactly (to the printed fourth decimal) for all six
catchment bins, including their km labels, and for eight of twelve establishment types; the
remaining four types agree within 0.0006 (the original assignment of a few marginal SafeGraph
categories to types is not recoverable, so the mapping follows Table G1's printed
descriptions). **Table C1** reproduces the segregation row exactly (+0.0050 over 100 draws)
and the income row at published rounding (−$1,070 → −$1,100, n.s.); homeownership matches to
0.1 pp (−1.29 vs −1.2, p = .015 vs p < .01), while the rent estimate is smaller here (−0.8%
vs −1.0%) and short of marginal significance — the original's exact rent transform and SD
weighting could not be recovered. Each script's header records this verification.

## Contact

Wenhao Jiang (wj93@duke.edu).
