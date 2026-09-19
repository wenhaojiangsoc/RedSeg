# Data staging

The scripts expect the following tree under the project root (this directory). Nothing in it is
committed to the repository: the mobility data are proprietary, and the public inputs are
redistributed by their original archives, which should be cited directly.

```
data/
├── flows_2019/                          <- output of extract_flows_2019.py, one file per month
│   └── 2019-01.csv.gz ... 2019-12.csv.gz
├── safegraph_2019_m.rdata               <- assembled place-level file (see below)
├── poi_seg_reweighted.rds               <- built by compute_seg_reweighted.R
├── poi_real_assign.rds                  <- built by holc_poi_assign.R
├── poi_bcd_pts.rds                      <- built by holc_poi_assign.R
├── cbg_race_2019.csv                    <- built by build_cbg_race_2019.R
├── cbg_race_frac.csv                    <- built by pull_cbg_race_frac.R
├── placepulse_location_scores_us.csv    <- built by placepulse_scores.R
├── POI_redline_distance.csv             <- POI-to-HOLC-boundary distances and grades (authors)
├── mappinginequality.json               <- Mapping Inequality HOLC polygons (see below)
└── ahm_staging/                         <- AHM (2021) replication .dta files (see below)
    ├── redlining_dataset_t2010_4th_new.dta
    ├── redlining_dataset_tgrid2010_4th_new.dta
    ├── redlining_dataset_t2010_8th_new.dta
    ├── redlining_dataset_tgrid2010_8th_new.dta
    └── HOLC_USA_4th_remap.dta
```

## Sources

**SafeGraph Monthly Patterns, 2019 backfill** (proprietary). Licensed via SafeGraph or Dewey
Data. Required fields: `safegraph_place_id`, `visitor_home_cbgs` (JSON of home-CBG visitor
counts), `raw_visit_counts`, place attributes (`sub_category`, latitude/longitude, `poi_cbg`).
`extract_flows_2019.py` converts each monthly directory into a flat `POI x CBG x visitors` file.

**`safegraph_2019_m.rdata`** is the authors' assembled 2019 place-level file
(`places_usa_2019`): one row per POI with identifiers, coordinates, `sub_category`, visit
counts, CBG/tract/MSA identifiers, CBG racial composition (`black_rate`, `tot_pop`), the
published segregation measures (`race_seg_msa_adj_mean`, `race_seg_unadj_mean`), income
segregation, catchment (mean same-MSA visitor travel distance), and the contemporary control
vector (neighborhood income, education, majority race, walkability, Wharton land-use index,
gentrification indicators). Researchers with SafeGraph access can rebuild it from Patterns +
ACS; contact the authors for the exact assembly code from the prior project.

**AHM replication data.** Aaronson, Hartley, and Mazumder (2021), "The Effects of the 1930s
HOLC 'Redlining' Maps" (AEJ: Economic Policy) — replication archive on openICPSR. Copy the
listed `.dta` files into `data/ahm_staging/`. These carry the tract-boundary assignments
(`boundary_id`, sides, grades, river/rail flags) and the 1910–1930 historical covariates.

**Mapping Inequality.** Nelson, Winling, et al., *Mapping Inequality: Redlining in New Deal
America*, Digital Scholarship Lab, University of Richmond. Download the full GeoJSON as
`mappinginequality.json`.

**ACS.** 2019 5-year tables via the Census API (`tidycensus`; set `CENSUS_API_KEY`): B03002
for fractional CBG race, CBG population for post-stratification.

**Place Pulse 2.0.** MIT Media Lab street-image perception collection (safety, wealth, beauty,
liveliness, depressing, boring). `placepulse_scores.R` converts pairwise votes to location
scores for U.S. cities.
