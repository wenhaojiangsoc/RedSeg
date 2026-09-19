# Publication numbers for the BROAD descriptive comparison C/D vs (B + ungraded):
# FE-only and full-contemporary-controls, with COUNTY-clustered SEs (as Table 1, Col 1).
suppressMessages({library(tidyverse); library(fixest); library(data.table)})
load("./data/safegraph_2019_m.rdata")
cat("county-ish cols:", paste(grep("county|countycode|fips|cbsa", names(places_usa_2019), value=TRUE, ignore.case=TRUE), collapse=", "), "\n")
rl <- fread("./data/POI_redline_distance.csv")[,.(safegraph_place_id, holc_grade)]
d <- places_usa_2019 %>% left_join(rl, by="safegraph_place_id") %>%
  mutate(seg_int = race_seg_msa_adj_mean,
         g = ifelse(is.na(holc_grade) | !(holc_grade %in% c("A","B","C","D")), "U", holc_grade)) %>%
  filter(!is.na(cbsa), !is.na(top_category), !is.na(seg_int), g %in% c("B","C","D","U")) %>%
  mutate(redCD = as.integer(g %in% c("C","D")))
cl <- if ("countycode" %in% names(d)) "countycode" else if ("county" %in% names(d)) "county" else grep("county|fips", names(d), value=TRUE, ignore.case=TRUE)[1]
cat("clustering by:", cl, "\n")
soc <- intersect(c("tot_pop","community","q_inc","pct_hs_above","pct_labor_force","age15_34",
         "foreign_born","public_transit","walkability","near_highway_1000","gentrify"), names(d))
sh <- function(m,t) cat(sprintf("  %-22s redCD=%+.5f (clustered se %.5f, p=%.3g)  n=%d\n",
                                t, coef(m)["redCD"], se(m)["redCD"], pvalue(m)["redCD"], m$nobs))
form <- function(rhs) as.formula(paste0("seg_int~redCD",rhs,"|cbsa+top_category"))
sh(feols(form(""), d, cluster=as.formula(paste0("~",cl))), "FE only")
sh(feols(form(paste0("+",paste(soc,collapse="+"))), d, cluster=as.formula(paste0("~",cl))), "+ full controls")
sh(feols(form(paste0("+",paste(c(soc,"black_rate","white_rate"),collapse="+"))), d, cluster=as.formula(paste0("~",cl))), "+ controls & race comp")
