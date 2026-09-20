# Oster (2019) delta* on the paper's ACTUAL descriptive spec: C/D vs (B + ungraded), full controls.
# Reports uncontrolled/controlled beta and R2, R_max=1.3*Rtilde, and delta* at R_max=1.3Rtilde and R_max=1.
# Integration (race_seg_msa_adj_mean, unscaled to match Table 1) and diversity (race_seg_unadj_mean).
suppressMessages({library(tidyverse); library(fixest); library(data.table)})
load("./data/safegraph_2019_m.rdata")
rl <- fread("./data/POI_redline_distance.csv")[,.(safegraph_place_id, holc_grade)]
FRC<-dplyr::as_tibble(readRDS("results/poi_ei_bounds_fraccorr.rds"))%>%dplyr::transmute(safegraph_place_id,seg_frac=seg_eco*5/8)
d <- places_usa_2019 %>% left_join(rl, by="safegraph_place_id") %>% left_join(FRC, by="safegraph_place_id") %>%
  mutate(g = ifelse(is.na(holc_grade) | !(holc_grade %in% c("A","B","C","D")), "U", holc_grade),
         redCD = as.integer(g %in% c("C","D"))) %>%
  filter(!is.na(cbsa), !is.na(top_category), g %in% c("B","C","D","U"))
ctrl <- intersect(c("tot_pop","community","q_inc","pct_hs_above","pct_labor_force","age15_34",
         "foreign_born","public_transit","walkability","near_highway_1000","gentrify"), names(d))
oster <- function(dv){
  dd <- d %>% filter(if_all(all_of(c(dv,ctrl)), ~!is.na(.)))
  f0 <- feols(as.formula(paste0(dv,"~redCD|cbsa+top_category")), dd)
  f1 <- feols(as.formula(paste0(dv,"~redCD+",paste(ctrl,collapse="+"),"|cbsa+top_category")), dd)
  b0<-coef(f0)["redCD"]; r0<-r2(f0,"r2"); b1<-coef(f1)["redCD"]; r1<-r2(f1,"r2")
  ds <- function(Rmax) unname((b1*(r1-r0))/((b0-b1)*(Rmax-r1)))
  cat(sprintf("%-22s beta0=%+.4f (R2=%.3f) | beta1=%+.4f (R2=%.3f) | Rmax(1.3R~)=%.3f | delta*=%.2f (Rmax=1.3R~), %.2f (Rmax=1) | n=%d\n",
    dv, b0,r0,b1,r1, min(1,1.3*r1), ds(min(1,1.3*r1)), ds(1), nrow(dd)))
}
oster("seg_frac")               # integration, FRACTIONAL attribution
oster("race_seg_msa_adj_mean")   # integration, majority (for reference)
oster("race_seg_unadj_mean")     # diversity
