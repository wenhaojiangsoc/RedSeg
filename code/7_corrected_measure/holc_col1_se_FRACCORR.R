suppressMessages({library(tidyverse); library(fixest); library(data.table)})
load("./data/safegraph_2019_m.rdata")
rl<-fread("./data/POI_redline_distance.csv")[,.(safegraph_place_id,holc_grade)]
FRC<-as_tibble(readRDS("results/poi_ei_bounds_fraccorr.rds"))%>%transmute(safegraph_place_id,seg_frac=seg_eco*5/8)
d<-places_usa_2019%>%left_join(rl,by="safegraph_place_id")%>%left_join(FRC,by="safegraph_place_id")%>%
  mutate(g=ifelse(is.na(holc_grade)|!(holc_grade%in%c("A","B","C","D")),"U",holc_grade),
         redCD=as.integer(g%in%c("C","D")),county=substr(sprintf("%05d",as.integer(fips)),1,5))%>%
  filter(!is.na(cbsa),!is.na(top_category),g%in%c("B","C","D","U"),!is.na(seg_frac))
ctrl<-intersect(c("tot_pop","community","q_inc","pct_hs_above","pct_labor_force","age15_34",
  "foreign_born","public_transit","walkability","near_highway_1000","gentrify"),names(d))
dd<-d%>%filter(if_all(all_of(ctrl),~!is.na(.)))
f0<-feols(seg_frac~redCD|cbsa+top_category,dd,cluster=~county)
f1<-feols(as.formula(paste0("seg_frac~redCD+",paste(ctrl,collapse="+"),"|cbsa+top_category")),dd,cluster=~county)
cat(sprintf("FE-only:    %+.5f (se %.5f, p=%.2g) N=%d\nControlled: %+.5f (se %.5f, p=%.2g) N=%d\n",
  coef(f0)["redCD"],se(f0)["redCD"],pvalue(f0)["redCD"],nobs(f0),
  coef(f1)["redCD"],se(f1)["redCD"],pvalue(f1)["redCD"],nobs(f1)))
