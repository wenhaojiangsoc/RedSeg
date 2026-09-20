# Appendix J: descriptive associations of exclusionary ZONING (Wharton WRLURI 2018) and
# GENTRIFICATION with fractional integration segregation. MSA + category FE, county-clustered;
# FE-only vs full contemporary controls. Associations only (no boundary variation exists).
suppressMessages({library(tidyverse); library(fixest)})
load("./data/safegraph_2019_m.rdata")
FRC<-as_tibble(readRDS("results/poi_ei_bounds_fraccorr.rds"))%>%transmute(safegraph_place_id,seg_frac=seg_eco*5/8)
d<-places_usa_2019%>%inner_join(FRC,by="safegraph_place_id")%>%
  mutate(county=substr(sprintf("%05d",as.integer(fips)),1,5),zWRLURI=as.numeric(scale(WRLURI18)))%>%
  filter(!is.na(cbsa),!is.na(top_category))
ctrl<-intersect(c("tot_pop","community","q_inc","pct_hs_above","pct_labor_force","age15_34",
  "foreign_born","public_transit","walkability","near_highway_1000"),names(d))
run<-function(v,controls){dd<-d%>%filter(!is.na(.data[[v]]),if_all(all_of(ctrl),~!is.na(.)))
  f<-if(controls) as.formula(paste0("seg_frac~",v,"+",paste(ctrl,collapse="+"),"|cbsa+top_category"))
     else as.formula(paste0("seg_frac~",v,"|cbsa+top_category"))
  m<-feols(f,dd,cluster=~county)
  for(nm in grep(v,names(coef(m)),value=TRUE))   # a factor (e.g. gentrify) prints one line per level
    cat(sprintf("%-28s ctrl=%d: b=%+.5f se=%.5f p=%.3g N=%s\n",nm,controls,coef(m)[nm],se(m)[nm],pvalue(m)[nm],format(nobs(m),big.mark=",")))}
for(v in c("zWRLURI","gentrify")) for(k in c(FALSE,TRUE)) run(v,k)
