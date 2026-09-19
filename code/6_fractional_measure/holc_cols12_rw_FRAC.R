# Columns 1 & 2 on the RECOMPUTED + REWEIGHTED segregation measure (seg_int_rw*5/8).
suppressMessages({library(tidyverse); library(haven); library(fixest)})
# run from the project root; # setwd("~/Library/CloudStorage/Dropbox/Demography-Rep-2026")
load("data/safegraph_2019_m.rdata"); rw<-readRDS("data/poi_seg_frac_as_rw.rds")
red<-read_csv("data/POI_redline_distance.csv",show_col_types=FALSE)%>%select(safegraph_place_id,holc_grade=any_of(c("holc_grade","grade","HOLCGrade")))
P<-places_usa_2019%>%left_join(rw%>%transmute(safegraph_place_id,seg_rw=seg_int_rw*5/8),by="safegraph_place_id")%>%
  left_join(red,by="safegraph_place_id")%>%
  mutate(redlining_CD=ifelse(!is.na(holc_grade)&holc_grade%in%c("C","D"),1L,ifelse(!is.na(holc_grade)&holc_grade%in%c("A","B"),0L,NA_integer_)),
         ltot=log(tot_pop+1),lvis=log(raw_visit_counts+1),county=substr(sprintf("%05d",as.integer(fips)),1,5))%>%
  filter(!is.na(seg_rw))
# ---- Col 1: national OLS, graded POIs only (B vs C/D), all available controls + MSA + category FE ----
ctrl<-intersect(c("ltot","community","q_inc","pct_hs_above","pct_labor_force","age15_34","foreign_born",
  "public_transit","walkability","suburban","near_highway_1000","gentrify","DRI18","AHI18"),names(P))
d1<-P%>%filter(!is.na(redlining_CD))
f1<-as.formula(paste0("seg_rw~redlining_CD+",paste(ctrl,collapse="+"),"|cbsa+top_category"))
m1<-feols(f1,d1,cluster=~county)
cat(sprintf("Col1 OLS (reweighted): est=%+.5f se=%.5f p=%.4g  N=%d  controls:%d\n",
  coef(m1)["redlining_CD"],se(m1)["redlining_CD"],pvalue(m1)["redlining_CD"],nobs(m1),length(ctrl)))
# ---- Col 2: naive border discontinuity on seg_rw ----
ahm_path<-"./data/ahm_staging/"; a2f<-function(x)paste0(substr(x,2,3),substr(x,5,7),substr(x,9,14)); BAR<-c("BUFFintrain","BUFFinsrivers","BUFFinbrivers")
poi<-P%>%mutate(tract11=substr(poi_cbg,1,11))%>%filter(!is.na(lvis))%>%select(tract11,sub_category,seg_rw,lvis)
rm2<-function(file,idcol,pfx){d<-read_dta(paste0(ahm_path,file))%>%filter(year==2010)%>%mutate(tract11=a2f(trctID2010))
  d%>%select(tract11,matches(paste0("^",idcol,"[0-9]+$")),matches("^boundary_side[0-9]+$"),matches("^holc_grade[0-9]+$"),
      matches("^boundary_grades[0-9]+$"),matches("^BUFFintrain[0-9]+$"),matches("^BUFFinsrivers[0-9]+$"),matches("^BUFFinbrivers[0-9]+$"))%>%
    rename_with(~sub(paste0("^",idcol),"boundary_id",.x))%>%
    pivot_longer(-tract11,names_to=c(".value","slot"),names_pattern="^(boundary_id|boundary_side|holc_grade|boundary_grades|BUFFintrain|BUFFinsrivers|BUFFinbrivers)([0-9]+)$")%>%
    drop_na(boundary_id)%>%mutate(boundary_id=paste0(pfx,boundary_id))}
tr<-rm2("redlining_dataset_t2010_4th_new.dta","boundary_id","A_")%>%filter(boundary_grades%in%c("2_to_3","2_to_4"))%>%
  mutate(ba=rowSums(across(all_of(BAR),~replace_na(.x,0))))
tr<-tr%>%semi_join(tr%>%group_by(boundary_id)%>%summarise(b=as.integer(max(ba,na.rm=T)>0.05))%>%filter(b==0),by="boundary_id")%>%
  mutate(lgs=as.integer(holc_grade%in%c("C","D")))%>%group_by(boundary_id)%>%filter(n_distinct(lgs)==2)%>%ungroup()
dat<-tr%>%distinct(boundary_id,boundary_side,tract11,lgs)%>%inner_join(poi,by="tract11",relationship="many-to-many")
m2<-feols(seg_rw~lgs|boundary_id+sub_category,dat,weights=~lvis,cluster=~boundary_id)
cat(sprintf("Col2 naive border (reweighted): est=%+.5f se=%.5f p=%.4g  boundaries=%d POIs=%d\n",
  coef(m2)["lgs"],se(m2)["lgs"],pvalue(m2)["lgs"],n_distinct(dat$boundary_id),nrow(dat)))
