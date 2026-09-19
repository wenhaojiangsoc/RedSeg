# REWEIGHTED-measure version of holc_cdb_ipw.R (outcome = seg_int_rw*5/8).
# CD-B boundary DiD with PROBIT-IPW weighting (robustness alternative to entropy balancing
# in holc_cdb_ebal.R). Same design: real 2->3/2->4 HOLC boundaries vs placebo grid boundaries,
# random low-grade-side orientation (RSUE) over 100 draws, visit-weighted feols.
# Weights: fit probit propensity of real-vs-placebo on the z-scored covariate GAPS (+ missingness
#   indicators), ATT odds weight placebo boundaries w = phat/(1-phat) [den. capped at .02, as in
#   holc_design_map.R], normalized to mean 1 within placebos so the scale matches the EB script
#   (real boundaries keep weight 1). All 1910/1920/1930 covariates enter the propensity model.
suppressMessages({library(tidyverse); library(haven); library(fixest)})
ahm_path<-"./data/ahm_staging/"; a2f<-function(x)paste0(substr(x,2,3),substr(x,5,7),substr(x,9,14))
BAR<-c("BUFFintrain","BUFFinsrivers","BUFFinbrivers")
covset<-c("black_","ownhome_","foreign_born_","read_write_","house_value_","rent_","radio_"); yrs<-c(1910,1920,1930)
load("./data/safegraph_2019_m.rdata")
rw<-readRDS("data/poi_seg_reweighted.rds")   # safegraph_place_id, seg_int_unw, seg_int_rw (raw index in [0,8/5])
poi<-places_usa_2019%>%left_join(rw,by="safegraph_place_id")%>%
  mutate(seg_unw=seg_int_unw*5/8,seg_rw=seg_int_rw*5/8,lvis=log(raw_visit_counts+1),tract11=substr(poi_cbg,1,11))%>%
  filter(!is.na(seg_rw),!is.na(seg_unw),!is.na(lvis))%>%select(tract11,sub_category,seg_unw,seg_rw,lvis)
read_map<-function(file,idcol,pfx,is_grid){
  d<-read_dta(paste0(ahm_path,file))%>%filter(year==2010)%>%mutate(tract11=a2f(trctID2010))
  if(is_grid) d<-d%>%rename_with(~sub("BUFF_CFin","BUFFin",.x))
  d%>%select(tract11,matches(paste0("^",idcol,"[0-9]+$")),matches("^boundary_side[0-9]+$|^side_grid[0-9]+$"),
      matches("^holc_grade[0-9]+$"),matches("^boundary_grades[0-9]+$"),
      matches("^BUFFintrain[0-9]+$"),matches("^BUFFinsrivers[0-9]+$"),matches("^BUFFinbrivers[0-9]+$"))%>%
    rename_with(~sub("^side_grid","boundary_side",.x))%>%rename_with(~sub(paste0("^",idcol),"boundary_id",.x))%>%
    pivot_longer(-tract11,names_to=c(".value","slot"),
      names_pattern="^(boundary_id|boundary_side|holc_grade|boundary_grades|BUFFintrain|BUFFinsrivers|BUFFinbrivers)([0-9]+)$")%>%
    drop_na(boundary_id)%>%mutate(boundary_id=paste0(pfx,boundary_id))
}
mA<-read_map("redlining_dataset_t2010_4th_new.dta","boundary_id","A_",FALSE)
mG<-read_map("redlining_dataset_tgrid2010_4th_new.dta","boundary_id_grid","G_",TRUE)
raw<-read_dta(paste0(ahm_path,"redlining_dataset_t2010_4th_new.dta"))%>%filter(year%in%yrs)%>%mutate(tract11=a2f(trctID2010))
cov_all<-raw%>%select(tract11,year,any_of(covset))%>%group_by(tract11,year)%>%summarise(across(everything(),~mean(.x,na.rm=T)),.groups="drop")%>%
  pivot_wider(names_from=year,values_from=any_of(covset),names_glue="{.value}{year}"); cvars<-setdiff(names(cov_all),"tract11")
tr<-mA%>%filter(boundary_grades%in%c("2_to_3","2_to_4"))%>%mutate(real=1L)
pl<-mG%>%filter(boundary_grades%in%c("2_to_2","3_to_3","4_to_4"))%>%mutate(real=0L)
bt<-bind_rows(tr,pl)%>%left_join(cov_all,by="tract11")%>%mutate(barrier_amt=rowSums(across(all_of(BAR),~replace_na(.x,0))))
bt<-bt%>%semi_join(bt%>%group_by(boundary_id)%>%summarise(b=as.integer(max(barrier_amt,na.rm=T)>0.05))%>%filter(b==0),by="boundary_id")

estimate<-function(seed,yv){
  side<-bt%>%group_by(boundary_id,real,boundary_side)%>%summarise(grade=first(holc_grade[holc_grade%in%c("A","B","C","D")]),
    across(all_of(cvars),~mean(.x,na.rm=T)),.groups="drop")%>%group_by(boundary_id)%>%filter(n()==2)%>%ungroup()
  set.seed(seed)
  side<-side%>%group_by(boundary_id)%>%mutate(lgs=if(first(real)==1L)as.integer(grade%in%c("C","D"))else{u<-runif(n());as.integer(u==max(u))})%>%ungroup()
  gap<-side%>%group_by(boundary_id,real)%>%summarise(across(all_of(cvars),~mean(.x[lgs==1],na.rm=T)-mean(.x[lgs==0],na.rm=T)),.groups="drop")
  # probit propensity of real vs placebo on z-scored covariate gaps (+ missingness indicators), ALL covariates
  g<-gap; rhs<-c()
  for(v in cvars){s<-as.numeric(scale(g[[v]])); g[[paste0("R_",v)]]<-ifelse(is.na(s),0,s); g[[paste0("MI_",v)]]<-as.integer(is.na(s)); rhs<-c(rhs,paste0("R_",v),paste0("MI_",v))}
  rhs<-rhs[sapply(rhs,function(c)sd(g[[c]],na.rm=T)>0)]
  ps<-suppressWarnings(glm(reformulate(rhs,"real"),data=g,family=binomial("probit"))); g$phat<-predict(ps,type="response")
  # ATT odds weights on placebos (den. capped at .02), normalized to mean 1 within placebos; real=1
  odds<-g$phat/pmax(1-g$phat,.02)
  g<-g%>%mutate(w_att=ifelse(real==1L,1,odds/mean(odds[real==0])))
  # balance check on covariate gaps under IPW weights (compare to EB's ~0)
  bvars<-cvars[sapply(cvars,function(v)sum(is.finite(as.numeric(scale(gap[[v]]))))>100)]
  bal<-max(sapply(bvars,function(v){z<-as.numeric(scale(gap[[v]]));ok<-is.finite(z);
    abs(mean(z[ok&g$real==1])-weighted.mean(z[ok&g$real==0],g$w_att[ok&g$real==0]))}),na.rm=T)
  dat<-side%>%select(boundary_id,boundary_side,lgs,real)%>%inner_join(g%>%select(boundary_id,w_att),by="boundary_id")%>%
    inner_join(bt%>%distinct(boundary_id,boundary_side,tract11),by=c("boundary_id","boundary_side"))%>%
    inner_join(poi,by="tract11",relationship="many-to-many")%>%mutate(w=w_att*lvis)
  # feols(seg ~ lgs + lgs:real): coef(lgs)=placebo (grid) discontinuity;
  # coef(lgs)+coef(lgs:real)=real-boundary discontinuity; coef(lgs:real)=DiD.
  m<-suppressMessages(feols(as.formula(paste0(yv,"~lgs+lgs:real")),dat,fixef=c("boundary_id","sub_category"),weights=~w,cluster=~boundary_id))
  tibble(seed=seed,placebo=unname(coef(m)["lgs"]),real=unname(coef(m)["lgs"]+coef(m)["lgs:real"]),
    did=unname(coef(m)["lgs:real"]),did_se=unname(se(m)["lgs:real"]),plac_se=unname(se(m)["lgs"]),
    p=unname(pvalue(m)["lgs:real"]),bal=bal)}
report<-function(r,lab){
  cat(sprintf("\n== %s ==\n",lab))
  cat("                       estimate    SE(asymptotic)   SE(randomization)\n")
  cat(sprintf("Real boundary effect  %+.5f        --            %.5f\n",mean(r$real),sd(r$real)))
  cat(sprintf("Placebo (grid) effect %+.5f     %.5f          %.5f\n",mean(r$placebo),mean(r$plac_se),sd(r$placebo)))
  cat(sprintf("DiD (real - placebo)  %+.5f     %.5f          %.5f\n",mean(r$did),mean(r$did_se),sd(r$did)))
  cat(sprintf("95%% randomization CI: [%+.5f, %+.5f] %s | %% draws did>0=%.0f%% | %% p<.05=%.0f%% | maxbal=%.4f\n",
    quantile(r$did,.025),quantile(r$did,.975),ifelse(quantile(r$did,.025)>0,"EXCLUDES 0","incl 0"),
    100*mean(r$did>0),100*mean(r$p<.05),max(r$bal)))}
cat("Probit-IPW CD-B difference-in-discontinuity on RECOMPUTED measures, 100 draws.\n")
cat("Difference between the two rows below is the SafeGraph post-stratification reweighting.\n")
r_unw<-map_dfr(1:100,estimate,yv="seg_unw"); report(r_unw,"Unweighted recompute (seg_int_unw)")
r_rw <-map_dfr(1:100,estimate,yv="seg_rw");  report(r_rw, "Reweighted (seg_int_rw)")
saveRDS(r_unw,"results/holc_cdb_ipw_unw_draws.rds"); saveRDS(r_rw,"results/holc_cdb_ipw_rw_draws.rds")
