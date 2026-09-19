# Appendix robustness: EB boundary DiD (reweighted measure) under alternative buffer width (1/8 mi)
# and alternative covariate set (key-4), plus category heterogeneity. Loads rdata once.
suppressMessages({library(tidyverse); library(haven); library(fixest)})
# run from the project root; # setwd("~/Library/CloudStorage/Dropbox/Demography-Rep-2026")
ahm<-"./data/ahm_staging/"; a2f<-function(x)paste0(substr(x,2,3),substr(x,5,7),substr(x,9,14))
BAR<-c("BUFFintrain","BUFFinsrivers","BUFFinbrivers")
covset<-c("black_","ownhome_","foreign_born_","read_write_","house_value_","rent_","radio_"); yrs<-c(1910,1920,1930)
load("./data/safegraph_2019_m.rdata"); rw<-readRDS("data/poi_seg_frac_as_rw.rds")
POI<-places_usa_2019%>%left_join(rw%>%transmute(safegraph_place_id,seg_rw=seg_int_rw*5/8),by="safegraph_place_id")%>%
  mutate(lvis=log(raw_visit_counts+1),tract11=substr(poi_cbg,1,11),
         daily=as.integer(grepl("Grocery|Restaurant|Food|Gasoline|Pharmac|Drug|Retail|Convenience|Full-Service|Limited-Service|Supermarket|Department",sub_category,ignore.case=TRUE)))%>%
  filter(!is.na(seg_rw),!is.na(lvis))
read_map<-function(file,idcol,pfx,is_grid){d<-read_dta(paste0(ahm,file))%>%filter(year==2010)%>%mutate(tract11=a2f(trctID2010))
  if(is_grid) d<-d%>%rename_with(~sub("BUFF_CFin","BUFFin",.x))
  d%>%select(tract11,matches(paste0("^",idcol,"[0-9]+$")),matches("^boundary_side[0-9]+$|^side_grid[0-9]+$"),matches("^holc_grade[0-9]+$"),
      matches("^boundary_grades[0-9]+$"),matches("^BUFFintrain[0-9]+$"),matches("^BUFFinsrivers[0-9]+$"),matches("^BUFFinbrivers[0-9]+$"))%>%
    rename_with(~sub("^side_grid","boundary_side",.x))%>%rename_with(~sub(paste0("^",idcol),"boundary_id",.x))%>%
    pivot_longer(-tract11,names_to=c(".value","slot"),names_pattern="^(boundary_id|boundary_side|holc_grade|boundary_grades|BUFFintrain|BUFFinsrivers|BUFFinbrivers)([0-9]+)$")%>%
    drop_na(boundary_id)%>%mutate(boundary_id=paste0(pfx,boundary_id))}
ebalance<-function(X,target,maxit=300,tol=1e-8){Z<-sweep(X,2,target,"-");n<-nrow(Z);base<-rep(1/n,n);lam<-rep(0,ncol(Z));w<-base
  for(it in 1:maxit){u<-as.vector(Z%*%lam);u<-u-max(u);w<-base*exp(u);w<-w/sum(w);g<-colSums(w*Z);if(max(abs(g))<tol)break
    H<-t(Z)%*%(w*Z)-outer(g,g);step<-tryCatch(solve(H+diag(1e-8,ncol(Z)),g),error=function(e)g);best<-max(abs(g));lam0<-lam
    for(s in c(1,.5,.25,.1,.03)){lt<-lam0-s*step;ut<-as.vector(Z%*%lt);ut<-ut-max(ut);wt<-base*exp(ut);wt<-wt/sum(wt);gt<-colSums(wt*Z);if(max(abs(gt))<best){lam<-lt;best<-max(abs(gt));break}};if(identical(lam,lam0))break};w}
run<-function(tf,gf,ebset,catf,ndraw=100){
  mA<-read_map(tf,"boundary_id","A_",FALSE); mG<-read_map(gf,"boundary_id_grid","G_",TRUE)
  raw<-read_dta(paste0(ahm,tf))%>%filter(year%in%yrs)%>%mutate(tract11=a2f(trctID2010))
  cov_all<-raw%>%select(tract11,year,any_of(covset))%>%group_by(tract11,year)%>%summarise(across(everything(),~mean(.x,na.rm=T)),.groups="drop")%>%
    pivot_wider(names_from=year,values_from=any_of(covset),names_glue="{.value}{year}");cvars<-setdiff(names(cov_all),"tract11")
  key<-intersect(paste0(c("black_","ownhome_","house_value_","foreign_born_"),"1930"),cvars)
  poi<-if(is.null(catf)) POI else POI%>%filter(daily==catf)
  poi<-poi%>%select(tract11,sub_category,seg_rw,lvis)
  tr<-mA%>%filter(boundary_grades%in%c("2_to_3","2_to_4"))%>%mutate(real=1L)
  pl<-mG%>%filter(boundary_grades%in%c("2_to_2","3_to_3","4_to_4"))%>%mutate(real=0L)
  bt<-bind_rows(tr,pl)%>%left_join(cov_all,by="tract11")%>%mutate(ba=rowSums(across(all_of(BAR),~replace_na(.x,0))))
  bt<-bt%>%semi_join(bt%>%group_by(boundary_id)%>%summarise(b=as.integer(max(ba,na.rm=T)>0.05))%>%filter(b==0),by="boundary_id")
  est<-function(seed){side<-bt%>%group_by(boundary_id,real,boundary_side)%>%summarise(grade=first(holc_grade[holc_grade%in%c("A","B","C","D")]),across(all_of(cvars),~mean(.x,na.rm=T)),.groups="drop")%>%group_by(boundary_id)%>%filter(n()==2)%>%ungroup()
    set.seed(seed);side<-side%>%group_by(boundary_id)%>%mutate(lgs=if(first(real)==1L)as.integer(grade%in%c("C","D"))else{u<-runif(n());as.integer(u==max(u))})%>%ungroup()
    gap<-side%>%group_by(boundary_id,real)%>%summarise(across(all_of(cvars),~mean(.x[lgs==1],na.rm=T)-mean(.x[lgs==0],na.rm=T)),.groups="drop")
    bvars<-if(ebset=="all") cvars[sapply(cvars,function(v)sum(is.finite(as.numeric(scale(gap[[v]]))))>100)] else key
    gs<-gap[,"real"];cons<-c()
    for(v in bvars){z<-as.numeric(scale(gap[[v]]));gs[[paste0("Z_",v)]]<-ifelse(is.finite(z),z,0);cons<-c(cons,paste0("Z_",v))
      if(ebset=="all"&&any(!is.finite(z))){gs[[paste0("MI_",v)]]<-as.numeric(!is.finite(z));cons<-c(cons,paste0("MI_",v))}}
    tgt<-colMeans(gs[gs$real==1,cons,drop=FALSE]);wp<-ebalance(as.matrix(gs[gs$real==0,cons,drop=FALSE]),tgt)*sum(gs$real==0)
    g<-gap%>%mutate(w_att=ifelse(real==1L,1,NA));g$w_att[g$real==0]<-wp
    dat<-side%>%select(boundary_id,boundary_side,lgs,real)%>%inner_join(g%>%select(boundary_id,w_att),by="boundary_id")%>%
      inner_join(bt%>%distinct(boundary_id,boundary_side,tract11),by=c("boundary_id","boundary_side"))%>%inner_join(poi,by="tract11",relationship="many-to-many")%>%mutate(w=w_att*lvis)
    m<-suppressMessages(feols(seg_rw~lgs+lgs:real|boundary_id+sub_category,dat,weights=~w,cluster=~boundary_id))
    c(did=unname(coef(m)["lgs:real"]),se=unname(se(m)["lgs:real"]),n=nobs(m))}
  R<-map_dfr(1:ndraw,function(s){v<-est(s);tibble(did=v["did"],se=v["se"],n=v["n"])})
  sprintf("did=%+.5f  se=%.5f  t=%.2f  p=%.3f  Nobs~%d",mean(R$did),mean(R$se),mean(R$did)/mean(R$se),2*pnorm(-abs(mean(R$did)/mean(R$se))),round(mean(R$n)))}
cat("Baseline (1/4 mi, all cov):   ",run("redlining_dataset_t2010_4th_new.dta","redlining_dataset_tgrid2010_4th_new.dta","all",NULL),"\n")
cat("Key-4 covariates (1/4 mi):    ",run("redlining_dataset_t2010_4th_new.dta","redlining_dataset_tgrid2010_4th_new.dta","key4",NULL),"\n")
cat("1/8-mile buffer (all cov):    ",run("redlining_dataset_t2010_8th_new.dta","redlining_dataset_tgrid2010_8th_new.dta","all",NULL),"\n")
cat("Daily-routine POIs (1/4,all): ",run("redlining_dataset_t2010_4th_new.dta","redlining_dataset_tgrid2010_4th_new.dta","all",1),"\n")
cat("Infrequent POIs (1/4,all):    ",run("redlining_dataset_t2010_4th_new.dta","redlining_dataset_tgrid2010_4th_new.dta","all",0),"\n")
