# Race-income INTERSECTIONALITY check on the entropy-balanced CD-B boundary DiD.
# Does redlining also raise INCOME activity-space segregation, and does the RACIAL effect survive
# netting out income? Two income measures (per the two imputations):
#   (1) visitor-income: inc_seg_msa_adj (income analog of the racial measure) + mean visitor quartile;
#   (2) POI-CBG income: median household income of the POI's neighborhood (nbhd_linc, standardized log).
# Same entropy-balanced boundary design/weights/clustering as the main estimator (holc_cdb_ebal_rw.R).
suppressMessages({library(tidyverse); library(haven); library(fixest)})
ahm_path<-"./data/ahm_staging/"; a2f<-function(x)paste0(substr(x,2,3),substr(x,5,7),substr(x,9,14))
BAR<-c("BUFFintrain","BUFFinsrivers","BUFFinbrivers")
covset<-c("black_","ownhome_","foreign_born_","read_write_","house_value_","rent_","radio_"); yrs<-c(1910,1920,1930)
load("./data/safegraph_2019_m.rdata")
poi<-places_usa_2019%>%
  left_join(readRDS("data/poi_seg_fraccorr_as_rw.rds")%>%transmute(safegraph_place_id,seg_int_rw),by="safegraph_place_id")%>%
  mutate(seg_race=seg_int_rw*5/8, seg_inc=inc_seg_msa_adj_mean,
         nbhd_linc=as.numeric(scale(log(median_hinc))),
         vinc=1*v_q1_rate+2*v_q2_rate+3*v_q3_rate+4*v_q4_rate,   # mean visitor income quartile (1--4)
         lvis=log(raw_visit_counts+1), tract11=substr(poi_cbg,1,11))%>%
  filter(!is.na(seg_race),!is.na(seg_inc),!is.na(lvis))%>%
  select(tract11,sub_category,seg_race,seg_inc,nbhd_linc,vinc,lvis)
read_map<-function(file,idcol,pfx,is_grid){d<-read_dta(paste0(ahm_path,file))%>%filter(year==2010)%>%mutate(tract11=a2f(trctID2010))
  if(is_grid) d<-d%>%rename_with(~sub("BUFF_CFin","BUFFin",.x))
  d%>%select(tract11,matches(paste0("^",idcol,"[0-9]+$")),matches("^boundary_side[0-9]+$|^side_grid[0-9]+$"),matches("^holc_grade[0-9]+$"),
      matches("^boundary_grades[0-9]+$"),matches("^BUFFintrain[0-9]+$"),matches("^BUFFinsrivers[0-9]+$"),matches("^BUFFinbrivers[0-9]+$"))%>%
    rename_with(~sub("^side_grid","boundary_side",.x))%>%rename_with(~sub(paste0("^",idcol),"boundary_id",.x))%>%
    pivot_longer(-tract11,names_to=c(".value","slot"),names_pattern="^(boundary_id|boundary_side|holc_grade|boundary_grades|BUFFintrain|BUFFinsrivers|BUFFinbrivers)([0-9]+)$")%>%
    drop_na(boundary_id)%>%mutate(boundary_id=paste0(pfx,boundary_id))}
mA<-read_map("redlining_dataset_t2010_4th_new.dta","boundary_id","A_",FALSE)
mG<-read_map("redlining_dataset_tgrid2010_4th_new.dta","boundary_id_grid","G_",TRUE)
raw<-read_dta(paste0(ahm_path,"redlining_dataset_t2010_4th_new.dta"))%>%filter(year%in%yrs)%>%mutate(tract11=a2f(trctID2010))
cov_all<-raw%>%select(tract11,year,any_of(covset))%>%group_by(tract11,year)%>%summarise(across(everything(),~mean(.x,na.rm=T)),.groups="drop")%>%pivot_wider(names_from=year,values_from=any_of(covset),names_glue="{.value}{year}");cvars<-setdiff(names(cov_all),"tract11")
tr<-mA%>%filter(boundary_grades%in%c("2_to_3","2_to_4"))%>%mutate(real=1L);pl<-mG%>%filter(boundary_grades%in%c("2_to_2","3_to_3","4_to_4"))%>%mutate(real=0L)
bt<-bind_rows(tr,pl)%>%left_join(cov_all,by="tract11")%>%mutate(barrier_amt=rowSums(across(all_of(BAR),~replace_na(.x,0))))
bt<-bt%>%semi_join(bt%>%group_by(boundary_id)%>%summarise(b=as.integer(max(barrier_amt,na.rm=T)>0.05))%>%filter(b==0),by="boundary_id")
ebalance<-function(X,target,maxit=300,tol=1e-8){Z<-sweep(X,2,target,"-");n<-nrow(Z);base<-rep(1/n,n);lam<-rep(0,ncol(Z));w<-base
  for(it in 1:maxit){u<-as.vector(Z%*%lam);u<-u-max(u);w<-base*exp(u);w<-w/sum(w);g<-colSums(w*Z);if(max(abs(g))<tol)break
    H<-t(Z)%*%(w*Z)-outer(g,g);step<-tryCatch(solve(H+diag(1e-8,ncol(Z)),g),error=function(e)g);best<-max(abs(g));lam0<-lam
    for(s in c(1,.5,.25,.1,.03)){lt<-lam0-s*step;ut<-as.vector(Z%*%lt);ut<-ut-max(ut);wt<-base*exp(ut);wt<-wt/sum(wt);gt<-colSums(wt*Z);if(max(abs(gt))<best){lam<-lt;best<-max(abs(gt));break}};if(identical(lam,lam0))break};w}
key<-intersect(paste0(c("black_","ownhome_","house_value_","foreign_born_"),"1930"),cvars)
mkdat<-function(seed){side<-bt%>%group_by(boundary_id,real,boundary_side)%>%summarise(grade=first(holc_grade[holc_grade%in%c("A","B","C","D")]),across(all_of(cvars),~mean(.x,na.rm=T)),.groups="drop")%>%group_by(boundary_id)%>%filter(n()==2)%>%ungroup()
  set.seed(seed);side<-side%>%group_by(boundary_id)%>%mutate(lgs=if(first(real)==1L)as.integer(grade%in%c("C","D"))else{u<-runif(n());as.integer(u==max(u))})%>%ungroup()
  gap<-side%>%group_by(boundary_id,real)%>%summarise(across(all_of(cvars),~mean(.x[lgs==1],na.rm=T)-mean(.x[lgs==0],na.rm=T)),.groups="drop")
  bvars<-cvars[sapply(cvars,function(v)sum(is.finite(as.numeric(scale(gap[[v]]))))>100)]
  gs<-gap[,"real"];cons<-c()
  for(v in bvars){z<-as.numeric(scale(gap[[v]]));gs[[paste0("Z_",v)]]<-ifelse(is.finite(z),z,0);cons<-c(cons,paste0("Z_",v))
    if(any(!is.finite(z))){gs[[paste0("MI_",v)]]<-as.numeric(!is.finite(z));cons<-c(cons,paste0("MI_",v))}}
  tgt<-colMeans(gs[gs$real==1,cons,drop=FALSE]);wp<-ebalance(as.matrix(gs[gs$real==0,cons,drop=FALSE]),tgt)*sum(gs$real==0)
  gg<-gap%>%mutate(w_att=ifelse(real==1L,1,NA));gg$w_att[gg$real==0]<-wp
  side%>%select(boundary_id,boundary_side,lgs,real)%>%inner_join(gg%>%select(boundary_id,w_att),by="boundary_id")%>%
    inner_join(bt%>%distinct(boundary_id,boundary_side,tract11),by=c("boundary_id","boundary_side"))%>%inner_join(poi,by="tract11",relationship="many-to-many")%>%mutate(w=w_att*lvis)}
g<-function(m)c(b=unname(coef(m)["lgs:real"]),s=unname(se(m)["lgs:real"]))
run<-function(seed){d<-mkdat(seed)
  race<-feols(seg_race~lgs+lgs:real|boundary_id+sub_category,d,weights=~w,cluster=~boundary_id)
  inc <-feols(seg_inc ~lgs+lgs:real|boundary_id+sub_category,d,weights=~w,cluster=~boundary_id)
  vinc<-feols(vinc   ~lgs+lgs:real|boundary_id+sub_category,d,weights=~w,cluster=~boundary_id)
  rcn <-feols(seg_race~lgs+lgs:real+nbhd_linc|boundary_id+sub_category,d,weights=~w,cluster=~boundary_id)  # ctrl POI-CBG income
  rci <-feols(seg_race~lgs+lgs:real+seg_inc  |boundary_id+sub_category,d,weights=~w,cluster=~boundary_id)  # ctrl income-seg
  tibble(seed=seed,
    race_b=g(race)["b"],race_s=g(race)["s"],
    inc_b=g(inc)["b"],inc_s=g(inc)["s"],
    vinc_b=g(vinc)["b"],vinc_s=g(vinc)["s"],
    rcn_b=g(rcn)["b"],rcn_s=g(rcn)["s"],
    rci_b=g(rci)["b"],rci_s=g(rci)["s"],
    n_race=nobs(race),n_rcn=nobs(rcn))}
R<-map_dfr(1:25,run); s<-R%>%summarise(across(-seed,mean))
pf<-function(b,se)sprintf("%+.5f (se %.5f, t=%.2f, p=%.4f)",b,se,b/se,2*pnorm(-abs(b/se)))
cat("=== Race-income intersectionality, entropy-balanced boundary DiD (avg over 25 placebo orientations) ===\n\n")
cat(sprintf("1. Redlining -> RACIAL integration segregation (baseline):   %s\n",pf(s$race_b,s$race_s)))
cat(sprintf("2. Redlining -> INCOME segregation (visitor income vs MSA):   %s\n",pf(s$inc_b,s$inc_s)))
cat(sprintf("3. Redlining -> mean visitor income quartile (1-4):           %s\n",pf(s$vinc_b,s$vinc_s)))
cat(sprintf("4. RACIAL effect, controlling POI-CBG income (nbhd_linc):     %s   [n=%.0f vs %.0f]\n",pf(s$rcn_b,s$rcn_s),s$n_rcn,s$n_race))
cat(sprintf("5. RACIAL effect, controlling income-segregation (seg_inc):   %s\n",pf(s$rci_b,s$rci_s)))
cat(sprintf("\n   Racial effect retained after income controls: %.0f%% (nbhd income), %.0f%% (income seg)\n",100*s$rcn_b/s$race_b,100*s$rci_b/s$race_b))
saveRDS(R,"results/holc_income_intersect_FRACCORR.rds")
