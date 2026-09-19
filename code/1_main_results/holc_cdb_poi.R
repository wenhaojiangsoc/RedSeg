# POI-LEVEL boundary difference-in-discontinuity (Steps 2-3).
# Real boundaries: adjacent B<->C/D polygon pairs; POIs assigned by their own distance to the
#   border within 1/4 mile (holc_poi_assign.R). Placebo: 1/2-mile grid, same-grade cell edges,
#   POIs within 1/4 mile split by side, random orientation (RSUE, 100 draws).
# Balancing: 1910-1930 covariates (attached per POI via its tract11), EB and IPW, EB_SET="all".
# Outcome: MY recomputed segregation, unweighted (seg_int_unw) and SafeGraph-reweighted (seg_int_rw).
# Reports EB & IPW x unweighted & reweighted, all as POI-level DiD.
suppressMessages({library(tidyverse); library(haven); library(fixest)})
# run from the project root; # setwd("~/Library/CloudStorage/Dropbox/Demography-Rep-2026")
GRID_M<-804.672; BUFFER_M<-402.336
covset<-c("black_","ownhome_","foreign_born_","read_write_","house_value_","rent_","radio_"); yrs<-c(1910,1920,1930)
a2f<-function(x)paste0(substr(x,2,3),substr(x,5,7),substr(x,9,14)); ahm_path<-"./data/ahm_staging/"

# ---- 1910-1930 covariates per 2010 tract ----
raw<-read_dta(paste0(ahm_path,"redlining_dataset_t2010_4th_new.dta"))%>%filter(year%in%yrs)%>%mutate(tract11=a2f(trctID2010))
cov_all<-raw%>%select(tract11,year,any_of(covset))%>%group_by(tract11,year)%>%summarise(across(everything(),~mean(.x,na.rm=T)),.groups="drop")%>%
  pivot_wider(names_from=year,values_from=any_of(covset),names_glue="{.value}{year}"); cvars<-setdiff(names(cov_all),"tract11")

# ---- real POIs (low side = C/D by construction) ----
real<-readRDS("data/poi_real_assign.rds")%>%
  transmute(safegraph_place_id,boundary_id=paste0("R",boundary_id),real=1L,lgs=as.integer(side=="low"),
            gside=NA_real_,grade,tract11,sub_category,lvis,seg_unw=seg_int_unw*5/8,seg_rw=seg_int_rw*5/8)

# ---- placebo POIs: nearest 1/2-mile grid cell-edge within 1/4 mi, same-grade, both sides present ----
pts<-readRDS("data/poi_bcd_pts.rds"); G<-GRID_M
kx<-round(pts$X/G)*G; dxv<-pts$X-kx; ky<-round(pts$Y/G)*G; dyh<-pts$Y-ky
vert<-abs(dxv)<=abs(dyh)
plac<-pts%>%mutate(
    boundary_id=paste0("P",ifelse(vert,"V","H"),"_",ifelse(vert,round(X/G),round(Y/G)),"_",ifelse(vert,floor(Y/G),floor(X/G))),
    gside=ifelse(vert,sign(dxv),sign(dyh)), gdist=ifelse(vert,abs(dxv),abs(dyh)))%>%
  filter(gdist<=BUFFER_M, gside!=0)
valid<-plac%>%group_by(boundary_id)%>%summarise(ns=n_distinct(gside),ng=n_distinct(grade),.groups="drop")%>%
  filter(ns==2,ng==1)%>%pull(boundary_id)
plac<-plac%>%filter(boundary_id%in%valid)%>%
  transmute(safegraph_place_id,boundary_id,real=0L,lgs=NA_integer_,gside,grade,tract11,sub_category,lvis,
            seg_unw=seg_int_unw*5/8,seg_rw=seg_int_rw*5/8)
cat(sprintf("Real: %d POIs / %d boundaries | Placebo (same-grade grid): %d POIs / %d boundaries\n",
    nrow(real),n_distinct(real$boundary_id),nrow(plac),n_distinct(plac$boundary_id)))

allpoi<-bind_rows(real,plac)%>%left_join(cov_all,by="tract11")
plac_bids<-unique(plac$boundary_id)

# ---- entropy balancing (same solver as holc_cdb_ebal.R) ----
ebalance<-function(X,target,maxit=300,tol=1e-8){Z<-sweep(X,2,target,"-");n<-nrow(Z);base<-rep(1/n,n);lam<-rep(0,ncol(Z));w<-base
  for(it in 1:maxit){u<-as.vector(Z%*%lam);u<-u-max(u);w<-base*exp(u);w<-w/sum(w);g<-colSums(w*Z);if(max(abs(g))<tol)break
    H<-t(Z)%*%(w*Z)-outer(g,g);step<-tryCatch(solve(H+diag(1e-8,ncol(Z)),g),error=function(e)g);best<-max(abs(g));lam0<-lam
    for(s in c(1,.5,.25,.1,.03)){lt<-lam0-s*step;ut<-as.vector(Z%*%lt);ut<-ut-max(ut);wt<-base*exp(ut);wt<-wt/sum(wt);gt<-colSums(wt*Z);if(max(abs(gt))<best){lam<-lt;best<-max(abs(gt));break}};if(identical(lam,lam0))break};w}

draw<-function(seed){
  set.seed(seed)
  po<-tibble(boundary_id=plac_bids,low_side=sample(c(-1,1),length(plac_bids),replace=TRUE))
  d<-allpoi%>%left_join(po,by="boundary_id")%>%
    mutate(lgs=ifelse(real==1L,lgs,as.integer(gside==low_side)))
  # per-boundary covariate gaps (low - high), keep boundaries with both sides
  bs<-d%>%group_by(boundary_id,real)%>%summarise(n0=sum(lgs==0),n1=sum(lgs==1),
        across(all_of(cvars),~mean(.x[lgs==1],na.rm=T)-mean(.x[lgs==0],na.rm=T)),.groups="drop")%>%
     filter(n0>0,n1>0)
  bvars<-cvars[sapply(cvars,function(v)sum(is.finite(as.numeric(scale(bs[[v]]))))>100)]
  gs<-bs[,c("boundary_id","real")]; cons<-c()
  for(v in bvars){z<-as.numeric(scale(bs[[v]]));gs[[paste0("Z_",v)]]<-ifelse(is.finite(z),z,0);cons<-c(cons,paste0("Z_",v))
    if(any(!is.finite(z))){gs[[paste0("MI_",v)]]<-as.numeric(!is.finite(z));cons<-c(cons,paste0("MI_",v))}}
  tgt<-colMeans(gs[gs$real==1,cons,drop=FALSE])
  Xp<-as.matrix(gs[gs$real==0,cons,drop=FALSE])
  wp<-ebalance(Xp,tgt)*sum(gs$real==0)                     # EB placebo weights
  ps<-suppressWarnings(glm(reformulate(cons,"real"),data=gs,family=binomial("probit")));ph<-predict(ps,type="response")
  odds<-ph/pmax(1-ph,.02); w_ipw<-ifelse(gs$real==1,1,odds/mean(odds[gs$real==0]))   # IPW ATT weights
  bw<-gs%>%transmute(boundary_id,w_eb=ifelse(real==1,1,NA_real_),w_ipw=w_ipw); bw$w_eb[gs$real==0]<-wp
  eb_bal<-max(sapply(bvars,function(v){z<-as.numeric(scale(bs[[v]]));ok<-is.finite(z)
    abs(mean(z[ok&gs$real==1])-weighted.mean(z[ok&gs$real==0],wp[ok[gs$real==0]]))}),na.rm=T)
  dd<-d%>%inner_join(bw,by="boundary_id")
  fit<-function(yv,wcol){m<-suppressMessages(feols(as.formula(paste0(yv,"~lgs+lgs:real")),dd,
        fixef=c("boundary_id","sub_category"),weights=as.formula(paste0("~I(",wcol,"*lvis)")),cluster=~boundary_id))
    tibble(placebo=unname(coef(m)["lgs"]),real=unname(coef(m)["lgs"]+coef(m)["lgs:real"]),
           did=unname(coef(m)["lgs:real"]),did_se=unname(se(m)["lgs:real"]),p=unname(pvalue(m)["lgs:real"]))}
  bind_rows(
    fit("seg_unw","w_eb") %>%mutate(scheme="EB", measure="Recomputed (unwtd)"),
    fit("seg_rw","w_eb")  %>%mutate(scheme="EB", measure="Recomputed (reweighted)"),
    fit("seg_unw","w_ipw")%>%mutate(scheme="IPW",measure="Recomputed (unwtd)"),
    fit("seg_rw","w_ipw") %>%mutate(scheme="IPW",measure="Recomputed (reweighted)"))%>%
    mutate(seed=seed,eb_bal=eb_bal)}

cat("running 100 random-orientation draws (POI-level DiD)...\n")
R<-map_dfr(1:100,draw)
saveRDS(R,"results/holc_cdb_poi_draws.rds")
tab<-R%>%group_by(scheme,measure)%>%summarise(
  placebo=mean(placebo),real=mean(real),did=mean(did),se_asy=mean(did_se),
  ci_lo=quantile(did,.025),ci_hi=quantile(did,.975),pos=mean(did>0)*100,p05=mean(p<.05)*100,.groups="drop")
cat("\n=== POI-LEVEL boundary difference-in-discontinuity ===\n")
cat(sprintf("%-5s %-24s %9s %9s %-22s %5s %6s\n","Sch","Measure","Real","DiD","95% rand CI","%>0","%p.05"))
for(i in seq_len(nrow(tab)))with(tab[i,],cat(sprintf("%-5s %-24s %+9.5f %+9.5f  [%+.5f,%+.5f] %5.0f %6.0f\n",
  scheme,measure,real,did,ci_lo,ci_hi,pos,p05)))
cat(sprintf("\nMean EB max covariate imbalance across draws: %.4f\n",mean(R$eb_bal)))
