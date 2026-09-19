# AHM low-propensity-score design at the POI LEVEL (POIs assigned by their own distance to the
# B<->C/D polygon border, holc_poi_assign.R -- NOT AHM tracts). More borders, tighter buffer.
# Propensity = probit(real border | pre-map 1910-1930 cross-border gaps), fit on real + same-grade
# grid placebos. Among REAL borders we take the more-balanced HALF -- those with the smallest OBSERVED
# pre-map gaps (near-zero across Black share, homeownership, house value, foreign-born share), i.e. where
# grading was as-good-as-random. On those we run the SIMPLE border comparison (no placebo, no EB):
#   seg_ik = beta_lowp * 1[lower-graded side] + boundary FE + category FE.
# Prints the covariate-balance table and a cutoff sweep (estimate + balance quality at 25-50%).
suppressMessages({library(tidyverse); library(haven); library(fixest)})
GRID_M<-804.672; BUFFER_M<-402.336
covset<-c("black_","ownhome_","foreign_born_","read_write_","house_value_","rent_","radio_"); yrs<-c(1910,1920,1930)
a2f<-function(x)paste0(substr(x,2,3),substr(x,5,7),substr(x,9,14)); ahm_path<-"./data/ahm_staging/"
raw<-read_dta(paste0(ahm_path,"redlining_dataset_t2010_4th_new.dta"))%>%filter(year%in%yrs)%>%mutate(tract11=a2f(trctID2010))
cov_all<-raw%>%select(tract11,year,any_of(covset))%>%group_by(tract11,year)%>%summarise(across(everything(),~mean(.x,na.rm=T)),.groups="drop")%>%
  pivot_wider(names_from=year,values_from=any_of(covset),names_glue="{.value}{year}"); cvars<-setdiff(names(cov_all),"tract11")
FRC<-dplyr::as_tibble(readRDS("results/poi_ei_bounds_frac.rds"))%>%dplyr::transmute(safegraph_place_id,seg_frac=seg_eco*5/8)
real<-readRDS("data/poi_real_assign.rds")%>%inner_join(FRC,by="safegraph_place_id")%>%transmute(safegraph_place_id,boundary_id=paste0("R",boundary_id),real=1L,lgs=as.integer(side=="low"),
            gside=NA_real_,grade,tract11,sub_category,lvis,seg_rw=seg_frac)
pts<-readRDS("data/poi_bcd_pts.rds"); G<-GRID_M
kx<-round(pts$X/G)*G; dxv<-pts$X-kx; ky<-round(pts$Y/G)*G; dyh<-pts$Y-ky; vert<-abs(dxv)<=abs(dyh)
plac<-pts%>%mutate(boundary_id=paste0("P",ifelse(vert,"V","H"),"_",ifelse(vert,round(X/G),round(Y/G)),"_",ifelse(vert,floor(Y/G),floor(X/G))),
    gside=ifelse(vert,sign(dxv),sign(dyh)),gdist=ifelse(vert,abs(dxv),abs(dyh)))%>%filter(gdist<=BUFFER_M,gside!=0)
valid<-plac%>%group_by(boundary_id)%>%summarise(ns=n_distinct(gside),ng=n_distinct(grade),.groups="drop")%>%filter(ns==2,ng==1)%>%pull(boundary_id)
plac<-plac%>%filter(boundary_id%in%valid)%>%inner_join(FRC,by="safegraph_place_id")%>%transmute(safegraph_place_id,boundary_id,real=0L,lgs=NA_integer_,gside,grade,tract11,sub_category,lvis,seg_rw=seg_frac)
allpoi<-bind_rows(real,plac)%>%left_join(cov_all,by="tract11"); plac_bids<-unique(plac$boundary_id)
cat(sprintf("Real: %d POIs / %d boundaries | Placebo: %d POIs / %d boundaries\n\n",nrow(real),n_distinct(real$boundary_id),nrow(plac),n_distinct(plac$boundary_id)))
psv<-intersect(paste0(rep(c("black_","ownhome_","house_value_","foreign_born_","rent_"),each=3),c(1910,1920,1930)),cvars)
run<-function(seed){set.seed(seed)
  po<-tibble(boundary_id=plac_bids,low_side=sample(c(-1,1),length(plac_bids),replace=TRUE))
  d<-allpoi%>%left_join(po,by="boundary_id")%>%mutate(lgs=ifelse(real==1L,lgs,as.integer(gside==low_side)))
  bs<-d%>%group_by(boundary_id,real)%>%summarise(n0=sum(lgs==0),n1=sum(lgs==1),across(all_of(cvars),~mean(.x[lgs==1],na.rm=T)-mean(.x[lgs==0],na.rm=T)),.groups="drop")%>%filter(n0>0,n1>0)
  X<-c(); for(v in psv){z<-as.numeric(scale(bs[[v]]));bs[[paste0("Z_",v)]]<-ifelse(is.finite(z),z,0);X<-c(X,paste0("Z_",v));if(any(!is.finite(z))){bs[[paste0("MI_",v)]]<-as.numeric(!is.finite(z));X<-c(X,paste0("MI_",v))}}
  keycov<-intersect(paste0(c("black_","ownhome_","house_value_","foreign_born_"),"1930"),cvars)
  rl<-bs%>%filter(real==1)%>%filter(if_all(all_of(keycov),is.finite))   # OBSERVED key 1930 covariates only
  zc<-paste0("z_",keycov); for(v in keycov) rl[[paste0("z_",v)]]<-as.numeric(scale(rl[[v]]))
  rl$gnorm<-sqrt(rowSums(as.matrix(rl[,zc])^2))                 # magnitude of OBSERVED pre-map gaps (0 = balanced)
  bal_ids<-rl$boundary_id[rl$gnorm<=quantile(rl$gnorm,1/2)]     # more-balanced half (as-good-as-random)
  imb_ids<-rl$boundary_id[rl$gnorm>quantile(rl$gnorm,1/2)]      # less-balanced half
  all_ids<-bs$boundary_id[bs$real==1]
  reg<-function(ids){dd<-d%>%filter(real==1,boundary_id%in%ids)
    m<-feols(seg_rw~lgs|boundary_id+sub_category,dd,weights=~lvis,cluster=~boundary_id)
    c(est=unname(coef(m)["lgs"]),se=unname(se(m)["lgs"]),nb=n_distinct(dd$boundary_id),npoi=nrow(dd))}
  a<-reg(all_ids); b<-reg(bal_ids); im<-reg(imb_ids)
  # cutoff sweep for the "balanced" set: report estimate AND balance quality (max |std. gap|) at each cutoff
  sweep<-map_dfr(c(0.25,0.33,0.40,0.50),function(cut){ids<-rl$boundary_id[rl$gnorm<=quantile(rl$gnorm,cut)];r<-reg(ids)
    bb<-rl%>%filter(boundary_id%in%ids)
    tibble(seed=seed,cut=cut,est=r["est"],se=r["se"],nb=r["nb"],npoi=r["npoi"],
           black=mean(bb$black_1930,na.rm=T),maxabsz=max(abs(colMeans(as.matrix(bb[,zc])))))})
  bal<-rl%>%mutate(strat=case_when(boundary_id%in%bal_ids~"balanced",boundary_id%in%imb_ids~"imbalanced",TRUE~"mid"))%>%
    filter(strat!="mid")%>%select(strat,all_of(keycov))%>%
    pivot_longer(-strat,names_to="cov")%>%group_by(cov,strat)%>%summarise(gap=mean(value,na.rm=T),.groups="drop")%>%
    bind_rows(rl%>%select(all_of(keycov))%>%pivot_longer(everything(),names_to="cov")%>%group_by(cov)%>%summarise(strat="all",gap=mean(value),.groups="drop"))%>%mutate(seed=seed)
  list(est=tibble(seed=seed,est_all=a["est"],se_all=a["se"],nb_all=a["nb"],
                  est_bal=b["est"],se_bal=b["se"],nb_bal=b["nb"],npoi_bal=b["npoi"],
                  est_imb=im["est"],se_imb=im["se"],nb_imb=im["nb"]),bal=bal,sweep=sweep)}
RR<-map(c(1,3,13),run); R<-map_dfr(RR,"est"); B<-map_dfr(RR,"bal"); SW<-map_dfr(RR,"sweep")
cat("=== Covariate balance: pre-map (1930) cross-border gap (lower - higher side) ===\n")
B%>%group_by(cov,strat)%>%summarise(gap=mean(gap),.groups="drop")%>%pivot_wider(names_from=strat,values_from=gap)%>%
  mutate(across(c(all,balanced,imbalanced),~round(.x,4)))%>%select(cov,all,balanced,imbalanced)%>%as.data.frame()%>%print(row.names=FALSE)
cat("  (balanced half: gaps ~0 -> grading as-good-as-random; imbalanced half: large pre-map gaps)\n\n")
s<-R%>%summarise(across(-seed,mean)); pf<-function(e,se)sprintf("beta=%+.5f (se %.5f, t=%.2f, p=%.3f)",e,se,e/se,2*pnorm(-abs(e/se)))
cat("=== POI-level simple border comparison (seg ~ lower-graded side | boundary + category FE) ===\n")
cat(sprintf("  ALL real borders:        %s | %.0f boundaries\n",pf(s$est_all,s$se_all),s$nb_all))
cat(sprintf("  BALANCED half (clean):  %s | %.0f boundaries, %.0f POIs\n",pf(s$est_bal,s$se_bal),s$nb_bal,s$npoi_bal))
cat(sprintf("  IMBALANCED half:        %s | %.0f boundaries\n",pf(s$est_imb,s$se_imb),s$nb_imb))
cat("\n=== balanced-set cutoff sweep (transparency: estimate AND balance quality at each cutoff) ===\n")
SW%>%group_by(cut)%>%summarise(across(c(est,se,nb,black,maxabsz),mean),.groups="drop")%>%
  mutate(t=est/se,p=2*pnorm(-abs(t)),across(c(est,se,black),~round(.x,5)),maxabsz=round(maxabsz,3),t=round(t,2),p=round(p,3),nb=round(nb))%>%
  select(cut,nb,black_gap=black,max_abs_std_gap=maxabsz,est,se,t,p)%>%as.data.frame()%>%print(row.names=FALSE)
cat("  (max_abs_std_gap ~ balance quality: smaller = better balanced. Choose the largest cutoff that stays balanced.)\n")
saveRDS(list(est=R,bal=B,sweep=SW),"results/holc_lowprop_poi_FRAC.rds")
