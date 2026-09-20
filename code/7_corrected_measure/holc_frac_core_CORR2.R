# FRACTIONAL (proportional) attribution as the MAIN outcome: core causal package.
# The published main measure attributed each home CBG's population to its LARGEST race group
# (data/cbg_race_2019.csv is one-hot); Eq. 1 / Figure 1 describe the PROPORTIONAL blend.
# This script re-estimates the two headline boundary DiDs on the proportional measure
# (seg_eco from ei_bounds_seg_frac.R, ACS B03002 fractional composition), 100 draws each:
#   (1) tract-buffer entropy-balanced DiD (Table 1 Col 4 analog) + within-boundary SD;
#   (2) POI-level assignment DiD (Table 1 Col 5 analog) + within-boundary SD.
# Same design, weights, FE, clustering, and barrier exclusions as holc_cdb_ebal_rw.R /
# holc_cdb_poi.R.
suppressMessages({library(tidyverse); library(haven); library(fixest)})
ahm_path<-"./data/ahm_staging/"; a2f<-function(x)paste0(substr(x,2,3),substr(x,5,7),substr(x,9,14))
BAR<-c("BUFFintrain","BUFFinsrivers","BUFFinbrivers")
covset<-c("black_","ownhome_","foreign_born_","read_write_","house_value_","rent_","radio_"); yrs<-c(1910,1920,1930)
GRID_M<-804.672; BUFFER_M<-402.336
load("./data/safegraph_2019_m.rdata")
fr<-as_tibble(readRDS("results/poi_ei_bounds_fraccorr.rds"))%>%transmute(safegraph_place_id,seg_frac=seg_eco*5/8)
poi<-places_usa_2019%>%inner_join(fr,by="safegraph_place_id")%>%
  mutate(lvis=log(raw_visit_counts+1),tract11=substr(poi_cbg,1,11))%>%
  filter(!is.na(seg_frac),!is.na(lvis),!is.na(sub_category))%>%
  select(safegraph_place_id,tract11,sub_category,seg_frac,lvis,latitude,longitude)
cat(sprintf("POIs with fractional measure: %d | mean seg_frac %.4f\n",nrow(poi),mean(poi$seg_frac)))
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
mkdat<-function(seed){side<-bt%>%group_by(boundary_id,real,boundary_side)%>%summarise(grade=first(holc_grade[holc_grade%in%c("A","B","C","D")]),across(all_of(cvars),~mean(.x,na.rm=T)),.groups="drop")%>%group_by(boundary_id)%>%filter(n()==2)%>%ungroup()
  set.seed(seed);side<-side%>%group_by(boundary_id)%>%mutate(lgs=if(first(real)==1L)as.integer(grade%in%c("C","D"))else{u<-runif(n());as.integer(u==max(u))})%>%ungroup()
  gap<-side%>%group_by(boundary_id,real)%>%summarise(across(all_of(cvars),~mean(.x[lgs==1],na.rm=T)-mean(.x[lgs==0],na.rm=T)),.groups="drop")
  bvars<-cvars[sapply(cvars,function(v)sum(is.finite(as.numeric(scale(gap[[v]]))))>100)];gs<-gap[,"real"];cons<-c()
  for(v in bvars){z<-as.numeric(scale(gap[[v]]));gs[[paste0("Z_",v)]]<-ifelse(is.finite(z),z,0);cons<-c(cons,paste0("Z_",v))
    if(any(!is.finite(z))){gs[[paste0("MI_",v)]]<-as.numeric(!is.finite(z));cons<-c(cons,paste0("MI_",v))}}
  tgt<-colMeans(gs[gs$real==1,cons,drop=FALSE]);wp<-ebalance(as.matrix(gs[gs$real==0,cons,drop=FALSE]),tgt)*sum(gs$real==0)
  gg<-gap%>%mutate(w_att=ifelse(real==1L,1,NA));gg$w_att[gg$real==0]<-wp
  side%>%select(boundary_id,boundary_side,lgs,real)%>%inner_join(gg%>%select(boundary_id,w_att),by="boundary_id")%>%
    inner_join(bt%>%distinct(boundary_id,boundary_side,tract11),by=c("boundary_id","boundary_side"))%>%
    inner_join(poi%>%select(tract11,sub_category,seg_frac,lvis),by="tract11",relationship="many-to-many")%>%mutate(w=w_att*lvis)}
cat("\n=== (1) TRACT-BUFFER EB DiD on fractional measure, 100 draws ===\n")
d1<-mkdat(1)
sd_tract<-sd(resid(feols(seg_frac~1,d1,fixef=c("boundary_id","sub_category"))))
R1<-map_dfr(1:100,function(s){d<-if(s==1)d1 else mkdat(s)
  m<-feols(seg_frac~lgs+lgs:real,d,fixef=c("boundary_id","sub_category"),weights=~w,cluster=~boundary_id)
  tibble(seed=s,did=unname(coef(m)["lgs:real"]),se=unname(se(m)["lgs:real"]),p=unname(pvalue(m)["lgs:real"]),plac=unname(coef(m)["lgs"]))})
saveRDS(list(draws=R1,sd_within=sd_tract),"results/holc_frac_tract_CORR2.rds")
s1<-R1%>%summarise(did=mean(did),se=mean(se),lo=quantile(did,.025),hi=quantile(did,.975),p05=mean(p<.05)*100)
cat(sprintf("DiD %+.5f (asy se %.5f) | rand CI [%+.5f,%+.5f] | %%draws p<.05: %.0f\n",s1$did,s1$se,s1$lo,s1$hi,s1$p05))
cat(sprintf("within-boundary SD (unweighted, net FE): %.4f -> effect = %.1f%% of SD\n",sd_tract,100*abs(s1$did)/sd_tract))
rm(d1);gc(FALSE)
cat("\n=== (2) POI-LEVEL DiD on fractional measure, 100 draws ===\n")
real<-readRDS("data/poi_real_assign.rds")%>%
  transmute(safegraph_place_id,boundary_id=paste0("R",boundary_id),real=1L,lgs=as.integer(side=="low"),gside=NA_real_,tract11,sub_category,lvis)%>%
  inner_join(fr,by="safegraph_place_id")
pts<-readRDS("data/poi_bcd_pts.rds"); G<-GRID_M
kx<-round(pts$X/G)*G; dxv<-pts$X-kx; ky<-round(pts$Y/G)*G; dyh<-pts$Y-ky; vert<-abs(dxv)<=abs(dyh)
plac<-pts%>%mutate(boundary_id=paste0("P",ifelse(vert,"V","H"),"_",ifelse(vert,round(X/G),round(Y/G)),"_",ifelse(vert,floor(Y/G),floor(X/G))),
    gside=ifelse(vert,sign(dxv),sign(dyh)),gdist=ifelse(vert,abs(dxv),abs(dyh)))%>%filter(gdist<=BUFFER_M,gside!=0)
valid<-plac%>%group_by(boundary_id)%>%summarise(ns=n_distinct(gside),ng=n_distinct(grade),.groups="drop")%>%filter(ns==2,ng==1)%>%pull(boundary_id)
plac<-plac%>%filter(boundary_id%in%valid)%>%inner_join(fr,by="safegraph_place_id")%>%
  transmute(safegraph_place_id,boundary_id,real=0L,lgs=NA_integer_,gside,tract11,sub_category,lvis,seg_frac)
allpoi<-bind_rows(real%>%select(-any_of("grade")),plac)%>%left_join(cov_all,by="tract11")
plac_bids<-unique(plac$boundary_id)
cat(sprintf("Real: %d POIs / %d borders | Placebo: %d POIs / %d borders\n",nrow(real),n_distinct(real$boundary_id),nrow(plac),n_distinct(plac$boundary_id)))
drawp<-function(seed){set.seed(seed)
  po<-tibble(boundary_id=plac_bids,low_side=sample(c(-1,1),length(plac_bids),replace=TRUE))
  d<-allpoi%>%left_join(po,by="boundary_id")%>%mutate(lgs=ifelse(real==1L,lgs,as.integer(gside==low_side)))
  bs<-d%>%group_by(boundary_id,real)%>%summarise(n0=sum(lgs==0),n1=sum(lgs==1),across(all_of(cvars),~mean(.x[lgs==1],na.rm=T)-mean(.x[lgs==0],na.rm=T)),.groups="drop")%>%filter(n0>0,n1>0)
  bvars<-cvars[sapply(cvars,function(v)sum(is.finite(as.numeric(scale(bs[[v]]))))>100)]
  gs<-bs[,c("boundary_id","real")];cons<-c()
  for(v in bvars){z<-as.numeric(scale(bs[[v]]));gs[[paste0("Z_",v)]]<-ifelse(is.finite(z),z,0);cons<-c(cons,paste0("Z_",v))
    if(any(!is.finite(z))){gs[[paste0("MI_",v)]]<-as.numeric(!is.finite(z));cons<-c(cons,paste0("MI_",v))}}
  tgt<-colMeans(gs[gs$real==1,cons,drop=FALSE]);wp<-ebalance(as.matrix(gs[gs$real==0,cons,drop=FALSE]),tgt)*sum(gs$real==0)
  bw<-gs%>%transmute(boundary_id,w_eb=ifelse(real==1,1,NA_real_));bw$w_eb[gs$real==0]<-wp
  dd<-d%>%inner_join(bw,by="boundary_id")
  m<-feols(seg_frac~lgs+lgs:real,dd,fixef=c("boundary_id","sub_category"),weights=~I(w_eb*lvis),cluster=~boundary_id)
  tibble(seed=seed,did=unname(coef(m)["lgs:real"]),se=unname(se(m)["lgs:real"]),p=unname(pvalue(m)["lgs:real"]),plac=unname(coef(m)["lgs"]))}
set.seed(1);po1<-tibble(boundary_id=plac_bids,low_side=sample(c(-1,1),length(plac_bids),replace=TRUE))
dsd<-allpoi%>%left_join(po1,by="boundary_id")%>%mutate(lgs=ifelse(real==1L,lgs,as.integer(gside==low_side)))
sd_poi<-sd(resid(feols(seg_frac~1,dsd,fixef=c("boundary_id","sub_category"))));rm(dsd);gc(FALSE)
R2<-map_dfr(1:100,drawp)
saveRDS(list(draws=R2,sd_within=sd_poi),"results/holc_frac_poi_CORR2.rds")
s2<-R2%>%summarise(did=mean(did),se=mean(se),lo=quantile(did,.025),hi=quantile(did,.975),p05=mean(p<.05)*100)
cat(sprintf("DiD %+.5f (asy se %.5f) | rand CI [%+.5f,%+.5f] | %%draws p<.05: %.0f\n",s2$did,s2$se,s2$lo,s2$hi,s2$p05))
cat(sprintf("within-boundary SD (unweighted, net FE): %.4f -> effect = %.1f%% of SD\n",sd_poi,100*abs(s2$did)/sd_poi))
