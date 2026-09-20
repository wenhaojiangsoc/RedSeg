# NON-LOCAL robustness (step 2 of 2): entropy-balanced boundary DiD on the recomputed segregation
# measures (all visitors / excl own CBG / excl own tract) from holc_nonlocal_seg.R. Tests whether the
# redlining effect survives when local residents are excluded. Also runs the composition mediation on
# the non-local measure. Same design/weights/clustering as the main estimator.
suppressMessages({library(tidyverse); library(haven); library(fixest)})
ahm_path<-"./data/ahm_staging/"; a2f<-function(x)paste0(substr(x,2,3),substr(x,5,7),substr(x,9,14))
BAR<-c("BUFFintrain","BUFFinsrivers","BUFFinbrivers")
covset<-c("black_","ownhome_","foreign_born_","read_write_","house_value_","rent_","radio_"); yrs<-c(1910,1920,1930)
S <- as_tibble(readRDS("results/poi_nonlocal_seg_FRACCORR.rds"))   # safegraph_place_id,tract11,seg_all,seg_excbg,seg_extract,share_cbg,share_tract,t
# join POI category + neighborhood Black-vs-MSA distance for the mediation
load("./data/safegraph_2019_m.rdata")
msa_black<-places_usa_2019%>%distinct(poi_cbg,.keep_all=TRUE)%>%group_by(msa)%>%summarise(mb=weighted.mean(black_rate,tot_pop,na.rm=TRUE),.groups="drop")
meta<-places_usa_2019%>%left_join(msa_black,by="msa")%>%transmute(safegraph_place_id,sub_category,absdev=abs(black_rate-mb))
poi<-S%>%left_join(meta,by="safegraph_place_id")%>%mutate(lvis=log(t+1))%>%
  filter(!is.na(seg_all),!is.na(lvis),!is.na(sub_category))%>%
  select(tract11,sub_category,seg_all,seg_excbg,seg_extract,absdev,lvis)
read_map<-function(file,idcol,pfx,is_grid){d<-read_dta(paste0(ahm_path,file))%>%filter(year==2010)%>%mutate(tract11=a2f(trctID2010))
  if(is_grid) d<-d%>%rename_with(~sub("BUFF_CFin","BUFFin",.x))
  d%>%select(tract11,matches(paste0("^",idcol,"[0-9]+$")),matches("^boundary_side[0-9]+$|^side_grid[0-9]+$"),matches("^holc_grade[0-9]+$"),
      matches("^boundary_grades[0-9]+$"),matches("^BUFFintrain[0-9]+$"),matches("^BUFFinsrivers[0-9]+$"),matches("^BUFFinbrivers[0-9]+$"))%>%
    rename_with(~sub("^side_grid","boundary_side",.x))%>%rename_with(~sub(paste0("^",idcol),"boundary_id",.x))%>%
    pivot_longer(-tract11,names_to=c(".value","slot"),names_pattern="^(boundary_id|boundary_side|holc_grade|boundary_grades|BUFFintrain|BUFFinsrivers|BUFFinbrivers)([0-9]+)$")%>%
    drop_na(boundary_id)%>%mutate(boundary_id=paste0(pfx,boundary_id))}
mA<-read_map("redlining_dataset_t2010_4th_new.dta","boundary_id","A_",FALSE);mG<-read_map("redlining_dataset_tgrid2010_4th_new.dta","boundary_id_grid","G_",TRUE)
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
  bvars<-cvars[sapply(cvars,function(v)sum(is.finite(as.numeric(scale(gap[[v]]))))>100)]
  gs<-gap[,"real"];cons<-c()
  for(v in bvars){z<-as.numeric(scale(gap[[v]]));gs[[paste0("Z_",v)]]<-ifelse(is.finite(z),z,0);cons<-c(cons,paste0("Z_",v));if(any(!is.finite(z))){gs[[paste0("MI_",v)]]<-as.numeric(!is.finite(z));cons<-c(cons,paste0("MI_",v))}}
  tgt<-colMeans(gs[gs$real==1,cons,drop=FALSE]);wp<-ebalance(as.matrix(gs[gs$real==0,cons,drop=FALSE]),tgt)*sum(gs$real==0)
  gg<-gap%>%mutate(w_att=ifelse(real==1L,1,NA));gg$w_att[gg$real==0]<-wp
  side%>%select(boundary_id,boundary_side,lgs,real)%>%inner_join(gg%>%select(boundary_id,w_att),by="boundary_id")%>%
    inner_join(bt%>%distinct(boundary_id,boundary_side,tract11),by=c("boundary_id","boundary_side"))%>%inner_join(poi,by="tract11",relationship="many-to-many")%>%mutate(w=w_att*lvis)}
g<-function(m)c(b=unname(coef(m)["lgs:real"]),s=unname(se(m)["lgs:real"]))
run<-function(seed){d<-mkdat(seed)
  f<-function(y)g(feols(as.formula(paste0(y,"~lgs+lgs:real")),d,fixef=c("boundary_id","sub_category"),weights=~w,cluster=~boundary_id))
  # composition mediation on the non-local (excl-CBG) measure
  mc<-feols(seg_excbg~lgs+lgs:real|boundary_id+sub_category,d,weights=~w,cluster=~boundary_id)
  ma<-feols(absdev~lgs+lgs:real|boundary_id+sub_category,d,weights=~w,cluster=~boundary_id)
  mb<-feols(seg_excbg~lgs+lgs:real+absdev|boundary_id+sub_category,d,weights=~w,cluster=~boundary_id)
  cc<-g(mc);a<-g(ma);b<-c(coef(mb)["absdev"],se(mb)["absdev"]);ind<-a["b"]*b[1];prop<-ind/cc["b"]
  tibble(seed=seed,all=f("seg_all")["b"],all_s=f("seg_all")["s"],excbg=f("seg_excbg")["b"],excbg_s=f("seg_excbg")["s"],
         extract=f("seg_extract")["b"],extract_s=f("seg_extract")["s"],med_excbg=unname(prop))}
R<-map_dfr(c(1,3,7,13,21),run);s<-summarise(R,across(-seed,mean))
pf<-function(b,se)sprintf("%+.5f (se %.5f, t=%.2f, p=%.4f)",b,se,b/se,2*pnorm(-abs(b/se)))
cat("=== NON-LOCAL boundary DiD (avg 5 orientations) ===\n")
cat(sprintf("all visitors (main DV):   %s\n",pf(s$all,s$all_s)))
cat(sprintf("excl. own-CBG residents:  %s\n",pf(s$excbg,s$excbg_s)))
cat(sprintf("excl. own-tract residents:%s\n",pf(s$extract,s$extract_s)))
cat(sprintf("composition mediation of excl-CBG effect: %.0f%%\n",100*s$med_excbg))
