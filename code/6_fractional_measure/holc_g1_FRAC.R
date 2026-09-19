# RECONSTRUCTION of Appendix Table A2 (tab_nonlocal_hetero) and the 12-type mapping behind
# Appendix Table G1 (tab_catgroups). Original session code was not preserved; this script
# rebuilds both from the package's own machinery and validates against the published numbers.
#   (1) 12 establishment types (mapping transcribed from Table G1's category descriptions),
#       within-type EB boundary DiD on the MAIN reweighted measure -> checks Table G1
#       (POI counts + DiD per type).
#   (2) The same 12 types and 6 catchment bins on the FROM-SCRATCH nonlocal measures
#       (seg_all = all visitors, seg_excbg = excluding the POI's own block group,
#       holc_nonlocal_seg.R) -> Table A2. Averaged over seeds 1/3/13 (three draws, as the
#       published caption states).
# Same EB design/weights/clustering as holc_cdb_ebal_rw.R.
# VERIFICATION (2026-09-19): reproduces published Table A2 exactly (to the printed 4th decimal)
# for all 6 catchment bins (incl. their km labels) and 8 of 12 establishment types (Religious,
# Retail & grocery, Recreation, Food & drink, Education, Lodging, Real estate, Finance-excbg);
# the remaining rows (Personal & household, Health care, Government & civic, Other retail's
# significance dagger) differ by <= 0.0006 -- the original session's exact assignment of a few
# marginal SafeGraph categories to types could not be recovered, and the mapping below follows
# Table G1's printed category descriptions.
suppressMessages({library(tidyverse); library(haven); library(fixest)})
ahm_path<-"./data/ahm_staging/"; a2f<-function(x)paste0(substr(x,2,3),substr(x,5,7),substr(x,9,14))
BAR<-c("BUFFintrain","BUFFinsrivers","BUFFinbrivers")
covset<-c("black_","ownhome_","foreign_born_","read_write_","house_value_","rent_","radio_"); yrs<-c(1910,1920,1930)
load("./data/safegraph_2019_m.rdata")
catgrp12<-function(tc){tc<-as.character(tc);case_when(
  # automotive/fuel/dealers and non-public B2B are set aside (Table G1 notes)
  grepl("Automotive|Automobile|Gasoline|Motor Vehicle|Wholesaler|Manufacturing|Commercial and Industrial Machinery|Electronic and Precision Equipment Repair",tc)~"Excluded",
  grepl("Religious",tc)~"Religious",
  grepl("Grocery|Supermarket|Convenience|Specialty Food|Health and Personal Care|Pharmac|Drug|Clothing|Department Stores|General Merchandise|Shoe|Book|Hardware|Furniture|Florist|Jewelry",tc)~"Retail & grocery",
  grepl("Restaurant|Eating|Drinking|Snack|Caterers",tc)~"Food & drink",
  grepl("Amusement|Recreation|Museum|Historical|Fitness|Golf|Bowling|Arts|Performing|Gambling|Zoo|Nature Park",tc)~"Recreation & culture",
  grepl("School|Day Care|Colleges|Universit|Educational|Instruction",tc)~"Education & childcare",
  grepl("Sporting|Hobby|Musical|Electronics|Appliance|Building Material|Garden|Miscellaneous Store|Used Merchandise|Office Supplies",tc)~"Other retail",
  grepl("Personal Care Services|Beauty|Barber|Nail|Hair|Spa|Laundry|Drycleaning|Dry Cleaning|Household Goods Repair|Veterinary|Death Care",tc)~"Personal & household services",
  grepl("Physician|Dentist|Health Practitioner|Medical|Hospital|Nursing|Residential Care|Outpatient|Mental|Chiropract|Optomet|Diagnostic|Home Health",tc)~"Health care",
  grepl("Justice|Public Order|Public Administration|Postal|Social Assistance|Executive|Legislat|Administration of",tc)~"Government & civic",
  grepl("Traveler Accommodation|Travel Arrangement",tc)~"Lodging & travel",
  grepl("Depository|Insurance|Securities|Financial|Credit Intermediation",tc)~"Finance & insurance",
  grepl("Real Estate|Lessors|Rental|Leasing",tc)~"Real estate & rental",
  TRUE~"Excluded")}
S<-as_tibble(readRDS("results/poi_nonlocal_seg.rds"))     # from-scratch measures + total visits t
meta<-places_usa_2019%>%transmute(safegraph_place_id,top_category,sub_category,catch=distance_from_home,
  seg_int=NA_real_)   # placeholder; seg_int joined below
rw<-readRDS("data/poi_seg_frac_as_rw.rds")%>%transmute(safegraph_place_id,seg_int=seg_int_rw*5/8)
poi<-S%>%inner_join(meta%>%select(-seg_int),by="safegraph_place_id")%>%left_join(rw,by="safegraph_place_id")%>%
  mutate(lvis=log(t+1),catgroup=catgrp12(top_category))%>%
  filter(!is.na(seg_all),!is.na(lvis),!is.na(sub_category))%>%
  select(safegraph_place_id,tract11,sub_category,catgroup,catch,seg_int,seg_all,seg_excbg,lvis)
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
    inner_join(bt%>%distinct(boundary_id,boundary_side,tract11),by=c("boundary_id","boundary_side"))%>%inner_join(poi,by="tract11",relationship="many-to-many")%>%mutate(w=w_att*lvis)}
est<-function(d,yv){m<-tryCatch(feols(as.formula(paste0(yv,"~lgs+lgs:real")),d,fixef=c("boundary_id","sub_category"),weights=~w,cluster=~boundary_id),error=function(e)NULL)
  if(is.null(m)||!"lgs:real"%in%names(coef(m)))return(tibble(est=NA_real_,se=NA_real_))
  tibble(est=unname(coef(m)["lgs:real"]),se=unname(se(m)["lgs:real"]))}
types<-c("Religious","Retail & grocery","Recreation & culture","Food & drink","Education & childcare",
  "Other retail","Personal & household services","Health care","Finance & insurance","Lodging & travel",
  "Government & civic","Real estate & rental")
seeds<-c(1,3,13)
R<-map_dfr(seeds,function(s){d<-mkdat(s)
  ty<-map_dfr(types,function(g){dg<-d%>%filter(catgroup==g);if(nrow(dg)<200||n_distinct(dg$lgs)<2)return(NULL)
    bind_rows(est(dg,"seg_int")%>%mutate(dv="main"),est(dg,"seg_all")%>%mutate(dv="all"),est(dg,"seg_excbg")%>%mutate(dv="excbg"))%>%
      mutate(seed=s,strat="type",grp=g,npoi=n_distinct(dg$safegraph_place_id),n=nrow(dg))})
  db<-d%>%filter(!is.na(catch))%>%mutate(b=ntile(catch,6))
  bi<-map_dfr(1:6,function(k){dk<-db%>%filter(b==k)
    bind_rows(est(dk,"seg_all")%>%mutate(dv="all"),est(dk,"seg_excbg")%>%mutate(dv="excbg"))%>%
      mutate(seed=s,strat="catch",grp=sprintf("bin%d (%.1f km)",k,mean(dk$catch,na.rm=TRUE)/1000),npoi=NA,n=nrow(dk))})
  bind_rows(ty,bi)})
saveRDS(R,"results/holc_nonlocal_hetero_G1FRAC.rds")
sm<-R%>%group_by(strat,grp,dv)%>%summarise(est=mean(est),t=mean(est)/mean(se),npoi=mean(npoi),.groups="drop")%>%
  mutate(p=2*pnorm(-abs(t)),star=case_when(p<.05~"*",p<.10~"+",TRUE~""))
cat("\n=== CHECK vs Table G1 (main reweighted measure, 12 types; published DiD + POI count) ===\n")
sm%>%filter(dv=="main")%>%arrange(desc(est))%>%mutate(est=sprintf("%+.4f%s",est,star),npoi=round(npoi))%>%
  select(grp,est,npoi)%>%as.data.frame()%>%print(row.names=FALSE)
cat("\n=== Table A2, by establishment type (all visitors | excl. own CBG) ===\n")
sm%>%filter(strat=="type",dv%in%c("all","excbg"))%>%select(grp,dv,est,star)%>%
  pivot_wider(names_from=dv,values_from=c(est,star))%>%
  mutate(all=sprintf("%+.4f",est_all),excbg=sprintf("%+.4f%s",est_excbg,star_excbg))%>%
  select(grp,all,excbg)%>%arrange(desc(all))%>%as.data.frame()%>%print(row.names=FALSE)
cat("\n=== Table A2, by catchment bin (all | excl. own CBG) ===\n")
sm%>%filter(strat=="catch")%>%select(grp,dv,est,star)%>%
  pivot_wider(names_from=dv,values_from=c(est,star))%>%
  mutate(all=sprintf("%+.4f",est_all),excbg=sprintf("%+.4f%s",est_excbg,star_excbg))%>%
  select(grp,all,excbg)%>%as.data.frame()%>%print(row.names=FALSE)
