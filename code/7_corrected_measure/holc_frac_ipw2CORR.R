# FRACTIONAL measure campaign, part 2: everything downstream of the headline DiDs.
# Outcome = seg_frac (proportional/fractional attribution, seg_eco*5/8 from ei_bounds_seg_frac.R).
# Reproduces, on the fractional measure: raw boundary discontinuity (Table 1 Col 2), probit-IPW
# DiD (Col 3, 100 draws), overall + by-group compositional mediation (Fig 4; seeds 1/3/13),
# amenity-mix decomposition (seeds 1/3/13), 7-type category heterogeneity and 6-bin catchment
# heterogeneity (Fig 5; seeds 1/3/13). Design identical to the holc_cdb_* scripts.
suppressMessages({library(tidyverse); library(haven); library(fixest)})
ahm_path<-"./data/ahm_staging/"; a2f<-function(x)paste0(substr(x,2,3),substr(x,5,7),substr(x,9,14))
BAR<-c("BUFFintrain","BUFFinsrivers","BUFFinbrivers")
covset<-c("black_","ownhome_","foreign_born_","read_write_","house_value_","rent_","radio_"); yrs<-c(1910,1920,1930)
load("./data/safegraph_2019_m.rdata")
places_usa_2019$other_rate<-pmax(0,1-(places_usa_2019$black_rate+places_usa_2019$white_rate+places_usa_2019$hispanic_rate+places_usa_2019$asian_rate))
grps<-c(black="black_rate",white="white_rate",hispanic="hispanic_rate",asian="asian_rate",other="other_rate")
msa_comp<-places_usa_2019%>%distinct(poi_cbg,.keep_all=TRUE)%>%group_by(msa)%>%
  summarise(across(all_of(unname(grps)),~weighted.mean(.x,tot_pop,na.rm=TRUE),.names="msa_{.col}"),.groups="drop")
catgrp<-function(tc){tc<-as.character(tc);case_when(
  grepl("Religious",tc)~"Religious",
  grepl("Grocery|Supermarket|Convenience|Specialty Food|Health and Personal Care|Pharmac|Drug|Clothing|Department|Merchandise|Shoe|Book|Hardware|Furniture|Florist|Jewelry",tc)~"Retail & grocery",
  grepl("Restaurant|Eating|Drinking|Snack|Caterers",tc)~"Food & drink",
  grepl("Personal Care Services|Beauty|Barber|Nail|Hair|Spa\b",tc)~"Personal care",
  grepl("Amusement|Recreation|Museum|Historical|Fitness|Golf|Bowling|Arts|Performing|Gambling|Zoo|Nature Park",tc)~"Recreation & culture",
  grepl("School|Day Care|Colleges|Universit|Educational|Instruction",tc)~"Education & childcare",
  grepl("Physician|Dentist|Health Practitioner|Medical|Hospital|Nursing|Outpatient|Mental|Chiropract|Optomet|Diagnostic",tc)~"Health care",
  TRUE~"Other")}
fr<-as_tibble(readRDS("results/poi_ei_bounds_fraccorr.rds"))%>%transmute(safegraph_place_id,seg_frac=seg_eco*5/8)
poi<-places_usa_2019%>%left_join(msa_comp,by="msa")%>%inner_join(fr,by="safegraph_place_id")%>%
  mutate(lvis=log(raw_visit_counts+1),tract11=substr(poi_cbg,1,11),
         catch=distance_from_home,catgroup=catgrp(top_category))
for(g in names(grps)){nb<-grps[[g]];ms<-paste0("msa_",nb)
  poi[[paste0("dev_",g)]]<-poi[[nb]]-poi[[ms]]; poi[[paste0("absdev_",g)]]<-abs(poi[[nb]]-poi[[ms]])}
poi$absdev_total<-rowSums(poi[,paste0("absdev_",names(grps))],na.rm=TRUE)
catm<-poi%>%filter(!is.na(seg_frac),!is.na(lvis),!is.na(sub_category))%>%group_by(sub_category)%>%
  summarise(m_cat=weighted.mean(seg_frac,lvis),.groups="drop")
poi<-poi%>%filter(!is.na(seg_frac),!is.na(lvis),!is.na(sub_category))%>%
  left_join(catm,by="sub_category")%>%mutate(res_cat=seg_frac-m_cat)%>%
  select(tract11,sub_category,catgroup,catch,seg_frac,m_cat,res_cat,lvis,starts_with("dev_"),starts_with("absdev_"))
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
mkdat<-function(seed,ipw=FALSE){side<-bt%>%group_by(boundary_id,real,boundary_side)%>%summarise(grade=first(holc_grade[holc_grade%in%c("A","B","C","D")]),across(all_of(cvars),~mean(.x,na.rm=T)),.groups="drop")%>%group_by(boundary_id)%>%filter(n()==2)%>%ungroup()
  set.seed(seed);side<-side%>%group_by(boundary_id)%>%mutate(lgs=if(first(real)==1L)as.integer(grade%in%c("C","D"))else{u<-runif(n());as.integer(u==max(u))})%>%ungroup()
  gap<-side%>%group_by(boundary_id,real)%>%summarise(across(all_of(cvars),~mean(.x[lgs==1],na.rm=T)-mean(.x[lgs==0],na.rm=T)),.groups="drop")
  bvars<-cvars[sapply(cvars,function(v)sum(is.finite(as.numeric(scale(gap[[v]]))))>100)];gs<-gap[,"real"];cons<-c()
  for(v in bvars){z<-as.numeric(scale(gap[[v]]));gs[[paste0("Z_",v)]]<-ifelse(is.finite(z),z,0);cons<-c(cons,paste0("Z_",v))
    if(any(!is.finite(z))){gs[[paste0("MI_",v)]]<-as.numeric(!is.finite(z));cons<-c(cons,paste0("MI_",v))}}
  tgt<-colMeans(gs[gs$real==1,cons,drop=FALSE]);wp<-ebalance(as.matrix(gs[gs$real==0,cons,drop=FALSE]),tgt)*sum(gs$real==0)
  gg<-gap%>%mutate(w_att=ifelse(real==1L,1,NA));gg$w_att[gg$real==0]<-wp
  if(ipw){ps<-suppressWarnings(glm(reformulate(cons,"real"),data=cbind(gs,real=gap$real),family=binomial("probit")))
    ph<-predict(ps,type="response");odds<-ph/pmax(1-ph,.02)
    gg$w_ipw<-ifelse(gap$real==1,1,odds/mean(odds[gap$real==0]))}
  side%>%select(boundary_id,boundary_side,lgs,real)%>%inner_join(gg%>%select(boundary_id,any_of(c("w_att","w_ipw"))),by="boundary_id")%>%
    inner_join(bt%>%distinct(boundary_id,boundary_side,tract11),by=c("boundary_id","boundary_side"))%>%
    inner_join(poi,by="tract11",relationship="many-to-many")%>%mutate(w=w_att*lvis)}
OUT<-list()
cat("\n=== (f) IPW DiD, 100 draws (Col 3 analog) ===\n")
R<-map_dfr(1:100,function(s){d<-mkdat(s,ipw=TRUE)
  m<-feols(seg_frac~lgs+lgs:real,d,fixef=c("boundary_id","sub_category"),weights=~I(w_ipw*lvis),cluster=~boundary_id)
  tibble(seed=s,plac=unname(coef(m)["lgs:real"])*0+unname(coef(m)["lgs"]),did=unname(coef(m)["lgs:real"]),se=unname(se(m)["lgs:real"]),p=unname(pvalue(m)["lgs:real"]))})
OUT$ipw<-R
cat(sprintf("IPW DiD %+.5f | rand CI [%+.5f,%+.5f] | %%p<.05 %.0f\n",mean(R$did),quantile(R$did,.025),quantile(R$did,.975),100*mean(R$p<.05)))
saveRDS(OUT,"results/holc_frac_ipw2_CORR.rds")
cat("Saved results/holc_frac_ipw2_CORR.rds\n")
