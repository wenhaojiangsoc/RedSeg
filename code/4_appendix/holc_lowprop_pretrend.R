# PRE-TREND for the SECOND identification (balanced / low-propensity boundaries, app:lowprop).
# AHM (2021) show their low-propensity borders have "no pre-trend -- the gap in 1910, 1920 and 1930 is
# essentially zero." We reproduce the balanced/imbalanced split of holc_lowprop.R and report, for each
# stratum, the mean within-boundary gap (lower- minus higher-graded side) of the key covariates at
# 1910, 1920, 1930 -- i.e., is the pre-map TRAJECTORY flat/zero on the balanced boundaries?
suppressMessages({library(tidyverse); library(haven)})
ahm_path<-"./data/ahm_staging/"; a2f<-function(x)paste0(substr(x,2,3),substr(x,5,7),substr(x,9,14))
BAR<-c("BUFFintrain","BUFFinsrivers","BUFFinbrivers")
covset<-c("black_","ownhome_","foreign_born_","read_write_","house_value_","rent_","radio_"); yrs<-c(1910,1920,1930)
read_map<-function(file,idcol,pfx,is_grid){d<-read_dta(paste0(ahm_path,file))%>%filter(year==2010)%>%mutate(tract11=a2f(trctID2010))
  if(is_grid) d<-d%>%rename_with(~sub("BUFF_CFin","BUFFin",.x))
  d%>%select(tract11,matches(paste0("^",idcol,"[0-9]+$")),matches("^boundary_side[0-9]+$|^side_grid[0-9]+$"),matches("^holc_grade[0-9]+$"),
      matches("^boundary_grades[0-9]+$"),matches("^BUFFintrain[0-9]+$|^BUFFinsrivers[0-9]+$|^BUFFinbrivers[0-9]+$"))%>%
    rename_with(~sub("^side_grid","boundary_side",.x))%>%rename_with(~sub(paste0("^",idcol),"boundary_id",.x))%>%
    pivot_longer(-tract11,names_to=c(".value","slot"),names_pattern="^(boundary_id|boundary_side|holc_grade|boundary_grades|BUFFintrain|BUFFinsrivers|BUFFinbrivers)([0-9]+)$")%>%
    drop_na(boundary_id)%>%mutate(boundary_id=paste0(pfx,boundary_id))}
mA<-read_map("redlining_dataset_t2010_4th_new.dta","boundary_id","A_",FALSE)
raw<-read_dta(paste0(ahm_path,"redlining_dataset_t2010_4th_new.dta"))%>%filter(year%in%yrs)%>%mutate(tract11=a2f(trctID2010))
cov_all<-raw%>%select(tract11,year,any_of(covset))%>%group_by(tract11,year)%>%summarise(across(everything(),~mean(.x,na.rm=T)),.groups="drop")%>%pivot_wider(names_from=year,values_from=any_of(covset),names_glue="{.value}{year}");cvars<-setdiff(names(cov_all),"tract11")
tr<-mA%>%filter(boundary_grades%in%c("2_to_3","2_to_4"))%>%mutate(real=1L)
bt<-tr%>%left_join(cov_all,by="tract11")%>%mutate(barrier_amt=rowSums(across(all_of(BAR),~replace_na(.x,0))))
bt<-bt%>%semi_join(bt%>%group_by(boundary_id)%>%summarise(b=as.integer(max(barrier_amt,na.rm=T)>0.05))%>%filter(b==0),by="boundary_id")
yv<-intersect(paste0(rep(c("black_","ownhome_","foreign_born_"),each=3),c(1910,1920,1930)),cvars)   # trajectory covariates
keycov<-intersect(paste0(c("black_","ownhome_","house_value_","foreign_born_"),"1930"),cvars)        # selection covariates (1930)
run<-function(seed){
  side<-bt%>%group_by(boundary_id,boundary_side)%>%summarise(grade=first(holc_grade[holc_grade%in%c("A","B","C","D")]),across(all_of(cvars),~mean(.x,na.rm=T)),.groups="drop")%>%
    group_by(boundary_id)%>%filter(n()==2)%>%mutate(lgs=as.integer(grade%in%c("C","D")))%>%ungroup()
  gap<-side%>%group_by(boundary_id)%>%summarise(across(all_of(cvars),~mean(.x[lgs==1],na.rm=T)-mean(.x[lgs==0],na.rm=T)),.groups="drop")
  rl<-gap%>%filter(if_all(all_of(keycov),is.finite))
  for(v in keycov) rl[[paste0("z_",v)]]<-as.numeric(scale(rl[[v]]))
  rl$gnorm<-sqrt(rowSums(as.matrix(rl[,paste0("z_",keycov)])^2))
  bal<-rl$boundary_id[rl$gnorm<=quantile(rl$gnorm,.5)]; imb<-rl$boundary_id[rl$gnorm>quantile(rl$gnorm,.5)]
  tab<-function(ids,lab) gap%>%filter(boundary_id%in%ids)%>%summarise(across(all_of(yv),~mean(.x,na.rm=T)))%>%mutate(strat=lab)
  bind_rows(tab(gap$boundary_id,"all"),tab(bal,"balanced"),tab(imb,"imbalanced"))%>%mutate(seed=seed)}
R<-map_dfr(c(1,3,13),run)%>%group_by(strat)%>%summarise(across(all_of(yv),mean),.groups="drop")
cat("Within-boundary gap (lower - higher side) by pre-map year, for the SECOND-identification strata:\n")
long<-R%>%pivot_longer(-strat,names_to="cv",values_to="gap")%>%separate(cv,c("cov","year"),sep="(?<=_)(?=[0-9])")%>%mutate(cov=sub("_$","",cov))
for(c0 in c("black","ownhome","foreign_born")){
  cat(sprintf("\n%-14s  1910     1920     1930\n",c0))
  for(s in c("balanced","imbalanced","all")){
    v<-long%>%filter(cov==c0,strat==s)%>%arrange(year)%>%pull(gap)
    cat(sprintf("  %-11s %+.4f  %+.4f  %+.4f\n",s,v[1],v[2],v[3]))}}
