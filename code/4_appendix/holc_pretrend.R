# PRE-TREND check for the boundary design (analog of Aaronson-Hartley-Mazumder 2021, who show the
# reweighted comparison boundaries reproduce the TREATED boundaries' pre-map covariate trajectories).
# Mobility is observed only in 2019, so we cannot event-study the outcome; instead we ask whether the
# entropy-reweighted placebo boundaries reproduce the real boundaries' within-boundary covariate GAPS
# at each pre-map census (1910, 1920, 1930). Genuine test = LEAVE-OUT: balance on 1930 gaps only, then
# check whether the held-out 1910 and 1920 gaps also match (parallel pre-trends not mechanically imposed).
suppressMessages({library(tidyverse); library(haven)})
ahm_path<-"./data/ahm_staging/"; a2f<-function(x)paste0(substr(x,2,3),substr(x,5,7),substr(x,9,14))
BAR<-c("BUFFintrain","BUFFinsrivers","BUFFinbrivers")
covset<-c("black_","ownhome_","house_value_","foreign_born_"); yrs<-c(1910,1920,1930)   # 4 principal covariates
read_map<-function(file,idcol,pfx,is_grid){d<-read_dta(paste0(ahm_path,file))%>%filter(year==2010)%>%mutate(tract11=a2f(trctID2010))
  if(is_grid) d<-d%>%rename_with(~sub("BUFF_CFin","BUFFin",.x))
  d%>%select(tract11,matches(paste0("^",idcol,"[0-9]+$")),matches("^boundary_side[0-9]+$|^side_grid[0-9]+$"),matches("^holc_grade[0-9]+$"),
      matches("^boundary_grades[0-9]+$"),matches("^BUFFintrain[0-9]+$|^BUFFinsrivers[0-9]+$|^BUFFinbrivers[0-9]+$"))%>%
    rename_with(~sub("^side_grid","boundary_side",.x))%>%rename_with(~sub(paste0("^",idcol),"boundary_id",.x))%>%
    pivot_longer(-tract11,names_to=c(".value","slot"),names_pattern="^(boundary_id|boundary_side|holc_grade|boundary_grades|BUFFintrain|BUFFinsrivers|BUFFinbrivers)([0-9]+)$")%>%
    drop_na(boundary_id)%>%mutate(boundary_id=paste0(pfx,boundary_id))}
mA<-read_map("redlining_dataset_t2010_4th_new.dta","boundary_id","A_",FALSE)
mG<-read_map("redlining_dataset_tgrid2010_4th_new.dta","boundary_id_grid","G_",TRUE)
raw<-read_dta(paste0(ahm_path,"redlining_dataset_t2010_4th_new.dta"))%>%filter(year%in%yrs)%>%mutate(tract11=a2f(trctID2010))
cov_all<-raw%>%select(tract11,year,any_of(covset))%>%group_by(tract11,year)%>%summarise(across(everything(),~mean(.x,na.rm=T)),.groups="drop")%>%
  pivot_wider(names_from=year,values_from=any_of(covset),names_glue="{.value}{year}"); cvars<-setdiff(names(cov_all),"tract11")
tr<-mA%>%filter(boundary_grades%in%c("2_to_3","2_to_4"))%>%mutate(real=1L); pl<-mG%>%filter(boundary_grades%in%c("2_to_2","3_to_3","4_to_4"))%>%mutate(real=0L)
bt<-bind_rows(tr,pl)%>%left_join(cov_all,by="tract11")%>%mutate(barrier_amt=rowSums(across(all_of(BAR),~replace_na(.x,0))))
bt<-bt%>%semi_join(bt%>%group_by(boundary_id)%>%summarise(b=as.integer(max(barrier_amt,na.rm=T)>0.05))%>%filter(b==0),by="boundary_id")
ebalance<-function(X,target,maxit=300,tol=1e-8){Z<-sweep(X,2,target,"-");n<-nrow(Z);base<-rep(1/n,n);lam<-rep(0,ncol(Z));w<-base
  for(it in 1:maxit){u<-as.vector(Z%*%lam);u<-u-max(u);w<-base*exp(u);w<-w/sum(w);g<-colSums(w*Z);if(max(abs(g))<tol)break
    H<-t(Z)%*%(w*Z)-outer(g,g);step<-tryCatch(solve(H+diag(1e-8,ncol(Z)),g),error=function(e)g);best<-max(abs(g));lam0<-lam
    for(s in c(1,.5,.25,.1,.03)){lt<-lam0-s*step;ut<-as.vector(Z%*%lt);ut<-ut-max(ut);wt<-base*exp(ut);wt<-wt/sum(wt);gt<-colSums(wt*Z);if(max(abs(gt))<best){lam<-lt;best<-max(abs(gt));break}};if(identical(lam,lam0))break};w}
draw<-function(seed,TARGETVARS){
  side<-bt%>%group_by(boundary_id,real,boundary_side)%>%summarise(grade=first(holc_grade[holc_grade%in%c("A","B","C","D")]),
    across(all_of(cvars),~mean(.x,na.rm=T)),.groups="drop")%>%group_by(boundary_id)%>%filter(n()==2)%>%ungroup()
  set.seed(seed); side<-side%>%group_by(boundary_id)%>%mutate(lgs=if(first(real)==1L)as.integer(grade%in%c("C","D"))else{u<-runif(n());as.integer(u==max(u))})%>%ungroup()
  gap<-side%>%group_by(boundary_id,real)%>%summarise(across(all_of(cvars),~mean(.x[lgs==1],na.rm=T)-mean(.x[lgs==0],na.rm=T)),.groups="drop")
  gs<-gap[,"real"];cons<-c(); for(v in TARGETVARS){z<-as.numeric(scale(gap[[v]]));gs[[v]]<-ifelse(is.finite(z),z,0);cons<-c(cons,v)}
  tgt<-colMeans(gs[gs$real==1,cons,drop=FALSE]); wp<-ebalance(as.matrix(gs[gs$real==0,cons,drop=FALSE]),tgt)*sum(gs$real==0)
  w<-ifelse(gap$real==1L,1,NA); w[gap$real==0]<-wp
  # return, for each cvar: real gap and reweighted-placebo gap
  sapply(cvars,function(v){c(real=mean(gap[[v]][gap$real==1],na.rm=T), plac=weighted.mean(gap[[v]][gap$real==0],w[gap$real==0],na.rm=T))})}
avg<-function(TARGETVARS){R<-lapply(c(1,3,7,13,21),draw,TARGETVARS=TARGETVARS); Reduce(`+`,R)/length(R)}
show<-function(M,lab){cat("\n==",lab,"==\n  covariate-year        real gap   reweighted-placebo gap\n")
  for(v in cvars) cat(sprintf("  %-20s  %+8.4f   %+8.4f\n",v,M["real",v],M["plac",v]))}
key1930<-paste0(c("black_","ownhome_","house_value_","foreign_born_"),"1930")
cat("Within-boundary covariate GAPS (lower- minus higher-graded side), real vs entropy-reweighted placebos.\n")
cat("If placebos match real at 1910/1920 too, pre-map trajectories are parallel.\n")
show(avg(cvars),      "(a) BALANCE ON ALL YEARS 1910-1930 (main spec: gaps matched by construction)")
show(avg(key1930),    "(b) LEAVE-OUT: balance on 1930 ONLY -> are held-out 1910 & 1920 gaps matched?")
