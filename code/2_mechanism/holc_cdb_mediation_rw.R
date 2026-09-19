# Mediation of the ENTROPY-BALANCED CD-B boundary DiD.
# Integration segregation is benchmarked to the METRO (MSA) composition, so the correct
# mediator is the DISTANCE of neighborhood composition from the MSA: M = |Black_nbhd - Black_MSA|.
# (Raw Black share would be wrong: raising Black share only raises segregation if it moves the
#  neighborhood FURTHER from the metro average. We verify below that it does.)
# DiD mediation, treatment = lgs:real (redline boundary net of placebo), entropy-balanced weights:
#   c  total  : seg ~ lgs + lgs:real
#   a         : |dev| ~ lgs + lgs:real
#   b, c' dir : seg ~ lgs + lgs:real + |dev|
#   indirect = a*b (Sobel SE); prop mediated = a*b/c
suppressMessages({library(tidyverse); library(haven); library(fixest)})
# EB_SET="all": baseline -- all 1910-1930 covariate gaps + missingness indicators (as in
# holc_cdb_ebal.R). "key4": original four key 1930 covariates.
EB_SET<-"all"
ahm_path<-"./data/ahm_staging/"; a2f<-function(x)paste0(substr(x,2,3),substr(x,5,7),substr(x,9,14))
BAR<-c("BUFFintrain","BUFFinsrivers","BUFFinbrivers")
covset<-c("black_","ownhome_","foreign_born_","read_write_","house_value_","rent_","radio_"); yrs<-c(1910,1920,1930)
load("./data/safegraph_2019_m.rdata")
# MSA Black share (population-weighted over distinct block groups), then distance from it
msa_comp<-places_usa_2019%>%distinct(poi_cbg,.keep_all=TRUE)%>%group_by(msa)%>%summarise(msa_black=weighted.mean(black_rate,tot_pop,na.rm=TRUE),.groups="drop")
poi<-places_usa_2019%>%left_join(msa_comp,by="msa")%>%left_join(readRDS("data/poi_seg_reweighted.rds")%>%transmute(safegraph_place_id,seg_int_rw),by="safegraph_place_id")%>%mutate(seg_int=seg_int_rw*5/8,lvis=log(raw_visit_counts+1),
  tract11=substr(poi_cbg,1,11),dev=black_rate-msa_black,absdev=abs(black_rate-msa_black))%>%
  filter(!is.na(seg_int),!is.na(lvis),!is.na(black_rate),!is.na(msa_black))%>%select(tract11,sub_category,seg_int,lvis,black_rate,msa_black,dev,absdev)
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
key<-intersect(paste0(c("black_","ownhome_","house_value_","foreign_born_"),"1930"),cvars)
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
  bvars<-if(EB_SET=="all") cvars[sapply(cvars,function(v)sum(is.finite(as.numeric(scale(gap[[v]]))))>100)] else key
  gs<-gap[,"real"];cons<-c()
  for(v in bvars){z<-as.numeric(scale(gap[[v]]));gs[[paste0("Z_",v)]]<-ifelse(is.finite(z),z,0);cons<-c(cons,paste0("Z_",v))
    if(EB_SET=="all" && any(!is.finite(z))){gs[[paste0("MI_",v)]]<-as.numeric(!is.finite(z));cons<-c(cons,paste0("MI_",v))}}
  tgt<-colMeans(gs[gs$real==1,cons,drop=FALSE]);wp<-ebalance(as.matrix(gs[gs$real==0,cons,drop=FALSE]),tgt)*sum(gs$real==0)
  gg<-gap%>%mutate(w_att=ifelse(real==1L,1,NA));gg$w_att[gg$real==0]<-wp
  side%>%select(boundary_id,boundary_side,lgs,real)%>%inner_join(gg%>%select(boundary_id,w_att),by="boundary_id")%>%
    inner_join(bt%>%distinct(boundary_id,boundary_side,tract11),by=c("boundary_id","boundary_side"))%>%inner_join(poi,by="tract11",relationship="many-to-many")%>%mutate(w=w_att*lvis)}
cf<-function(m,t)c(b=unname(coef(m)[t]),s=unname(se(m)[t]))
mediate<-function(dat){mc<-feols(seg_int~lgs+lgs:real|boundary_id+sub_category,dat,weights=~w,cluster=~boundary_id)
  ma<-feols(absdev~lgs+lgs:real|boundary_id+sub_category,dat,weights=~w,cluster=~boundary_id)
  mb<-feols(seg_int~lgs+lgs:real+absdev|boundary_id+sub_category,dat,weights=~w,cluster=~boundary_id)
  c<-cf(mc,"lgs:real");a<-cf(ma,"lgs:real");b<-cf(mb,"absdev");cp<-cf(mb,"lgs:real")
  ind<-a["b"]*b["b"];se_ind<-sqrt(b["b"]^2*a["s"]^2+a["b"]^2*b["s"]^2)
  list(c=c,a=a,b=b,cp=cp,ind=unname(ind),z=unname(ind/se_ind),prop=unname(ind/c["b"]))}
cat("Mediator = |neighborhood Black share - MSA Black share| (distance from metro composition).\n")
d1<-mkdat(1)
cat("\n== Composition relative to MSA, real boundaries (redlined lgs=1 vs better lgs=0) ==\n")
print(d1%>%filter(real==1)%>%group_by(lgs)%>%summarise(nbhd_black=weighted.mean(black_rate,w),msa_black=weighted.mean(msa_black,w),
  signed_dev=weighted.mean(dev,w),abs_dev=weighted.mean(absdev,w),.groups="drop")%>%mutate(across(-lgs,~round(.x,4)))%>%as.data.frame(),row.names=FALSE)
cat("(redlined side sits further ABOVE the MSA Black share -> larger distance -> more segregation)\n")
for(sd in c(1,3,13)){r<-mediate(mkdat(sd));cat(sprintf("\n-- seed %d --  c=%+.5f | a=%+.5f(se%.5f) | b=%+.5f(se%.5f) | c'=%+.5f | indirect=%+.5f (z=%.2f, p=%.4f) | prop mediated=%.0f%%\n",
  sd,r$c["b"],r$a["b"],r$a["s"],r$b["b"],r$b["s"],r$cp["b"],r$ind,r$z,2*pnorm(-abs(r$z)),100*r$prop))}
