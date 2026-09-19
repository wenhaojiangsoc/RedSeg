# AMENITY-COMPOSITION channel of the entropy-balanced CD-B boundary DiD.
# Question: does redlining shift the MIX of establishment types toward categories that are
# intrinsically more segregated, rather than (only) making places more segregated within type?
# The main spec's sub_category FE absorbs this channel entirely, so decompose
#   seg_int = m_cat + (seg_int - m_cat),  m_cat = national visit-weighted mean integration
# segregation of the POI's sub_category, and run the boundary DiD WITHOUT category FE on each
# piece. By linearity (same sample/weights/FE) the two DiDs sum exactly to the total:
#   total (no cat FE) = amenity-mix part (DiD on m_cat) + within-category part (DiD on resid).
# Compare with the racial-composition channel (~2/5 of the within effect, holc_cdb_mediation_rw.R).
# Design (boundaries, barriers, EB weights, seeds 1/3/13) identical to holc_cdb_mediation_rw.R.
suppressMessages({library(tidyverse); library(haven); library(fixest)})
EB_SET<-"all"
ahm_path<-"./data/ahm_staging/"; a2f<-function(x)paste0(substr(x,2,3),substr(x,5,7),substr(x,9,14))
BAR<-c("BUFFintrain","BUFFinsrivers","BUFFinbrivers")
covset<-c("black_","ownhome_","foreign_born_","read_write_","house_value_","rent_","radio_"); yrs<-c(1910,1920,1930)
load("./data/safegraph_2019_m.rdata")
poi<-places_usa_2019%>%left_join(readRDS("data/poi_seg_reweighted.rds")%>%transmute(safegraph_place_id,seg_int_rw),by="safegraph_place_id")%>%
  mutate(seg_int=seg_int_rw*5/8,lvis=log(raw_visit_counts+1),tract11=substr(poi_cbg,1,11))%>%
  filter(!is.na(seg_int),!is.na(lvis),!is.na(sub_category))%>%select(tract11,sub_category,seg_int,lvis)
catm<-poi%>%group_by(sub_category)%>%summarise(m_cat=weighted.mean(seg_int,lvis),n_cat=n(),.groups="drop")
cat(sprintf("categories: %d | national m_cat range: %.3f - %.3f\n",nrow(catm),min(catm$m_cat),max(catm$m_cat)))
poi<-poi%>%left_join(catm%>%select(sub_category,m_cat),by="sub_category")%>%mutate(resid=seg_int-m_cat)
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
cf<-function(m,t)c(b=unname(coef(m)[t]),s=unname(se(m)[t]),p=unname(pvalue(m)[t]))
decomp<-function(dat){
  mt<-feols(seg_int~lgs+lgs:real|boundary_id,dat,weights=~w,cluster=~boundary_id)              # total, no category FE
  mm<-feols(m_cat  ~lgs+lgs:real|boundary_id,dat,weights=~w,cluster=~boundary_id)              # amenity-mix part
  mr<-feols(resid  ~lgs+lgs:real|boundary_id,dat,weights=~w,cluster=~boundary_id)              # within-category part
  mf<-feols(seg_int~lgs+lgs:real|boundary_id+sub_category,dat,weights=~w,cluster=~boundary_id) # paper spec (cat FE)
  list(tot=cf(mt,"lgs:real"),amen=cf(mm,"lgs:real"),with=cf(mr,"lgs:real"),paper=cf(mf,"lgs:real"))}
d1<-mkdat(1)
cat("\n== Amenity mix on REAL boundaries (redlined lgs=1 vs better lgs=0), visit-weighted ==\n")
print(d1%>%filter(real==1)%>%group_by(lgs)%>%summarise(exp_seg_of_mix=weighted.mean(m_cat,w),.groups="drop")%>%
  mutate(exp_seg_of_mix=round(exp_seg_of_mix,4))%>%as.data.frame(),row.names=FALSE)
cat("\n== Top category share shifts on real boundaries (redlined minus better side, pp of visits) ==\n")
sh<-d1%>%filter(real==1)%>%group_by(lgs,sub_category)%>%summarise(v=sum(w),.groups="drop")%>%group_by(lgs)%>%mutate(s=v/sum(v))%>%
  select(-v)%>%pivot_wider(names_from=lgs,values_from=s,values_fill=0,names_prefix="s")%>%mutate(d=s1-s0)%>%
  left_join(catm%>%select(sub_category,m_cat),by="sub_category")
print(bind_rows(sh%>%slice_max(d,n=5),sh%>%slice_min(d,n=5))%>%mutate(across(c(s0,s1,d,m_cat),~round(.x,4)))%>%as.data.frame(),row.names=FALSE)
for(sd in c(1,3,13)){r<-decomp(if(sd==1)d1 else mkdat(sd))
  cat(sprintf("\n-- seed %d --\n  total (no cat FE)    %+.5f (se %.5f, p=%.3f)\n  amenity-mix part     %+.5f (se %.5f, p=%.3f)  -> share of total %.0f%%\n  within-category part %+.5f (se %.5f, p=%.3f)\n  paper spec (cat FE)  %+.5f (se %.5f, p=%.3f)\n",
    sd,r$tot["b"],r$tot["s"],r$tot["p"],r$amen["b"],r$amen["s"],r$amen["p"],100*r$amen["b"]/r$tot["b"],
    r$with["b"],r$with["s"],r$with["p"],r$paper["b"],r$paper["s"],r$paper["p"]))}
