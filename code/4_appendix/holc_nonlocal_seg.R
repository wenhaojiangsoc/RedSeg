# NON-LOCAL robustness (step 1 of 2): recompute integration segregation EXCLUDING local residents
# (visitors from the POI's own CBG, and own tract) to test whether the redlining effect is mechanical
# (local residents visiting local places) or reflects genuine cross-neighborhood sorting.
# Per POI, from raw flows: visitor race composition for all / excl-own-CBG / excl-own-tract; then
# segregation = (5/8) * sum_g |pi_g - delta_g^MSA|. Saves results/poi_nonlocal_seg.rds. DiD is step 2.
suppressMessages({library(data.table)})
grp <- c("black","white","hispanic","asian","other")
cr <- fread("data/cbg_race_2019.csv", colClasses=list(character="GEOID"))
cr[, s := white+black+hispanic+asian+other]
for(g in grp) cr[, (paste0("r_",g)) := fifelse(s>0, get(g)/s, NA_real_)]
cr <- cr[s>0, c("GEOID","tot_pop",paste0("r_",grp)), with=FALSE]; setnames(cr,"GEOID","cbg"); setkey(cr,cbg)
load("data/safegraph_2019_m.rdata")
lk <- unique(as.data.table(places_usa_2019)[,.(safegraph_place_id, poi_cbg=as.character(poi_cbg), msa)], by="safegraph_place_id")
lk[, poi_tract := substr(poi_cbg,1,11)]; setkey(lk,safegraph_place_id)
cbg2msa <- unique(lk[,.(cbg=poi_cbg,msa)], by="cbg")
msad <- merge(cr, cbg2msa, by="cbg")[, lapply(.SD,function(x) weighted.mean(x,tot_pop,na.rm=TRUE)), by=msa, .SDcols=paste0("r_",grp)]
setnames(msad, paste0("r_",grp), paste0("d_",grp)); setkey(msad,msa)
vcols <- paste0("v_",grp)
addup <- function(a,b){ if(is.null(a)) return(b); rbindlist(list(a,b))[, lapply(.SD,sum), by=safegraph_place_id] }
accA<-accO<-accT<-NULL
for(f in list.files("data/flows_2019","\\.csv\\.gz$",full.names=TRUE)){
  dt <- fread(cmd=paste("gunzip -c",shQuote(f)), colClasses=list(character=c("safegraph_place_id","cbg")))
  dt <- lk[dt, on="safegraph_place_id", nomatch=0L]
  dt <- cr[dt, on="cbg", nomatch=0L]
  for(g in grp) dt[, (paste0("v_",g)) := visitors*get(paste0("r_",g))]
  sc <- c("visitors",vcols)
  accA <- addup(accA, dt[, lapply(.SD,sum), by=safegraph_place_id, .SDcols=sc])
  accO <- addup(accO, dt[cbg==poi_cbg,               lapply(.SD,sum), by=safegraph_place_id, .SDcols=sc])
  accT <- addup(accT, dt[substr(cbg,1,11)==poi_tract, lapply(.SD,sum), by=safegraph_place_id, .SDcols=sc])
  cat("done",basename(f),"POIs:",nrow(accA),"\n"); flush.console()
}
setnames(accA,c("visitors",vcols),c("t",paste0("a_",grp)))
setnames(accO,c("visitors",vcols),c("toc",paste0("oc_",grp)))
setnames(accT,c("visitors",vcols),c("tot",paste0("ot_",grp)))
P <- merge(merge(accA,accO,by="safegraph_place_id",all.x=TRUE),accT,by="safegraph_place_id",all.x=TRUE)
for(c in names(P)) if(is.numeric(P[[c]])) set(P, which(is.na(P[[c]])), c, 0)
P <- merge(P, lk[,.(safegraph_place_id,poi_cbg,msa)], by="safegraph_place_id")
P <- merge(P, msad, by="msa")
D <- as.matrix(P[, paste0("d_",grp), with=FALSE])
segfun <- function(numer, denom){ pi <- numer/denom; rowSums(abs(pi - D))*5/8 }
Aall<-as.matrix(P[,paste0("a_",grp),with=FALSE]); Aoc<-as.matrix(P[,paste0("oc_",grp),with=FALSE]); Aot<-as.matrix(P[,paste0("ot_",grp),with=FALSE])
P[, `:=`(seg_all=segfun(Aall,t), seg_excbg=segfun(Aall-Aoc,t-toc), seg_extract=segfun(Aall-Aot,t-tot),
         share_cbg=toc/t, share_tract=tot/t, tract11=substr(poi_cbg,1,11))]
P <- P[t>0 & is.finite(seg_all)]
saveRDS(P[,.(safegraph_place_id,tract11,seg_all,seg_excbg,seg_extract,share_cbg,share_tract,t)], "results/poi_nonlocal_seg.rds")
val <- as.data.table(places_usa_2019)[P, on="safegraph_place_id", race_seg_msa_adj_mean]
cat(sprintf("\n=== NON-LOCAL recompute over %d POIs ===\n", nrow(P)))
cat(sprintf("same-CBG visit share:   mean=%.3f median=%.3f\n", mean(P$share_cbg),median(P$share_cbg)))
cat(sprintf("same-tract visit share: mean=%.3f median=%.3f\n", mean(P$share_tract),median(P$share_tract)))
cat(sprintf("seg all=%.3f | excl own CBG=%.3f | excl own tract=%.3f\n", mean(P$seg_all),mean(P$seg_excbg,na.rm=T),mean(P$seg_extract,na.rm=T)))
cat(sprintf("VALIDATION corr(seg_all, main race_seg_msa_adj_mean) = %.3f\n", cor(P$seg_all, val, use="complete.obs")))
