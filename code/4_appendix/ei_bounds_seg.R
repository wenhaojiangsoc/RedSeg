# ECOLOGICAL-INFERENCE robustness (step 1 of 2): recompute the POI integration-segregation index under
# (a) HOMOGENEITY restriction -- use only flows from source CBGs whose largest race share >= tau, where the
#     ecological attribution is (near-)exact and within-CBG differential POI-sorting cannot bias it; and
# (b) DUNCAN-DAVIS method-of-bounds -- per-POI seg_lower/seg_upper from the deterministic bounds on race-g
#     counts, assuming CBG i's race-g outflow = r_gi * O_i (equal per-capita travel by race within a CBG),
#     O_i = total observed outflow of CBG i (pass 1). Saves results/poi_ei_bounds.rds. DiD is step 2.
suppressMessages({library(data.table)})
# run from the project root; # setwd("~/Library/CloudStorage/Dropbox/Demography-Rep-2026")
grp <- c("black","white","hispanic","asian","other")
cr <- fread("data/cbg_race_2019.csv", colClasses=list(character="GEOID"))
cr[, s := white+black+hispanic+asian+other]
for(g in grp) cr[, (paste0("r_",g)) := fifelse(s>0, get(g)/s, NA_real_)]
cr <- cr[s>0, c("GEOID","tot_pop",paste0("r_",grp)), with=FALSE]; setnames(cr,"GEOID","cbg")
cr[, maxrace := pmax(r_black,r_white,r_hispanic,r_asian,r_other)]; setkey(cr,cbg)
flist <- list.files("data/flows_2019","\\.csv\\.gz$",full.names=TRUE)

## PASS 1: total outflow O_i per source CBG
O <- NULL
for(f in flist){
  dt <- fread(cmd=paste("gunzip -c",shQuote(f)), select=c("cbg","visitors"), colClasses=list(character="cbg"))
  o <- dt[, .(O=sum(visitors)), by=cbg]; O <- if(is.null(O)) o else rbindlist(list(O,o))[, .(O=sum(O)), by=cbg]
  cat("pass1",basename(f),"\n"); flush.console()
}
setkey(O,cbg); cr <- O[cr, on="cbg"]; cr[is.na(O), O:=0]
for(g in grp) cr[, (paste0("cap_",g)) := get(paste0("r_",g))*O]     # race-g outflow cap

## POI -> MSA lookup + MSA racial benchmark d_g (population-weighted), as in holc_nonlocal_seg.R
load("data/safegraph_2019_m.rdata")
lk <- unique(as.data.table(places_usa_2019)[,.(safegraph_place_id, poi_cbg=as.character(poi_cbg), msa)], by="safegraph_place_id"); setkey(lk,safegraph_place_id)
cbg2msa <- unique(lk[,.(cbg=poi_cbg,msa)], by="cbg")
msad <- merge(cr, cbg2msa, by="cbg")[, lapply(.SD,function(x) weighted.mean(x,tot_pop,na.rm=TRUE)), by=msa, .SDcols=paste0("r_",grp)]
setnames(msad, paste0("r_",grp), paste0("d_",grp)); setkey(msad,msa)

## PASS 2: accumulate per POI -- eco a_g, bound max/min, homogeneity-restricted a_g and t at tau in {0.7,0.9}
taus <- c(0.7,0.9); acc <- NULL
for(f in flist){
  dt <- fread(cmd=paste("gunzip -c",shQuote(f)), colClasses=list(character=c("safegraph_place_id","cbg")))
  dt <- lk[dt, on="safegraph_place_id", nomatch=0L]; dt <- cr[dt, on="cbg", nomatch=0L]
  for(g in grp){
    dt[, (paste0("v_",g)) := visitors*get(paste0("r_",g))]
    dt[, (paste0("mx_",g)) := pmin(visitors, get(paste0("cap_",g)))]
    dt[, (paste0("mn_",g)) := pmax(0, visitors - (O - get(paste0("cap_",g))))]
  }
  for(tau in taus){ hm <- as.numeric(dt$maxrace >= tau)
    for(g in grp) dt[, (paste0("h",tau,"_",g)) := get(paste0("v_",g))*hm]
    dt[, (paste0("th",tau)) := visitors*hm] }
  cols <- c("visitors", paste0("v_",grp), paste0("mx_",grp), paste0("mn_",grp),
            unlist(lapply(taus,function(tau) c(paste0("h",tau,"_",grp), paste0("th",tau)))))
  a <- dt[, lapply(.SD,sum), by=safegraph_place_id, .SDcols=cols]
  acc <- if(is.null(acc)) a else rbindlist(list(acc,a))[, lapply(.SD,sum), by=safegraph_place_id]
  cat("pass2",basename(f),"POIs:",nrow(acc),"\n"); flush.console()
}
P <- merge(acc, lk[,.(safegraph_place_id,poi_cbg,msa)], by="safegraph_place_id"); P <- merge(P, msad, by="msa")
P[, tract11 := substr(poi_cbg,1,11)]; setnames(P,"visitors","t")
D <- as.matrix(P[, paste0("d_",grp), with=FALSE])
P[, seg_eco := rowSums(abs(as.matrix(.SD)/t - D))*5/8, .SDcols=paste0("v_",grp)]
Mx <- as.matrix(P[,paste0("mx_",grp),with=FALSE])/P$t; Mn <- as.matrix(P[,paste0("mn_",grp),with=FALSE])/P$t
P[, seg_upper := rowSums(pmax(abs(Mn-D),abs(Mx-D)))*5/8]
P[, seg_lower := rowSums(ifelse(D>=Mn & D<=Mx, 0, pmin(abs(Mn-D),abs(Mx-D))))*5/8]
for(tau in taus){ th<-P[[paste0("th",tau)]]; Ah<-as.matrix(P[,paste0("h",tau,"_",grp),with=FALSE])
  P[[paste0("seg_h",tau)]] <- fifelse(th>0, rowSums(abs(Ah/th-D))*5/8, NA_real_)
  P[[paste0("share_h",tau)]] <- th/P$t }
P <- P[t>0 & is.finite(seg_eco)]
saveRDS(P[, c("safegraph_place_id","tract11","t","seg_eco","seg_lower","seg_upper",
              paste0("seg_h",taus),paste0("share_h",taus)), with=FALSE], "results/poi_ei_bounds.rds")
cat(sprintf("\nPOIs %d | seg_eco=%.4f  [seg_lower=%.4f, seg_upper=%.4f]\n", nrow(P),mean(P$seg_eco),mean(P$seg_lower),mean(P$seg_upper)))
for(tau in taus) cat(sprintf("tau=%.1f homogeneity: seg=%.4f | share of visits from homog CBGs (mean)=%.3f\n",
    tau, mean(P[[paste0("seg_h",tau)]],na.rm=T), mean(P[[paste0("share_h",tau)]])))
if(file.exists("results/poi_nonlocal_seg.rds")){ nl<-as.data.table(readRDS("results/poi_nonlocal_seg.rds"))
  m<-merge(P[,.(safegraph_place_id,seg_eco)], nl[,.(safegraph_place_id,seg_all)], by="safegraph_place_id")
  cat(sprintf("VALIDATION corr(seg_eco, published seg_all)=%.4f (should be ~1.00)\n", cor(m$seg_eco,m$seg_all))) }
