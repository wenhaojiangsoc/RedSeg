# ECOLOGICAL-INFERENCE robustness on the FRACTIONAL CBG racial composition (ACS B03002).
# Memory-efficient: all per-POI sums computed INLINE in the by-aggregation (no per-row intermediate columns).
# Outputs seg_eco (fractional ecological measure), seg_lower/seg_upper (Duncan-Davis bounds under
# race-g outflow = r_gi*O_i), seg_h08 (homogeneity tau>=0.8). Saves results/poi_ei_bounds_frac.rds.
suppressMessages({library(data.table)}); try(mem.maxVSize(1.2e6), silent=TRUE)   # raise ceiling if allowed
# run from the project root; # setwd("~/Library/CloudStorage/Dropbox/Demography-Rep-2026")
grp <- c("black","white","hispanic","asian","other")
cr <- fread("data/cbg_race_frac_2019.csv", colClasses=list(character="GEOID"))
cr[, s := rowSums(.SD), .SDcols=grp]
for(g in grp) cr[, (paste0("r_",g)) := fifelse(s>0, get(g)/s, NA_real_)]
cr <- cr[s>0, c("GEOID","tot_pop",paste0("r_",grp)), with=FALSE]; setnames(cr,"GEOID","cbg")
cr[, maxrace := pmax(r_black,r_white,r_hispanic,r_asian,r_other)]; setkey(cr,cbg)
flist <- list.files("data/flows_2019","\\.csv\\.gz$",full.names=TRUE)

## PASS 1: total outflow O_i per source CBG
O <- NULL
for(f in flist){ dt <- fread(cmd=paste("gunzip -c",shQuote(f)), select=c("cbg","visitors"), colClasses=list(character="cbg"))
  o <- dt[, .(O=sum(visitors)), by=cbg]; O <- if(is.null(O)) o else rbindlist(list(O,o))[, .(O=sum(O)), by=cbg]
  rm(dt,o); gc(); cat("pass1",basename(f),"\n"); flush.console() }
setkey(O,cbg); cr <- O[cr, on="cbg"]; cr[is.na(O), O:=0]
for(g in grp) cr[, (paste0("cap_",g)) := get(paste0("r_",g))*O]

## POI -> MSA lookup + MSA benchmark d_g
load("data/safegraph_2019_m.rdata")
lk <- unique(as.data.table(places_usa_2019)[,.(safegraph_place_id, poi_cbg=as.character(poi_cbg), msa)], by="safegraph_place_id"); setkey(lk,safegraph_place_id)
cbg2msa <- unique(lk[,.(cbg=poi_cbg,msa)], by="cbg")
msad <- merge(cr, cbg2msa, by="cbg")[, lapply(.SD,function(x) weighted.mean(x,tot_pop,na.rm=TRUE)), by=msa, .SDcols=paste0("r_",grp)]
setnames(msad, paste0("r_",grp), paste0("d_",grp)); setkey(msad,msa)
rm(places_usa_2019); gc()

## PASS 2: inline aggregation (a_g, mx_g, mn_g, h_g, th, t) per POI, tau=0.8
parts <- "t=sum(visitors)"
for(g in grp) parts <- c(parts,
  sprintf("a_%s=sum(visitors*r_%s)",g,g),
  sprintf("mx_%s=sum(pmin(visitors, cap_%s))",g,g),
  sprintf("mn_%s=sum(pmax(0, visitors-(O-cap_%s)))",g,g),
  sprintf("h_%s=sum(visitors*r_%s*(maxrace>=0.8))",g,g))
parts <- c(parts, "th=sum(visitors*(maxrace>=0.8))")
J <- parse(text=paste0("list(",paste(parts,collapse=", "),")"))
crsub <- cr[, c("cbg","maxrace","O",paste0("r_",grp),paste0("cap_",grp)), with=FALSE]; setkey(crsub,cbg)
acc <- NULL
for(f in flist){
  dt <- fread(cmd=paste("gunzip -c",shQuote(f)), colClasses=list(character=c("safegraph_place_id","cbg")))
  dt <- crsub[dt, on="cbg", nomatch=0L]
  a <- dt[, eval(J), by=safegraph_place_id]
  acc <- if(is.null(acc)) a else rbindlist(list(acc,a))[, lapply(.SD,sum), by=safegraph_place_id]
  rm(dt,a); gc(); cat("pass2",basename(f),"POIs:",nrow(acc),"\n"); flush.console()
}
P <- merge(acc, lk[,.(safegraph_place_id,poi_cbg,msa)], by="safegraph_place_id"); P <- merge(P, msad, by="msa")
P[, tract11 := substr(poi_cbg,1,11)]
D <- as.matrix(P[, paste0("d_",grp), with=FALSE])
P[, seg_eco := rowSums(abs(as.matrix(.SD)/t - D))*5/8, .SDcols=paste0("a_",grp)]
Mx <- as.matrix(P[,paste0("mx_",grp),with=FALSE])/P$t; Mn <- as.matrix(P[,paste0("mn_",grp),with=FALSE])/P$t
P[, seg_upper := rowSums(pmax(abs(Mn-D),abs(Mx-D)))*5/8]
P[, seg_lower := rowSums(ifelse(D>=Mn & D<=Mx, 0, pmin(abs(Mn-D),abs(Mx-D))))*5/8]
Ah <- as.matrix(P[,paste0("h_",grp),with=FALSE])
P[, seg_h08 := fifelse(th>0, rowSums(abs(Ah/th - D))*5/8, NA_real_)]
P[, share_h08 := th/t]
P <- P[t>0 & is.finite(seg_eco)]
saveRDS(P[,.(safegraph_place_id,tract11,t,seg_eco,seg_lower,seg_upper,seg_h08,share_h08)], "results/poi_ei_bounds_frac.rds")
cat(sprintf("\nFRACTIONAL: POIs %d | seg_eco=%.4f  [lower=%.4f, upper=%.4f] | seg_h08=%.4f (homog visit share=%.3f)\n",
    nrow(P),mean(P$seg_eco),mean(P$seg_lower),mean(P$seg_upper),mean(P$seg_h08,na.rm=T),mean(P$share_h08)))
nl<-as.data.table(readRDS("results/poi_nonlocal_seg.rds"))
m<-merge(P[,.(safegraph_place_id,seg_eco)], nl[,.(safegraph_place_id,seg_all)], by="safegraph_place_id")
cat(sprintf("corr(fractional seg_eco, majority-race seg_all)=%.4f\n", cor(m$seg_eco,m$seg_all)))
