# ECOLOGICAL INFERENCE via GOODMAN (1953) ecological regression -- the model underlying King's (1997) EI --
# applied ESTABLISHMENT BY ESTABLISHMENT on the fractional CBG racial composition (ACS B03002).
#
# Neighborhood model (the paper's proportional measure): every race in source CBG i visits POI j at the same
#   per-capita rate, so race-g visitors = sum_i n_ij * r_gi.
# Goodman/constancy model (this script): each race g has its OWN visitation propensity beta_gj to POI j, constant
#   across the POI's source CBGs; identified from cross-CBG variation among senders:
#   n_ijt / O_it = sum_g beta_gj * r_gi   (unit = sender CBG-month, weight O_it = CBG i's total outflow in month t)
#   WLS normal equations:  S beta = a,  S_gh = sum O_it r_gi r_hi,  a_g = sum r_gi n_ijt (= the proportional count).
#   Solved with beta >= 0 (NNLS / QP). Implied race-g visitor count = beta_g X_g, X_g = sum r_gi O_it over senders;
#   sum_g beta_g X_g = t exactly at the unconstrained WLS solution (regressors sum to 1); normalized to shares.
# Falls back to the proportional composition when the POI has < MINSM sender-months (too few observations).
# Output: results/poi_ei_goodman.rds with seg_goodman (fallback-filled) and seg_goodman_strict (NA if fallback).
suppressMessages({library(data.table); library(quadprog)}); try(mem.maxVSize(1.2e6), silent=TRUE)
# run from the project root; # setwd("~/Library/CloudStorage/Dropbox/Demography-Rep-2026")
MINSM <- 24
grp <- c("black","white","hispanic","asian","other")
cr <- fread("data/cbg_race_frac_2019.csv", colClasses=list(character="GEOID"))
cr[, s := rowSums(.SD), .SDcols=grp]
for(g in grp) cr[, (paste0("r_",g)) := fifelse(s>0, get(g)/s, NA_real_)]
cr <- cr[s>0, c("GEOID","tot_pop",paste0("r_",grp)), with=FALSE]; setnames(cr,"GEOID","cbg"); setkey(cr,cbg)
flist <- list.files("data/flows_2019","\\.csv\\.gz$",full.names=TRUE)

## POI -> MSA lookup + MSA benchmark d_g (as in ei_bounds_seg_frac.R)
load("data/safegraph_2019_m.rdata")
lk <- unique(as.data.table(places_usa_2019)[,.(safegraph_place_id, poi_cbg=as.character(poi_cbg), msa)], by="safegraph_place_id"); setkey(lk,safegraph_place_id)
cbg2msa <- unique(lk[,.(cbg=poi_cbg,msa)], by="cbg")
msad <- merge(cr, cbg2msa, by="cbg")[, lapply(.SD,function(x) weighted.mean(x,tot_pop,na.rm=TRUE)), by=msa, .SDcols=paste0("r_",grp)]
setnames(msad, paste0("r_",grp), paste0("d_",grp)); setkey(msad,msa)
rm(places_usa_2019); gc()

## ONE PASS: per POI accumulate t, nsm, a_g (5), X_g (5), S_gh (15) -- monthly O_it computed within each file
pairs <- combn(grp, 2, simplify=FALSE)
parts <- c("t=sum(visitors)", "nsm=.N")
for(g in grp) parts <- c(parts, sprintf("a_%s=sum(visitors*r_%s)",g,g), sprintf("X_%s=sum(O*r_%s)",g,g),
                         sprintf("S_%s_%s=sum(O*r_%s*r_%s)",g,g,g,g))
for(p in pairs) parts <- c(parts, sprintf("S_%s_%s=sum(O*r_%s*r_%s)",p[1],p[2],p[1],p[2]))
J <- parse(text=paste0("list(",paste(parts,collapse=", "),")"))
crsub <- cr[, c("cbg",paste0("r_",grp)), with=FALSE]; setkey(crsub,cbg)
acc <- NULL
for(f in flist){
  dt <- fread(cmd=paste("gunzip -c",shQuote(f)), colClasses=list(character=c("safegraph_place_id","cbg")))
  dt[, O := sum(visitors), by=cbg]                       # monthly outflow of the source CBG
  dt <- crsub[dt, on="cbg", nomatch=0L]
  a <- dt[, eval(J), by=safegraph_place_id]
  acc <- if(is.null(acc)) a else rbindlist(list(acc,a))[, lapply(.SD,sum), by=safegraph_place_id]
  rm(dt,a); gc(); cat("pass",basename(f),"POIs:",nrow(acc),"\n"); flush.console()
}

## Solve the 5x5 non-negative WLS per POI
Sn <- function(g,h){ n1<-paste0("S_",g,"_",h); if(n1 %in% names(acc)) n1 else paste0("S_",h,"_",g) }
Smat <- sapply(grp, function(g) sapply(grp, function(h) Sn(g,h)))           # 5x5 matrix of column names
A <- as.matrix(acc[, paste0("a_",grp), with=FALSE]); X <- as.matrix(acc[, paste0("X_",grp), with=FALSE])
Sarr <- as.matrix(acc[, unique(as.vector(Smat)), with=FALSE])
n <- nrow(acc); B <- matrix(NA_real_, n, 5); ok <- rep(FALSE, n)
idx <- matrix(match(Smat, colnames(Sarr)), 5, 5)
cat("solving", n, "QPs...\n"); flush.console()
for(i in seq_len(n)){
  if(acc$nsm[i] < MINSM) next
  S <- matrix(Sarr[i, idx], 5, 5); eps <- 1e-8*max(diag(S)); if(!is.finite(eps) || eps<=0) next
  sol <- tryCatch(solve.QP(Dmat=S+diag(eps,5), dvec=A[i,], Amat=diag(5), bvec=rep(0,5))$solution, error=function(e) NULL)
  if(is.null(sol)) next
  B[i,] <- pmax(sol,0); ok[i] <- TRUE
  if(i %% 500000 == 0){ cat("  ", i, "\n"); flush.console() }
}
V <- B*X                                                  # implied race-g visitor counts (Goodman)
V <- V/rowSums(V)                                         # shares
acc[, ok := ok]
P <- merge(acc[, .(safegraph_place_id, t, nsm, ok)], lk[,.(safegraph_place_id,poi_cbg,msa)], by="safegraph_place_id")
P <- merge(P, msad, by="msa"); P[, tract11 := substr(poi_cbg,1,11)]
m <- match(P$safegraph_place_id, acc$safegraph_place_id)
D <- as.matrix(P[, paste0("d_",grp), with=FALSE])
Peco <- A[m,]/acc$t[m]; Pg <- V[m,]
P[, seg_eco := rowSums(abs(Peco - D))*5/8]
sg <- rowSums(abs(Pg - D))*5/8
P[, seg_goodman_strict := fifelse(ok & is.finite(sg), sg, NA_real_)]
P[, seg_goodman := fifelse(is.na(seg_goodman_strict), seg_eco, seg_goodman_strict)]
P <- P[t>0 & is.finite(seg_eco)]
saveRDS(P[,.(safegraph_place_id,tract11,t,nsm,ok,seg_eco,seg_goodman,seg_goodman_strict)], "results/poi_ei_goodman.rds")
cat(sprintf("\nGOODMAN EI: POIs %d | solved %.1f%% | seg_eco=%.4f | seg_goodman=%.4f | seg_goodman_strict=%.4f (n=%d)\n",
    nrow(P), 100*mean(P$ok), mean(P$seg_eco), mean(P$seg_goodman), mean(P$seg_goodman_strict,na.rm=TRUE), sum(P$ok)))
cat(sprintf("corr(seg_goodman_strict, seg_eco)=%.4f\n", cor(P$seg_goodman_strict, P$seg_eco, use="complete.obs")))
