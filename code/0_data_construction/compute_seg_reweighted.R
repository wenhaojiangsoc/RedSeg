# Recompute integration segregation with SafeGraph's CBG post-stratification correction.
# Method (SafeGraph, "Measuring and correcting sampling bias in SafeGraph Patterns"):
#   adjust_factor_g = (pop_g / total_pop) / (devices_g / total_devices)
#   i.e. visitors from over-represented home CBGs are down-weighted. Because the
#   segregation measure uses visitor SHARES, the two normalizing totals cancel, so
#   w_g = pop_g / devices_g gives identical shares (we use the full ratio anyway for
#   transparency and so weighted visit counts remain population-scaled).
#
# REQUIRED INPUTS (not in this repo -- set paths below):
#   FLOWS_PATH : POI x visitor-home-CBG table (the upstream `sf_cbg_poi_2019.rdata`
#                object `cbg_poi_2019`, or monthly Patterns visitor_home_cbgs long file)
#                cols: safegraph_place_id, visitor_home_cbgs, visitor_count [, month]
#   PANEL_PATH : home_panel_summary.csv (2019; monthly or annualized)
#                cols: census_block_group, number_devices_residing [, month]
#   CBGRACE_PATH: CBG-level race counts used in the original seg computation
#                cols: GEOID (12-digit CBG), tot_pop, and counts/shares for the same
#                five groups used originally (white, black, hispanic, asian, other)
# OUTPUT: ./data/poi_seg_reweighted.rds
#   safegraph_place_id, seg_int_rw (reweighted), seg_int_unw (unweighted recomputation,
#   validation against places_usa_2019$race_seg_msa_adj_mean), n_cbg, cov_rw
# To rerun the boundary DiD on the corrected measure: in holc_cdb_ebal.R (and the
# mediation scripts), left_join this rds by safegraph_place_id and set
# seg_int = seg_int_rw * 5/8 in place of race_seg_msa_adj_mean * 5/8.
suppressMessages({library(tidyverse); library(data.table)})
FLOWS_DIR   <- "./data/flows_2019"                  # monthly flows_*.csv.gz (safegraph_place_id,cbg,visitors)
PANEL_PATH  <- "./data/home_panel_summary_2019.csv" # 2019 annualized home panel (built from 52 weekly files)
CBGRACE_PATH<- "./data/cbg_race_2019.csv"           # GEOID,tot_pop,white,black,hispanic,asian,other (counts)
if(!dir.exists(FLOWS_DIR)||length(list.files(FLOWS_DIR,"\\.csv\\.gz$"))==0) stop("No flow files in ",FLOWS_DIR)
for(f in c(PANEL_PATH,CBGRACE_PATH)) if(!file.exists(f))
  stop("Missing input: ",f," -- see header for what this file must contain.")
# ---- 1. CBG adjustment factors ----
panel <- read_csv(PANEL_PATH,show_col_types=FALSE) %>%
  mutate(cbg=str_pad(as.character(census_block_group),12,"left","0")) %>%
  group_by(cbg) %>% summarise(devices=mean(number_devices_residing,na.rm=TRUE),.groups="drop")
cbg_race <- read_csv(CBGRACE_PATH,show_col_types=FALSE) %>%
  mutate(cbg=str_pad(as.character(GEOID),12,"left","0"))
adj <- cbg_race %>% select(cbg,tot_pop) %>% inner_join(panel,by="cbg") %>%
  filter(devices>0,tot_pop>0) %>%
  mutate(adjust_factor=(tot_pop/sum(tot_pop))/(devices/sum(devices)))
cat(sprintf("CBGs with adjustment factor: %d; factor quartiles: %s\n",
  nrow(adj),paste(round(quantile(adj$adjust_factor,c(.25,.5,.75)),2),collapse="/")))
# ---- 2. CBG lookup: race shares (+ tot_pop) and adjustment factor ----
grp <- c("sh_white","sh_black","sh_hispanic","sh_asian","sh_other")
shares <- cbg_race %>%
  transmute(cbg,tot_pop,across(c(white,black,hispanic,asian,other),~.x/tot_pop,.names="sh_{.col}"))
lk <- as.data.table(shares)[as.data.table(adj)[,.(cbg,adjust_factor)],on="cbg",nomatch=0]  # cbg,tot_pop,sh_*,adjust_factor
# ---- 3. Stream monthly flows -> POI-level weighted sums (memory-safe: one month at a time) ----
# Race shares depend only on home CBG (not month), so summing w*sh over all monthly rows and
# dividing by summed w gives the annual visitor composition -- no need to materialize all flows.
sumcols <- c("w","wa",paste0("uwsum_",grp),paste0("rwsum_",grp))
COMP_CACHE <- "./data/poi_seg_comp_intermediate.rds"   # cache the (slow) streamed POI sums
if(file.exists(COMP_CACHE)){
  comp <- readRDS(COMP_CACHE); cat("Loaded cached streamed sums from ",COMP_CACHE,"\n",sep="")
}else{
  ff <- list.files(FLOWS_DIR,"\\.csv\\.gz$",full.names=TRUE)
  cat(sprintf("Streaming %d monthly flow files...\n",length(ff)))
  acc <- NULL
  for(p in ff){
    m <- fread(cmd=paste("gzcat",shQuote(p)),colClasses=list(character="cbg"))
    m[,cbg:=str_pad(as.character(cbg),12,"left","0")]
    m <- lk[m,on="cbg",nomatch=0]                       # attach sh_* + adjust_factor
    m[,`:=`(w=visitors,wa=visitors*adjust_factor)]
    for(g in grp){set(m,j=paste0("uwsum_",g),value=m$w*m[[g]]); set(m,j=paste0("rwsum_",g),value=m$wa*m[[g]])}
    part <- m[,c(list(n_cbg=.N),lapply(.SD,sum,na.rm=TRUE)),by=safegraph_place_id,.SDcols=sumcols]
    acc <- if(is.null(acc)) part else
      rbindlist(list(acc,part))[,c(list(n_cbg=sum(n_cbg)),lapply(.SD,sum,na.rm=TRUE)),by=safegraph_place_id,.SDcols=sumcols]
    rm(m,part); gc(FALSE)
    cat("  ",basename(p),"done; POIs so far:",nrow(acc),"\n")
  }
  # POI-level unweighted & reweighted visitor shares (n_cbg = contributing POI x home-CBG-month rows)
  comp <- as_tibble(acc) %>% mutate(cov_rw=wa/w)
  for(g in grp){comp[[paste0("unw_",g)]] <- comp[[paste0("uwsum_",g)]]/comp$w
                comp[[paste0("rw_", g)]] <- comp[[paste0("rwsum_",g)]]/comp$wa}
  saveRDS(comp,COMP_CACHE); cat("Cached streamed sums to ",COMP_CACHE,"\n",sep="")
}
# ---- 4. Integration segregation: sum |visitor share - MSA share| over 5 groups ----
# (raw index in [0, 8/5]; the analysis scripts apply the *5/8 normalization)
load("./data/safegraph_2019_m.rdata")
msa_lk <- places_usa_2019 %>% distinct(safegraph_place_id,poi_cbg,msa,race_seg_msa_adj_mean)
msa_comp <- places_usa_2019 %>% distinct(poi_cbg,msa) %>%   # select only keys (places has its own tot_pop)
  mutate(cbg=str_pad(as.character(poi_cbg),12,"left","0")) %>%
  inner_join(shares,by="cbg") %>% group_by(msa) %>%
  summarise(across(all_of(grp),~weighted.mean(.x,tot_pop,na.rm=TRUE),.names="msa_{.col}"),.groups="drop")
out <- comp %>% inner_join(msa_lk,by="safegraph_place_id") %>% inner_join(msa_comp,by="msa")
out$seg_int_unw <- rowSums(sapply(grp,function(g) abs(out[[paste0("unw_",g)]]-out[[paste0("msa_",g)]])))
out$seg_int_rw  <- rowSums(sapply(grp,function(g) abs(out[[paste0("rw_", g)]]-out[[paste0("msa_",g)]])))
# ---- 5. Validation: unweighted recomputation should track the stored measure ----
cat(sprintf("Validation: cor(seg_int_unw, stored race_seg_msa_adj_mean) = %.3f (n=%d)\n",
  cor(out$seg_int_unw,out$race_seg_msa_adj_mean,use="complete.obs"),sum(complete.cases(out[,c("seg_int_unw","race_seg_msa_adj_mean")]))))
cat(sprintf("Reweighting shift: mean(seg_int_rw - seg_int_unw) = %+.5f; cor = %.3f\n",
  mean(out$seg_int_rw-out$seg_int_unw,na.rm=TRUE),cor(out$seg_int_rw,out$seg_int_unw,use="complete.obs")))
saveRDS(out%>%select(safegraph_place_id,seg_int_unw,seg_int_rw,n_cbg,cov_rw),
  "./data/poi_seg_reweighted.rds")
cat("Saved ./data/poi_seg_reweighted.rds\n")
