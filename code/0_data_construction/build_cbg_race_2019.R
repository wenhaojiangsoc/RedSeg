# Assemble data/cbg_race_2019.csv for compute_seg_reweighted.R.
# Race SHARES come from the RA covariate file (acs_cov_wj.RData$acs) -- the same CBG race
# composition that fed the original POI segregation measure, so the reweighting validates
# against the stored race_seg_msa_adj_mean. Total population (needed for the post-stratification
# adjustment factor and MSA-composition weighting) is pulled fresh from ACS 2019 5-yr B01003.
# Output cols: GEOID, tot_pop, white, black, hispanic, asian, other  (the last five are COUNTS
# = share*tot_pop, so the downstream count/tot_pop recovers the original shares).
suppressMessages({library(tidyverse); library(tidycensus)})
census_api_key(Sys.getenv("CENSUS_API_KEY"), install=FALSE)
RA <- "data/acs_cov_wj.RData"  # CBG ACS covariate file (authors)
e<-new.env(); load(RA, envir=e)
shares <- e$acs %>% transmute(GEOID=str_pad(as.character(GEOID),12,"left","0"),
  sh_white=cbg_white, sh_black=cbg_black, sh_hispanic=cbg_hispanic, sh_asian=cbg_asian, sh_other=cbg_other) %>%
  filter(!is.na(GEOID), GEOID!="NA")
# state FIPS for the 50 states + DC (exclude territories)
states <- unique(tidycensus::fips_codes$state_code)
states <- states[as.integer(states) <= 56]
pop <- map_dfr(states, function(s){
  message("pop: state ", s)
  tryCatch(get_acs(geography="block group", variables="B01003_001", state=s, year=2019, survey="acs5") %>%
             transmute(GEOID, tot_pop=estimate),
           error=function(err){message("  skip ",s,": ",conditionMessage(err)); tibble()})
})
cat(sprintf("pop CBGs: %d | share CBGs: %d\n", nrow(pop), nrow(shares)))
out <- pop %>% inner_join(shares, by="GEOID") %>% filter(tot_pop>0) %>%
  transmute(GEOID, tot_pop,
    white=sh_white*tot_pop, black=sh_black*tot_pop, hispanic=sh_hispanic*tot_pop,
    asian=sh_asian*tot_pop, other=sh_other*tot_pop)
cat(sprintf("matched CBGs written: %d (%.1f%% of pop CBGs)\n", nrow(out), 100*nrow(out)/nrow(pop)))
write_csv(out, "data/cbg_race_2019.csv")
cat("Saved data/cbg_race_2019.csv\n")
