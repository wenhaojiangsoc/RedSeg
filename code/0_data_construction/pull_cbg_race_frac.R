suppressMessages({library(tidyverse); library(tidycensus)})
census_api_key(Sys.getenv("CENSUS_API_KEY"), install=FALSE)
vars <- c(tot="B03002_001", white="B03002_003", black="B03002_004", asian="B03002_006", hisp="B03002_012")
sts <- unique(tidycensus::fips_codes$state_code); sts <- sts[as.integer(sts)<=56]
R <- map_dfr(sts, function(s) tryCatch(
  get_acs("block group", variables=vars, state=s, year=2019, survey="acs5", output="wide") %>%
    transmute(GEOID, tot_pop=totE, white=whiteE, black=blackE, hispanic=hispE, asian=asianE,
              other=pmax(0, totE-(whiteE+blackE+asianE+hispE))),
  error=function(e){message("skip ",s); tibble()}))
R <- R %>% filter(!is.na(tot_pop), tot_pop>0)
readr::write_csv(R, "data/cbg_race_frac_2019.csv")
mr <- with(R, pmax(white,black,hispanic,asian,other)/pmax(tot_pop,1))
cat(sprintf("wrote %d CBGs | mean max-race share=%.3f | share one-hot=%.3f | median=%.3f\n",
    nrow(R), mean(mr), mean(mr>0.999), median(mr)))
