library(tidyverse)

path <- "./data/place-pulse-2.0/"

studies   <- read_tsv(paste0(path, "studies.tsv"))
locations <- read_tsv(paste0(path, "locations.tsv")) %>%
  rename(lat = `loc.0`, lon = `loc.1`)
places    <- read_tsv(paste0(path, "places.tsv")) %>%
  select(`_id`, city = place_name)
qscores   <- read_tsv(paste0(path, "qscores.tsv"))

# clean dimension labels from study names
dim_labels <- studies %>%
  select(`_id`, study_name) %>%
  mutate(dim = study_name %>%
           str_remove("_multicity") %>%
           str_replace_all(" ", "_") %>%
           str_replace("cleaner", "boring"))  # cleaner_multicity question = "more boring"

scores_long <- qscores %>%
  left_join(dim_labels,  by = c("study_id" = "_id")) %>%
  left_join(locations %>% select(`_id`, lat, lon),
            by = c("location_id" = "_id")) %>%
  left_join(places, by = c("place_id" = "_id"))

# one row per location, one column per dimension score
scores_wide <- scores_long %>%
  select(location_id, lat, lon, city, dim,
         score = trueskill.score,
         score_sd = `trueskill.stds.-1`,
         num_votes) %>%
  pivot_wider(
    names_from  = dim,
    values_from = c(score, score_sd, num_votes)
  )

us_cities <- c("Atlanta", "Boston", "Chicago", "Denver", "Houston",
               "Los Angeles", "Minneapolis", "New York", "Philadelphia",
               "Portland", "San Francisco", "Seattle", "Washington DC")

scores_us <- scores_wide %>%
  filter(city %in% us_cities)

write_csv(scores_us, "./data/placepulse_location_scores_us.csv")
cat("Done:", nrow(scores_us), "US locations written across", n_distinct(scores_us$city), "cities\n")

# ---- NYC safety map ----
library(ggmap)

nyc_pp <- scores_us %>%
  filter(city == "New York")

dim(nyc_pp[which(!is.na(nyc_pp$score_safer)),])

pad <- 0.06
bbox <- c(left   = min(nyc_pp$lon, na.rm = TRUE) - pad,
          bottom = min(nyc_pp$lat, na.rm = TRUE) - pad,
          right  = max(nyc_pp$lon, na.rm = TRUE) + pad,
          top    = max(nyc_pp$lat, na.rm = TRUE) + pad)

nyc_map <- get_map(location = bbox, source = "google", maptype = "terrain")

ggmap(nyc_map) +
  geom_point(data = nyc_pp,
             aes(x = lon, y = lat, color = score_safer),
             alpha = 0.7, size = 0.3) +
  scale_color_viridis_c(option = "plasma", direction = -1) +
  theme_void() +
  labs(title = "Perceived Safety, NYC (Place Pulse 2.0)",
       x = "", y = "",
       color = "Safety\nScore") +
  theme(legend.key.size = unit(0.5, "cm"))

ggsave("./results_race_new/nyc_placepulse_safety.png",
       dpi = 300, width = 79, height = 79, units = "mm")

# ---- Robustness: nearest Place Pulse score per restaurant ----
# Requires places_usa_2019 already in environment
library(sf)
library(fixest)
library(modelsummary)

dims <- c("score_safer", "score_wealthy", "score_more_beautiful",
          "score_livelier", "score_depressing", "score_boring")

pp_sf <- scores_us %>%
  filter(!is.na(lat), !is.na(lon)) %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326)

# test on NYC yelp restaurants first
nyc_yelp <- places_usa_2019 %>%
  filter(!is.na(yelp_rating), !is.na(longitude), !is.na(latitude))

poi_sf <- nyc_yelp %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326)

pp_sf_nyc <- scores_us %>%
  filter(city == "New York", !is.na(lat), !is.na(lon)) %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326)

idx <- st_nearest_feature(poi_sf, pp_sf_nyc)

places_pp <- nyc_yelp %>%
  bind_cols(st_drop_geometry(pp_sf_nyc)[idx, dims])

ctrl <- paste(dims, collapse = " + ")

model_pp <- list(
  "PRS:Ideal-Integration" = feols(
    race_seg_unadj_mean ~ redlining_CD + yelp_rating + yelp_price | poi_cbg + sub_category,
    data = places_pp
  ),
  "PRS:MSA-Benchmark" = feols(
    race_seg_msa_adj_mean ~ redlining_CD + yelp_rating + yelp_price | poi_cbg + sub_category,
    data = places_pp
  ),
  "PRS:Ideal-Integration+PP" = feols(
    as.formula(paste("race_seg_unadj_mean ~ redlining_CD + yelp_rating + yelp_price +", ctrl, "| poi_cbg + sub_category")),
    data = places_pp
  ),
  "PRS:MSA-Benchmark+PP" = feols(
    as.formula(paste("race_seg_msa_adj_mean ~ redlining_CD + yelp_rating + yelp_price +", ctrl, "| poi_cbg + sub_category")),
    data = places_pp
  )
)

modelsummary(
  model_pp,
  coef_map = c("redlining_CD" = "Redlining (C/D)",
               "yelp_rating"  = "Yelp Rating",
               "yelp_price"   = "Yelp Price"),
  fmt = 5,
  stars = TRUE,
  output = "./results_race_new/robustness_placepulse.html"
)
