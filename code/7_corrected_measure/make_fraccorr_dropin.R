# Corrected-fractional drop-in: puts the post-stratified fractional measure on the
# published seg_eco scale (x5/8 here; downstream scripts apply the second x5/8),
# stored under the seg_int_rw column name expected by the frac_variants scripts.
suppressMessages(library(dplyr))
a <- as_tibble(readRDS("data/poi_seg_frac_corrected.rds"))
out <- a %>% transmute(safegraph_place_id,
                       seg_int_rw  = seg_int_rw  * 5/8,
                       seg_int_unw = seg_int_unw * 5/8,
                       n_cbg, cov_rw)
saveRDS(out, "data/poi_seg_fraccorr_as_rw.rds")
b <- as_tibble(readRDS("results/poi_ei_bounds_frac.rds")) %>% transmute(safegraph_place_id, seg_eco)
m <- inner_join(out, b, by = "safegraph_place_id")
cat(sprintf("n matched %d | cor(corrected, uncorrected) = %.5f\n", nrow(m),
            cor(m$seg_int_rw, m$seg_eco, use = "complete.obs")))
cat(sprintf("means (seg_eco scale): corrected %.5f vs uncorrected %.5f | mean |diff| %.5f\n",
            mean(m$seg_int_rw, na.rm = TRUE), mean(m$seg_eco, na.rm = TRUE),
            mean(abs(m$seg_int_rw - m$seg_eco), na.rm = TRUE)))
# Also emit the corrected analog of results/poi_ei_bounds_frac.rds (seg_eco on its
# published scale, i.e. raw x 5/8) so scripts reading that file need no edits.
saveRDS(out %>% transmute(safegraph_place_id, seg_eco = seg_int_rw),
        "results/poi_ei_bounds_fraccorr.rds")
cat("Saved results/poi_ei_bounds_fraccorr.rds and data/poi_seg_fraccorr_as_rw.rds\n")
