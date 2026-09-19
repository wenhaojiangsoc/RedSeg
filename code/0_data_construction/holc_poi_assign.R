# POI-POINT-LEVEL boundary assignment (Step 1 of the POI-level boundary DiD).
# For each POI: point-in-polygon HOLC grade, then assign to the nearest REAL B-vs-(C/D) polygon
# border within 1/4 mile (402.336 m), side = B (high) or C/D (low). Boundary = adjacent
# B<->C/D polygon pair (the fixed-effect unit). No tract-level assignment.
# Output: data/poi_real_assign.rds  (safegraph_place_id, boundary_id, low_grade, side, dist_m,
#         grade, tract11, sub_category, lvis, seg_int_unw, seg_int_rw)
suppressMessages({library(tidyverse); library(sf)})
sf_use_s2(FALSE)
# run from the project root; # setwd("~/Library/CloudStorage/Dropbox/Demography-Rep-2026")
CRS_M <- 5070; BUFFER_M <- 402.336

# ---- HOLC polygons (all cities) ----
h  <- st_read("data/mappinginequality.json",quiet=TRUE)
poly<- h%>%filter(grade%in%c("A","B","C","D"))%>%st_make_valid()%>%st_transform(CRS_M)%>%
  mutate(poly_id=row_number())%>%select(poly_id,grade)
Bp <- poly%>%filter(grade=="B"); CDp<- poly%>%filter(grade%in%c("C","D"))

# ---- real B<->C/D boundaries = adjacent B / (C or D) polygon pairs ----
cat("finding adjacent B-(C/D) polygon pairs...\n")
tou<- st_touches(Bp,CDp)                                   # sparse adjacency
pairs<- tibble(bi=rep(seq_len(nrow(Bp)),lengths(tou)), cj=unlist(tou))
cat("touching B-(C/D) pairs:",nrow(pairs),"\n")
bd_line<-function(i){
  ln<-tryCatch(suppressWarnings(st_intersection(st_geometry(Bp[pairs$bi[i],]),st_geometry(CDp[pairs$cj[i],]))),error=function(e)NULL)
  if(is.null(ln)||length(ln)==0) return(NULL)
  tp<-as.character(st_geometry_type(ln))
  if(any(tp=="GEOMETRYCOLLECTION")) ln<-suppressWarnings(st_collection_extract(ln,"LINESTRING"))
  ln<-ln[grepl("LINESTRING",as.character(st_geometry_type(ln)))]                # drop point-touches
  if(length(ln)==0) return(NULL)
  st_sf(boundary_id=i, b_poly=Bp$poly_id[pairs$bi[i]], cd_poly=CDp$poly_id[pairs$cj[i]],
        low_grade=CDp$grade[pairs$cj[i]], geometry=st_union(ln))}
bl<- bind_rows(lapply(seq_len(nrow(pairs)),bd_line))
cat("real boundaries with a shared border:",nrow(bl),"\n")

# ---- POIs -> points, join reweighted seg, point-in-polygon grade ----
load("data/safegraph_2019_m.rdata")
rw<- readRDS("data/poi_seg_reweighted.rds")
pts<- places_usa_2019%>%transmute(safegraph_place_id,latitude,longitude,sub_category,
        tract11=substr(poi_cbg,1,11),lvis=log(raw_visit_counts+1))%>%
      inner_join(rw%>%select(safegraph_place_id,seg_int_unw,seg_int_rw),by="safegraph_place_id")%>%
      filter(!is.na(latitude),!is.na(longitude),!is.na(seg_int_rw),!is.na(lvis))
poi<- st_as_sf(pts,coords=c("longitude","latitude"),crs=4326)%>%st_transform(CRS_M)
cat("POIs with coords+seg:",nrow(poi),"\n")
poi<- st_join(poi,poly,join=st_within)%>%filter(grade%in%c("B","C","D"))   # grade via point-in-polygon
cat("POIs in B/C/D polygons:",nrow(poi),"\n")
# save full B/C/D point set (with projected coords) for the placebo-grid step
poi_xy<- poi%>%mutate(X=st_coordinates(.)[,1],Y=st_coordinates(.)[,2])%>%st_drop_geometry()%>%
  select(safegraph_place_id,X,Y,grade,tract11,sub_category,lvis,seg_int_unw,seg_int_rw)
saveRDS(poi_xy,"data/poi_bcd_pts.rds"); cat("saved data/poi_bcd_pts.rds (",nrow(poi_xy),"pts )\n")

# ---- assign each POI to nearest real boundary of its own polygon, within 1/4 mile ----
poi<- poi%>%mutate(rid=row_number())
poiByPoly<- split(poi$rid, poi$poly_id)
assign1<-function(k){                                   # per boundary: high POIs in b_poly, low POIs in cd_poly
  hi<- poiByPoly[[as.character(bl$b_poly[k])]]; lo<- poiByPoly[[as.character(bl$cd_poly[k])]]
  out<-list()
  for(grp in list(list(ids=hi,side="high"),list(ids=lo,side="low"))){
    if(length(grp$ids)==0) next
    d<- as.numeric(st_distance(poi[match(grp$ids,poi$rid),], bl[k,]))
    keep<- which(d<=BUFFER_M); if(!length(keep)) next
    out[[grp$side]]<- tibble(rid=grp$ids[keep], boundary_id=bl$boundary_id[k],
                             low_grade=bl$low_grade[k], side=grp$side, dist_m=d[keep])}
  bind_rows(out)}
cat("assigning POIs to",nrow(bl),"boundaries...\n")
asg<- bind_rows(lapply(seq_len(nrow(bl)),assign1))
asg<- asg%>%group_by(rid)%>%slice_min(dist_m,n=1,with_ties=FALSE)%>%ungroup()   # nearest boundary per POI
res<- poi%>%st_drop_geometry()%>%inner_join(asg,by="rid")%>%
  select(safegraph_place_id,boundary_id,low_grade,side,dist_m,grade,tract11,sub_category,lvis,seg_int_unw,seg_int_rw)
saveRDS(res,"data/poi_real_assign.rds")
cat(sprintf("\nASSIGNED: %d POIs across %d real boundaries (both sides present).\n",
    nrow(res), res%>%group_by(boundary_id)%>%filter(n_distinct(side)==2)%>%pull(boundary_id)%>%n_distinct()))
cat(sprintf("median POIs per boundary: %.0f | median dist: %.0f m\n",
    median(table(res$boundary_id)), median(res$dist_m)))
