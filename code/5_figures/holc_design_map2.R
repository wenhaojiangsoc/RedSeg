# =====================================================================================
# Design figure v2 (Los Angeles). Two equal-size square panels, side by side:
#   A: LA HOLC map (Mapping Inequality polygons) over a CartoDB Positron basemap,
#      with a black box marking the Panel-B zoom window.
#   B: zoomed border-design panel -- real B-C (solid) and B-D (dashed) borders (black),
#      the 1/2-mile counterfactual grid (grey), and ALL same-grade placebo grid borders
#      that enter the POI-level analysis (holc_cdb_poi.R) within the window (red).
#      NO buffer shading (removed at author's request).
# Grid + placebos reproduce the analysis exactly: grid anchored to 1/2-mi multiples of
#   EPSG:5070 (Conus Albers, hence slightly rotated vs north in LA); a cell edge is a valid
#   placebo if POIs sit on both sides within 1/4 mi (402.336 m) and all share one grade.
# Aesthetics (PAL, labels, colors, border linewidth 0.65, grid 0.22) match holc_design_map.R,
#   which is left untouched. Balance figure lives in holc_balance_fig.R.
# =====================================================================================
suppressMessages({library(tidyverse);library(sf);library(maptiles);library(tidyterra);library(patchwork)})
sf_use_s2(FALSE)
# run from the project root; # setwd("~/Library/CloudStorage/Dropbox/Demography-Rep-2026")

# ------------------------------- CONFIG ----------------------------------
CITY      <- "Los Angeles"
HOLC_JSON <- "data/mappinginequality.json"
WINDOW    <- c(xmin=-118.390,ymin=34.030,xmax=-118.315,ymax=34.085)  # zoom window (lon/lat), as in holc_design_map.R
PROJ_M    <- 3310                  # CA Albers (m) -- distance-correct
DISP      <- 3857                  # web mercator -- matches basemap tiles
GRID_M    <- 804.672               # 1/2-mile counterfactual grid (AHM)
BUFFER_M  <- 402.336               # 1/4-mile buffer (AHM "4th") -- placebo validity test
AN        <- 5070                  # analysis CRS (Conus Albers) -- grid anchored here (holc_cdb_poi.R)
PAL       <- c(A="#7E9E6E",B="#7FA3C0",C="#E2C368",D="#C98080")
LAB       <- c(A="A (best)",B="B (still desirable)",C="C (declining)",D="D (hazardous)")
COL_REAL  <- "black"; COL_PLACEBO <- "#C0392B"
LW_BORDER <- 0.85; LW_GRID <- 0.22; LW_POLY <- 0.25   # author's linewidths (real+placebo borders bumped 0.65->0.85)
BASE      <- 13.5                   # base font size (larger for Demography print)
OUT       <- "figures/fig0_design_map2"; OUT_W <- 12.5; OUT_H <- 7.4
# --------------------------------------------------------------------------

sq<-function(bb,pad=1){cx<-mean(bb[c("xmin","xmax")]);cy<-mean(bb[c("ymin","ymax")])
  hw<-max(bb["xmax"]-bb["xmin"],bb["ymax"]-bb["ymin"])/2*pad
  st_as_sfc(st_bbox(c(xmin=cx-hw,ymin=cy-hw,xmax=cx+hw,ymax=cy+hw),crs=st_crs(PROJ_M)))}
clip<-function(x,w) if(is.null(x)) NULL else suppressWarnings(st_intersection(x,w))

h  <- st_read(HOLC_JSON,quiet=TRUE)
la <- h%>%filter(grepl(CITY,city),grade%in%c("A","B","C","D"))%>%st_make_valid()%>%st_transform(PROJ_M)
Bp<-la%>%filter(grade=="B"); Cp<-la%>%filter(grade=="C"); Dp<-la%>%filter(grade=="D")
bl<-function(x,y){z<-suppressWarnings(st_intersection(st_geometry(st_boundary(x)),st_geometry(st_boundary(y))))
  z<-st_collection_extract(z,"LINESTRING"); if(length(z)==0) NULL else st_sf(geometry=z)}
BC<-bl(Bp,Cp); BD<-bl(Bp,Dp)

winB <- sq(st_bbox(st_transform(st_as_sfc(st_bbox(WINDOW,crs=4326)),PROJ_M)))   # square zoom window
# citywide window: zoom to the dense HOLC core (10-90% of polygon centroids), but keep winB inside
cc<-st_coordinates(st_centroid(la)); bw<-st_bbox(winB)
qb<-c(xmin=min(quantile(cc[,1],.10),bw["xmin"]),ymin=min(quantile(cc[,2],.10),bw["ymin"]),
      xmax=max(quantile(cc[,1],.90),bw["xmax"]),ymax=max(quantile(cc[,2],.90),bw["ymax"]))
winA <- sq(setNames(qb,c("xmin","ymin","xmax","ymax")),pad=1.02)

# ---- ANALYSIS 1/2-mile grid over the zoom window (anchored to 1/2-mi multiples of EPSG:5070) ----
winB_an<-st_transform(winB,AN); ba<-st_bbox(winB_an)
gx<-seq(floor(ba["xmin"]/GRID_M)*GRID_M, ceiling(ba["xmax"]/GRID_M)*GRID_M, by=GRID_M)
gy<-seq(floor(ba["ymin"]/GRID_M)*GRID_M, ceiling(ba["ymax"]/GRID_M)*GRID_M, by=GRID_M)
gridlines<-c(lapply(gx,function(x)st_linestring(rbind(c(x,min(gy)),c(x,max(gy))))),   # full spanning lines
             lapply(gy,function(y)st_linestring(rbind(c(min(gx),y),c(max(gx),y)))))
gridsf<-st_transform(st_segmentize(st_sfc(gridlines,crs=AN),100),PROJ_M)
# ---- ALL valid placebo cell edges, exactly as in the analysis (holc_cdb_poi.R lines 25-33):
#      POIs assigned to their nearest edge; edge valid if both sides present within 1/4 mi
#      and every assigned POI carries the same grade. Then keep edges in the zoom window.
pts<-readRDS("data/poi_bcd_pts.rds"); G<-GRID_M
kx<-round(pts$X/G)*G; dxv<-pts$X-kx; ky<-round(pts$Y/G)*G; dyh<-pts$Y-ky
vert<-abs(dxv)<=abs(dyh)
plac<-pts%>%mutate(bid=paste0(ifelse(vert,"V","H"),"_",ifelse(vert,round(X/G),round(Y/G)),"_",
    ifelse(vert,floor(Y/G),floor(X/G))),
  gside=ifelse(vert,sign(dxv),sign(dyh)),gdist=ifelse(vert,abs(dxv),abs(dyh)))%>%
  filter(gdist<=BUFFER_M,gside!=0)
validb<-plac%>%group_by(bid)%>%summarise(ns=n_distinct(gside),ng=n_distinct(grade),.groups="drop")%>%
  filter(ns==2,ng==1)%>%pull(bid)
mk<-function(s){p<-strsplit(s,"_")[[1]];a<-as.numeric(p[2]);b<-as.numeric(p[3])
  if(p[1]=="V")st_linestring(rbind(c(a*G,b*G),c(a*G,(b+1)*G))) else st_linestring(rbind(c(b*G,a*G),c((b+1)*G,a*G)))}
esf<-st_transform(st_segmentize(st_sfc(lapply(validb,mk),crs=AN),100),PROJ_M)
pick<-esf[lengths(st_intersects(esf,winB))>0]
cat(sprintf("valid analysis placebo borders: %d nationwide; %d in zoom window (all drawn)\n",
    length(validb),length(pick)))

# ================================ Panel A: LA citywide ================================
til<-get_tiles(st_transform(winA,DISP),provider="Esri.WorldGrayCanvas",zoom=11,crop=TRUE,cachedir=tempdir())  # CartoDB.Positron now watermarks "API KEY REQUIRED"
laA<-st_transform(clip(la,winA),DISP); winB_d<-st_transform(winB,DISP); bA<-st_bbox(st_transform(winA,DISP))
pA<-ggplot()+
  geom_spatraster_rgb(data=til,maxcell=6e6)+
  geom_sf(data=laA,aes(fill=grade),color="white",linewidth=.10,alpha=.55)+
  geom_sf(data=winB_d,fill=NA,color="#7B1113",linewidth=.55)+   # zoom-window box: dark red, thinner
  scale_fill_manual(values=PAL,name="HOLC grade",labels=LAB,limits=names(PAL),drop=FALSE)+
  coord_sf(crs=DISP,xlim=bA[c("xmin","xmax")],ylim=bA[c("ymin","ymax")],expand=FALSE)+
  labs(title=paste0("HOLC redlining map: ",CITY),tag="A")+
  theme_bw(base_size=BASE)+
  theme(plot.title=element_text(size=BASE+2,face="bold",hjust=.5),panel.grid=element_blank(),
        axis.text=element_blank(),axis.ticks=element_blank(),axis.title=element_blank(),
        aspect.ratio=1)

# ================================ Panel B: zoom border design ================================
laz<-clip(la,winB); BCz<-clip(BC,winB); BDz<-clip(BD,winB); bxB<-st_bbox(winB)
sb_x<-bxB["xmin"]+250; sb_y<-bxB["ymin"]+260
pB<-ggplot()+
  geom_sf(data=laz,aes(fill=grade),color="white",linewidth=LW_POLY,alpha=.75)+
  geom_sf(data=suppressWarnings(st_intersection(st_sf(geometry=gridsf),winB)),color="grey45",linewidth=LW_GRID,alpha=.8)+
  geom_sf(data=st_sf(geometry=pick),color=COL_PLACEBO,linewidth=LW_BORDER)+
  {if(!is.null(BCz)) geom_sf(data=BCz,color=COL_REAL,linewidth=LW_BORDER)}+
  {if(!is.null(BDz)) geom_sf(data=BDz,color=COL_REAL,linewidth=LW_BORDER,linetype="21")}+
  scale_fill_manual(values=PAL,name="HOLC grade",labels=LAB,limits=names(PAL),drop=FALSE,guide="none")+
  annotate("segment",x=sb_x,xend=sb_x+GRID_M,y=sb_y,yend=sb_y,linewidth=.7,color="black")+
  annotate("segment",x=sb_x,xend=sb_x,y=sb_y-45,yend=sb_y+45,linewidth=.5,color="black")+
  annotate("segment",x=sb_x+GRID_M,xend=sb_x+GRID_M,y=sb_y-45,yend=sb_y+45,linewidth=.5,color="black")+
  annotate("text",x=sb_x+GRID_M/2,y=sb_y+150,label="1/2 mile (grid)",size=3.4)+
  coord_sf(crs=PROJ_M,xlim=bxB[c("xmin","xmax")],ylim=bxB[c("ymin","ymax")],datum=NA,expand=FALSE)+
  labs(title="Boundary design (zoom)",tag="B")+
  theme_bw(base_size=BASE)+
  theme(plot.title=element_text(size=BASE+2,face="bold",hjust=.5),panel.grid=element_blank(),
        axis.text=element_blank(),axis.ticks=element_blank(),axis.title=element_blank(),
        aspect.ratio=1)
# line legend (borders + placebo) via dummy layer, shared at bottom -- B-D shown only if present
has_bd <- !is.null(BDz) && nrow(BDz)>0
lv  <- c("Real B-C border", if(has_bd) "Real B-D border", "Placebo grid border")
cv  <- c(`Real B-C border`=COL_REAL,`Real B-D border`=COL_REAL,`Placebo grid border`=COL_PLACEBO)[lv]
lty <- c(`Real B-C border`="solid",`Real B-D border`="21",`Placebo grid border`="solid")[lv]
leg<-tibble(x=NA_real_,y=NA_real_,what=factor(lv,levels=lv))
pB<-pB+geom_line(data=leg,aes(x,y,color=what,linetype=what),linewidth=.7,na.rm=TRUE)+
  scale_color_manual(values=cv,name=NULL)+scale_linetype_manual(values=lty,name=NULL)

fig<-(pA|pB)+plot_layout(widths=c(1,1),guides="collect")&
  theme(plot.tag=element_text(size=BASE+3,face="bold"),legend.position="bottom",legend.box="horizontal",
        legend.title=element_text(size=BASE),legend.text=element_text(size=BASE-1.5),
        legend.key.width=unit(1,"cm"))
ggsave(paste0(OUT,".pdf"),fig,width=OUT_W,height=OUT_H,useDingbats=FALSE)
ggsave(paste0(OUT,".png"),fig,width=OUT_W,height=OUT_H,dpi=300)
cat("Saved",paste0(OUT,".pdf/.png"),"\n")
