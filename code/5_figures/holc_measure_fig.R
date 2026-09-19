# Figure explaining the integration-based place-segregation measure (redesign of the original poi.png).
# Two clean steps: (A) a location's visitor mix = flow-weighted blend of its visitors' home-CBG mixes;
# (B) segregation = summed distance of that mix from the metropolitan benchmark (integration).
suppressMessages({library(tidyverse); library(patchwork); library(grid); library(png)})
# Render the equation with REAL LaTeX (Computer Modern) and read it as a raster, so the figure's
# formula matches the manuscript typography rather than R's plotmath.
latex_grob<-function(math,w=2.6,dpi=900){
  td<-tempfile("eq"); dir.create(td)
  writeLines(sprintf("\\documentclass[border=1pt]{standalone}\\usepackage{amsmath}\\begin{document}$\\displaystyle %s$\\end{document}",math),
             file.path(td,"eq.tex"))
  system2("pdflatex",c("-interaction=nonstopmode","-output-directory",td,file.path(td,"eq.tex")),stdout=FALSE,stderr=FALSE)
  system2("pdftoppm",c("-r",dpi,"-png",file.path(td,"eq.pdf"),file.path(td,"eq")),stdout=FALSE,stderr=FALSE)
  rasterGrob(readPNG(file.path(td,"eq-1.png")),width=unit(w,"cm"),interpolate=TRUE)
}
eqA<-latex_grob("\\pi_{ik}=\\frac{\\sum_{j}\\gamma_{ji}\\,v_{jk}}{\\sum_{j}v_{jk}}")
GRP<-c("White","Black","Hispanic","Asian","Other")   # five groups, matching the paper (M = 5)
PAL<-c(White="#4E79A7",Black="#59A14F",Hispanic="#E15759",Asian="#EDC948",Other="#BAB0AC")
# illustrative shares (percent); each row sums to 100
cbg<-tribble(~cbg,~V,~White,~Black,~Hispanic,~Asian,~Other,
  "cbg 1",0.6,18,55,15,5,7, "cbg 2",0.3,65,10,13,7,5, "cbg 3",0.1,50,18,20,7,5)
poi<-c(White=35.3,Black=37.8,Hispanic=14.9,Asian=5.8,Other=6.2)  # = flow-weighted blend of the three CBGs
msa<-c(White=55,Black=12,Hispanic=18,Asian=8,Other=7)            # metropolitan benchmark (delta)
long<-function(v,lab)tibble(group=factor(GRP,levels=GRP),share=as.numeric(v[GRP]),who=lab)

## ---- Panel A: visitor mix is a flow-weighted blend ----
cbgL<-cbg%>%pivot_longer(all_of(GRP),names_to="group",values_to="share")%>%mutate(group=factor(group,levels=GRP),
   y=recode(cbg,`cbg 1`=3,`cbg 2`=2,`cbg 3`=1))
bar<-function(d,x0,x1,y,h)d%>%arrange(group)%>%mutate(cum=cumsum(share),xr=x0+(x1-x0)*cum/100,xl=lag(xr,default=x0),ymin=y-h,ymax=y+h)
cbgbars<-cbgL%>%group_by(cbg,y)%>%group_modify(~bar(.x,0,1.5,.y$y,.32))%>%ungroup()
poibars<-bar(long(poi,"poi")%>%rename(group=group),5.4,6.9,2,.42)
pA<-ggplot()+
  geom_rect(data=cbgbars,aes(xmin=xl,xmax=xr,ymin=ymin,ymax=ymax,fill=group),color="white",linewidth=.3)+
  geom_rect(data=poibars,aes(xmin=xl,xmax=xr,ymin=ymin,ymax=ymax,fill=group),color="white",linewidth=.35)+
  geom_text(data=cbg,aes(x=-0.08,y=recode(cbg,`cbg 1`=3,`cbg 2`=2,`cbg 3`=1),label=cbg),hjust=1,size=3.1,color="grey20")+
  geom_curve(data=cbg,aes(x=1.62,y=recode(cbg,`cbg 1`=3,`cbg 2`=2,`cbg 3`=1),xend=5.3,yend=2,linewidth=V),
             curvature=-0.18,color="grey55",arrow=arrow(length=unit(6,"pt"),type="closed"),alpha=.8,lineend="round")+
  scale_linewidth(range=c(.4,2.4),guide="none")+
  geom_text(data=cbg,aes(x=3.2,y=c(2.95,2.34,1.5),label=sprintf("v = %.1f",V)),size=2.8,color="grey35",fontface="italic")+
  annotate("text",x=6.15,y=2.72,label="a location's\nvisitor mix",size=3,color="grey20",lineheight=.9)+
  annotation_custom(eqA,xmin=2.35,xmax=4.35,ymin=-0.25,ymax=1.05)+
  scale_fill_manual(values=PAL,name=NULL)+coord_cartesian(xlim=c(-1,7.1),ylim=c(-.15,3.45),clip="off")+
  guides(fill=guide_legend(nrow=1,keywidth=unit(10,"pt"),keyheight=unit(10,"pt")))+
  labs(title="A   A location's visitor mix blends its visitors' home neighborhoods")+
  theme_void(base_size=11)+theme(plot.title=element_text(face="bold",size=10.5,hjust=0.5,margin=margin(b=10)),
    legend.position="top",legend.justification="left",legend.text=element_text(size=8.6),
    legend.margin=margin(t=4,0,0,14),legend.box.spacing=unit(6,"pt"),plot.margin=margin(3,8,4,20))
## ---- Panel B: distance from the metropolitan benchmark ----
cmp<-bind_rows(long(poi,"Location's visitor mix"),long(msa,"Metropolitan benchmark"))%>%
  mutate(who=factor(who,levels=c("Location's visitor mix","Metropolitan benchmark")))
dev<-tibble(group=factor(GRP,levels=GRP),lo=pmin(poi[GRP],msa[GRP]),hi=pmax(poi[GRP],msa[GRP]))
seg<-sum(abs(poi[GRP]-msa[GRP]))*(5/8)/100   # 5/8 normalizes the five-group index to [0,1]
pB<-ggplot()+
  geom_linerange(data=dev,aes(x=group,ymin=lo,ymax=hi),color="grey70",linewidth=3.2,alpha=.5)+
  geom_point(data=cmp,aes(group,share,color=who,shape=who),size=3.1,stroke=1)+
  scale_color_manual(values=c("Location's visitor mix"="#B0392B","Metropolitan benchmark"="grey25"),name=NULL)+
  scale_shape_manual(values=c("Location's visitor mix"=16,"Metropolitan benchmark"=1),name=NULL)+
  scale_y_continuous("Share of visitors / residents (%)",limits=c(0,70),expand=expansion(mult=c(0,.05)))+
  labs(x=NULL,title="B   Distance from the metropolitan benchmark")+
  annotate("text",x=2.7,y=60,hjust=.5,size=3.3,color="grey20",
    label=sprintf("Integration-based segregation = %.2f",seg))+
  theme_classic(base_size=11)+theme(plot.title=element_text(face="bold",size=10.5,hjust=0.5,margin=margin(b=6)),
    legend.position=c(.83,.93),legend.text=element_text(size=8.4),legend.key.height=unit(11,"pt"),
    axis.line=element_line(color="grey35",linewidth=.4),axis.title.y=element_text(size=8.6,color="grey30"),
    axis.text=element_text(color="grey25"),plot.margin=margin(6,10,4,6))
fig<-pA/pB+plot_layout(heights=c(1,1.05))
ggsave("figures/fig_measure.pdf",fig,width=7.4,height=6.2)
cat(sprintf("saved fig_measure.pdf | check: integration seg = %.2f\n",seg))
