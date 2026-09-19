# Figure 4: heterogeneity of the boundary DiD -- two matched coefficient panels, side by side.
#  A. Establishment type (7 groups).  B. Catchment size (N bins).
# Professional dot-whisker style: dual-width intervals (thick=90%, thin=95%, no end caps),
# centered panel titles, shared x-axis, minimal theme.
suppressMessages({library(tidyverse); library(patchwork)})
INK<-"#26374a"; POINT<-"#b0392b"; GREY<-"grey62"
NB<-6  # catchment bins (matched to 7 categories)
qa<-readRDS("results/holc_catch_multi.rds")%>%filter(nb==NB)%>%group_by(bin)%>%
  summarise(est=mean(est),se=mean(se),mc=mean(meancatch),.groups="drop")%>%mutate(lab=sprintf("%.0f km",mc/1000))
qa$lab<-factor(qa$lab,levels=qa$lab[order(qa$mc,decreasing=TRUE)])
cat7<-readRDS("results/holc_category7.rds")%>%group_by(catgroup)%>%summarise(est=mean(est),se=mean(se),.groups="drop")%>%
  arrange(est)%>%mutate(catgroup=factor(catgroup,levels=catgroup))
xr<-range(c(qa$est-1.96*qa$se,qa$est+1.96*qa$se,cat7$est-1.96*cat7$se,cat7$est+1.96*cat7$se)); xr<-xr+c(-1,1)*diff(xr)*.04
panel<-function(d,yv,ttl){ggplot(d,aes(est,.data[[yv]]))+
  geom_vline(xintercept=0,linetype="22",color=GREY,linewidth=.4)+
  geom_linerange(aes(xmin=est-1.96*se,xmax=est+1.96*se),color=INK,linewidth=.5)+
  geom_linerange(aes(xmin=est-1.645*se,xmax=est+1.645*se),color=INK,linewidth=1.5)+
  geom_point(shape=21,fill=POINT,color="white",size=3,stroke=.8)+
  scale_x_continuous(limits=xr,labels=scales::number_format(accuracy=.001))+
  labs(x="Boundary effect on integration segregation",y=NULL,title=ttl)+
  theme_classic(base_size=11)+
  theme(plot.title=element_text(hjust=.5,face="bold",size=11,margin=margin(b=9)),
    axis.line=element_line(color="grey35",linewidth=.4),axis.ticks=element_line(color="grey55",linewidth=.35),
    axis.title.x=element_text(size=9,color="grey30",margin=margin(t=7)),
    axis.text.y=element_text(size=9.5,color="grey12"),axis.text.x=element_text(size=8,color="grey40"),
    plot.margin=margin(8,12,6,8))}
pA<-panel(qa,"lab","By catchment size")+labs(y="Mean visitor travel")+
  theme(axis.title.y=element_text(size=8.5,color="grey40",angle=90,margin=margin(r=4)))
pB<-panel(cat7,"catgroup","By establishment type")
fig<-(pB|pA)+plot_layout(widths=c(1.35,1))
ggsave("figures/fig_hetero.pdf",fig,width=8.6,height=3.5)
cat("saved fig_hetero.pdf with",NB,"catchment bins / 7 categories\n")
