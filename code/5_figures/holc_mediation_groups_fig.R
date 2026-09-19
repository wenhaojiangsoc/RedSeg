# Figure 6: mechanism decomposition by racial group -- two matched coefficient panels (Fig-5 style).
#  A. Share of the boundary effect mediated by each group's distance from the metro (+ overall).
#  B. Redlining's effect on the SIGNED distance from the metro, by group (direction: Black above, white below).
# Dual-width dot-whisker (thick=90%, thin=95%, no caps), centered titles, minimal theme -- matches fig_hetero.
suppressMessages({library(tidyverse); library(patchwork)})
INK<-"#26374a"; POINT<-"#b0392b"; AGG<-"#26374a"; GREY<-"grey62"
L<-readRDS("results/holc_mediation_groups.rds"); med<-L$med; dev<-L$dev
nm<-c(overall="Overall (all groups)",black="Black",white="White",hispanic="Hispanic",asian="Asian")

# Panel A: % of effect mediated, with CI from the Sobel SE of the indirect effect
A<-med%>%transmute(group,prop=100*prop_mediated,
    se=abs(indirect/z)/c_total*100,
    kind=ifelse(group=="overall","Overall","Group"),
    lab=factor(nm[group],levels=rev(nm)))
# Panel B: signed a-path (effect of redlining on dev_g = share_nbhd - share_MSA)
B<-dev%>%transmute(group,prop=est,se=se,kind="Group",
    lab=factor(nm[group],levels=rev(nm[c("black","white","hispanic","asian")])))

dw<-function(d,xlab,ttl,xfmt,vline=TRUE){
  g<-ggplot(d,aes(prop,lab))
  if(vline) g<-g+geom_vline(xintercept=0,linetype="22",color=GREY,linewidth=.4)
  g+geom_linerange(aes(xmin=prop-1.96*se,xmax=prop+1.96*se),color=INK,linewidth=.5)+
    geom_linerange(aes(xmin=prop-1.645*se,xmax=prop+1.645*se),color=INK,linewidth=1.5)+
    geom_point(aes(fill=kind,shape=kind),color="white",size=3,stroke=.8)+
    scale_fill_manual(values=c(Group=POINT,Overall=AGG),guide="none")+
    scale_shape_manual(values=c(Group=21,Overall=23),guide="none")+
    scale_x_continuous(labels=xfmt)+
    labs(x=xlab,y=NULL,title=ttl)+
    theme_classic(base_size=11)+
    theme(plot.title=element_text(hjust=.5,face="bold",size=11,margin=margin(b=9)),
      axis.line=element_line(color="grey35",linewidth=.4),axis.ticks=element_line(color="grey55",linewidth=.35),
      axis.title.x=element_text(size=9,color="grey30",margin=margin(t=7)),
      axis.text.y=element_text(size=9.5,color="grey12"),axis.text.x=element_text(size=8,color="grey40"),
      plot.margin=margin(8,12,6,8))}

pA<-dw(A,"Share of boundary effect mediated (%)","How much each group's composition mediates",
       scales::number_format(accuracy=1,suffix="%"),vline=TRUE)
pB<-dw(B,"Effect of redlining on distance from metro","How much redlining displaces each group",
       scales::number_format(accuracy=.01),vline=TRUE)
fig<-(pA|pB)+plot_layout(widths=c(1,1))
ggsave("figures/fig_mediation_groups.pdf",fig,width=8.8,height=3.1,useDingbats=FALSE)
ggsave("figures/fig_mediation_groups.png",fig,width=8.8,height=3.1,dpi=300)
cat("saved fig_mediation_groups.pdf (Panel A: % mediated; Panel B: signed a-path)\n")
