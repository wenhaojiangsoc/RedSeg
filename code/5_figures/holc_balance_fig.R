# Standalone covariate-balance figure (split out of holc_design_map.R at the author's request).
# Aesthetics (theme_classic, colors, shapes) are copied verbatim from holc_design_map.R Panel B
# -- do not restyle. Reads the balance table that holc_design_map.R already writes.
suppressMessages({library(tidyverse)})
# run from the project root; # setwd("~/Library/CloudStorage/Dropbox/Demography-Rep-2026")
bal <- read.csv("results_race_new/design_balance_three_schemes.csv") %>%
  mutate(scheme=factor(scheme,levels=c("Unweighted","Probit IPW","Entropy balancing")))

pB<-ggplot(bal,aes(bal,fct_rev(cov),color=scheme,shape=scheme))+
  geom_vline(xintercept=.05,linetype=3,color="grey40")+
  geom_point(size=2.8,stroke=1)+
  scale_color_manual(values=c(`Unweighted`="grey25",`Probit IPW`="#8C8C8C",`Entropy balancing`="#C0392B"),name=NULL)+
  scale_shape_manual(values=c(`Unweighted`=1,`Probit IPW`=17,`Entropy balancing`=16),name=NULL)+
  labs(x="Absolute standardized difference",y=NULL,title="Covariate balance")+
  theme_classic(base_size=11)+
  theme(panel.grid.minor=element_blank(),legend.position=c(1,.78),legend.justification=c(1,0),
        legend.background=element_rect(fill = "transparent", color = NA),legend.margin=margin(2,6,2,4),
        legend.key = element_blank(),
        plot.title=element_text(size=11,face="bold",hjust=0.5),axis.text=element_text(color="black"))

ggsave("figures/fig_covariate_balance.pdf",pB,width=5.4,height=4.6,useDingbats=FALSE)
ggsave("figures/fig_covariate_balance.png",pB,width=5.4,height=4.6,dpi=300)
cat("Saved fig_covariate_balance.pdf/.png\n")
