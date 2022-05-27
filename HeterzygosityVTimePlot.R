library(tidyverse)
library(gridExtra)

source("~/Desktop/Mideo.lab2/Coevolution/StochasticDF.R")

colors = c(rgb(0.227, 0.525, 1),rgb(0.224,0,0.6),rgb(1,0.741,0),rgb(1,0,0.33))

#### Two-patch ####
colors = c(rgb(0.224,0,0.6),rgb(1,0.741,0),rgb(1,0,0.33))

stoch.two = stoch.df %>% filter(.,i%in%c(1,2,3,5,6,7),time<=500) %>% group_by(i,time,model,type) %>% summarize(meanHetH=mean(HetH),sdHetH=sd(HetH))

strong.sel = ggplot()+
  geom_ribbon(data=stoch.two %>% filter(.,i%in%c(1,2,3)),aes(x=time,ymin=meanHetH-2*sdHetH/sqrt(1000),ymax=meanHetH+2*sdHetH/sqrt(1000),fill=factor(i),linetype=type),alpha=0.25)+
  geom_line(data=stoch.two %>% filter(.,i%in%c(1,2,3)),aes(x=time,y=meanHetH,col=factor(i),linetype=type))+
  scale_color_manual(values=colors,labels=c("Migration = 0.05","Migration = 0.01","Migration = 0.001"))+
  scale_fill_manual(values=colors,labels=c("Migration = 0.05","Migration = 0.01","Migration = 0.001"))+
  labs(x="Host generations",y="Host heterozygosity",linetype=NULL,color=NULL,tag="b)",fill=NULL)+
  scale_linetype_manual(values=c("solid","dashed"),breaks=c("Neutral","Coevolutionary"))+
  ggtitle("X = 0.7, Y = 0.1")+
  theme_bw()+
  theme(strip.background = element_blank(),strip.text = element_text(size=14),panel.grid=element_blank(),axis.title=element_text(size=14),axis.text=element_text(size=12),legend.position=c(0.7,0.75),legend.text=element_text(size=12),plot.title = element_text(size=16,hjust=0.5))


weak.sel = ggplot()+
  geom_ribbon(data=stoch.two %>% filter(.,i%in%c(5,6,7)),aes(x=time,ymin=meanHetH-2*sdHetH/sqrt(1000),ymax=meanHetH+2*sdHetH/sqrt(1000),fill=factor(i),linetype=type),alpha=0.25)+
  geom_line(data=stoch.two %>% filter(.,i%in%c(5,6,7)),aes(x=time,y=meanHetH,col=factor(i),linetype=type))+
  scale_color_manual(values=colors,labels=c("Migration = 0.05","Migration = 0.01","Migration = 0.001"))+
  scale_fill_manual(values=colors,labels=c("Migration = 0.05","Migration = 0.01","Migration = 0.001"))+
  labs(x="Host generations",y="Host heterozygosity",linetype=NULL,color=NULL,tag="a)",fill=NULL)+
  scale_linetype_manual(values=c("solid","dashed"),breaks=c("Neutral","Coevolutionary"))+
  ggtitle("X = 0.3, Y = 0.1")+
  theme_bw()+
  theme(strip.background = element_blank(),strip.text = element_text(size=14),panel.grid=element_blank(),axis.title=element_text(size=14),axis.text=element_text(size=12),legend.position="none",legend.text=element_text(size=12),plot.title = element_text(size=16,hjust=0.5))

g = grid.arrange(weak.sel,strong.sel,nrow=1)

ggsave("~/Desktop/Mideo.lab2/Coevolution/Data sets/FIgures/TwoPatch_Continuous_OGSeries.png",plot=g,width=10,height=5,units=c("in"),dpi=300)

#### Single patch ####
stoch.single = stoch.df %>% filter(.,i%in%c(4,8,54),time<=500) %>% group_by(i,time,model,type) %>% summarize(meanHetH=mean(HetH1),sdHetH=sd(HetH1))

single.cont.series = ggplot()+geom_ribbon(data=stoch.single,aes(x=time,ymin=meanHetH-2*sdHetH/sqrt(1000),ymax=meanHetH+2*sdHetH/sqrt(1000),fill=factor(i),linetype=type),alpha=0.25)+
  geom_line(data=stoch.single,aes(x=time,y=meanHetH,col=factor(i),linetype=type))+
  scale_color_manual(values=colors,labels=c("X=0.7, Y=0.1","X=0.3, Y=0.1","X=0.15, Y=0.1"))+
  scale_fill_manual(values=colors,labels=c("X=0.7, Y=0.1","X=0.3, Y=0.1","X=0.15, Y=0.1"))+
  labs(x="Host generations",y="Host heterozygosity",linetype=NULL,color=NULL,fill=NULL,tags="a)")+
  scale_linetype_manual(values=c("solid","dashed"),breaks=c("Neutral","Coevolutionary"))+
  theme_bw()+
  theme(strip.background = element_blank(),strip.text = element_text(size=14),panel.grid=element_blank(),axis.title=element_text(size=14),axis.text=element_text(size=12),legend.position=c(0.7,0.7),legend.text=element_text(size=12),plot.title = element_text(size=16,hjust=0.5))

stoch.single.tile = stoch.single %>% filter(.,i%in%c(4,8,54)) %>% select(.,-sdHetH) %>% filter(.,time==150) %>% pivot_wider(names_from="type",values_from=meanHetH)
stoch.single.tile$i = stoch.single.tile$i %>% as.character %>% as.numeric
stoch.single.tile$XY = 0.6
stoch.single.tile$XY[stoch.single.tile$i%in%c(5,6,7,8)] = 0.2
stoch.single.tile$XY[stoch.single.tile$i%in%c(54)] = 0.05

single.cont.tile = ggplot()+geom_tile(data=stoch.single.tile,aes(x=factor(`XY`),y=1,fill=Coevolutionary-Neutral))+labs(x="X-Y",y="",fill=expression(H[coev]-H[neut]))+
  theme_bw()+
  theme(strip.background = element_blank(), strip.text = element_text(size=14),panel.grid=element_blank(),axis.title=element_text(size=14),axis.text=element_text(size=12),legend.position="right",legend.text=element_text(size=10),plot.title = element_text(size=16,hjust=0.5),axis.text.y=element_blank(),axis.ticks.y=element_blank())+ labs(tag = "b)")

g = grid.arrange(single.cont.series,single.cont.tile,nrow=2,heights=c(2,1))

ggsave("~/Desktop/Mideo.lab2/Coevolution/Data sets/FIgures/OnePatch_Continuous_OGTile.png",g,width=7,height=8.5,units=c("in"),dpi=300)


#### Island-mainland ####
stoch.island = stoch.df2 %>% filter(.,i%in%c(1,2,3,5,6,7),time<=800) %>% group_by(i,time,model,type) %>% summarize(meanHetH=mean(HetH1),sdHetH=sd(HetH1))

colors = c(rgb(0.224,0,0.6),rgb(1,0.741,0),rgb(1,0,0.33))
strong.sel = ggplot()+
  geom_ribbon(data=stoch.island %>% filter(.,i%in%c(1,2,3)),aes(x=time,ymin=meanHetH-2*sdHetH/sqrt(1000),ymax=meanHetH+2*sdHetH/sqrt(1000),fill=factor(i),linetype=type),alpha=0.25)+
  geom_line(data=stoch.island %>% filter(.,i%in%c(1,2,3)),aes(x=time,y=meanHetH,col=factor(i),linetype=type))+
  scale_color_manual(values=colors,labels=c("Migration = 0.05","Migration = 0.01","Migration = 0.001"))+
  scale_fill_manual(values=colors,labels=c("Migration = 0.05","Migration = 0.01","Migration = 0.001"))+
  labs(x="Host generations",y="Host heterozygosity",linetype=NULL,color=NULL,tag="f)",fill=NULL)+
  scale_linetype_manual(values=c("solid","dashed"),breaks=c("Neutral","Coevolutionary"))+
  ggtitle("X = 0.7, Y = 0.1")+
  theme_bw()+
  theme(strip.background = element_blank(),strip.text = element_text(size=14),panel.grid=element_blank(),axis.title=element_text(size=14),axis.text=element_text(size=12),legend.position=c(0.7,0.75),legend.text=element_text(size=12),plot.title = element_text(size=16,hjust=0.5))


weak.sel = ggplot()+
  geom_ribbon(data=stoch.island %>% filter(.,i%in%c(5,6,7)),aes(x=time,ymin=meanHetH-2*sdHetH/sqrt(1000),ymax=meanHetH+2*sdHetH/sqrt(1000),fill=factor(i),linetype=type),alpha=0.25)+
  geom_line(data=stoch.island %>% filter(.,i%in%c(5,6,7)),aes(x=time,y=meanHetH,col=factor(i),linetype=type))+
  scale_color_manual(values=colors,labels=c("Migration = 0.05","Migration = 0.01","Migration = 0.001"))+
  scale_fill_manual(values=colors,labels=c("Migration = 0.05","Migration = 0.01","Migration = 0.001"))+
  labs(x="Host generations",y="Host heterozygosity",linetype=NULL,color=NULL,tag="e)",fill=NULL)+
  scale_linetype_manual(values=c("solid","dashed"),breaks=c("Neutral","Coevolutionary"))+
  ggtitle("X = 0.3, Y = 0.1")+
  theme_bw()+
  theme(strip.background = element_blank(),strip.text = element_text(size=14),panel.grid=element_blank(),axis.title=element_text(size=14),axis.text=element_text(size=12),legend.position="none",legend.text=element_text(size=12),plot.title = element_text(size=16,hjust=0.5))

g=grid.arrange(weak.sel,strong.sel,nrow=1)
ggsave("~/Desktop/Mideo.lab2/Coevolution/Data sets/FIgures/IslandMainland_Continuous_OGSeries.png",plot=g,width=10,height=5,units=c("in"),dpi=300)


stoch.island.tile = stoch.island %>% select(.,-sdHetH) %>% filter(.,time==700) %>% pivot_wider(names_from="type",values_from=meanHetH)
stoch.island.tile$i = stoch.island.tile$i %>% as.character %>% as.numeric
stoch.island.tile$XY = 0.6
stoch.island.tile$XY[stoch.island.tile$i%in%c(5,6,7)] = 0.2
stoch.island.tile$phi = c(0.05,0.01,0.001,0.05,0.01,0.001)

ggplot()+geom_tile(data=stoch.island.tile,aes(x=factor(phi),y=factor(`XY`),fill=Coevolutionary-Neutral))+labs(x=expression(phi),y="X-Y",fill=expression(H[coev]-H[neut]))+scale_fill_viridis_c(option="C")+
  theme_bw()+
  theme(strip.background = element_blank(),strip.text = element_text(size=14),panel.grid=element_blank(),axis.title=element_text(size=14),axis.text=element_text(size=12),legend.position="right",legend.text=element_text(size=10),plot.title = element_text(size=16,hjust=0.5))+ labs(tag = "g)")
ggsave("~/Desktop/Mideo.lab2/Coevolution/Data sets/FIgures/IslandMainland_Continuous_OGTile.png",width=7,height=5,units=c("in"),dpi=300)
