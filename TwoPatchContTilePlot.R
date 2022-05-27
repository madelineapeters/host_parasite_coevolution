library(gridExtra)
        
stoch.df$phi = 0.05
stoch.df$phi[stoch.df$i%in%c(2,6,10)] = 0.01
stoch.df$phi[stoch.df$i%in%c(3,7,11)] = 0.001
stoch.df$phi[stoch.df$i%in%c(4,8)] = 0

stoch.df$Z = 0.6
stoch.df$Z[stoch.df$i%in%c(5,6,7,8)] = 0.2
stoch.df$Z[stoch.df$i%in%c(9,10,11)] = 0.4

stoch.tile = stoch.df %>% filter(.,time==400,phi!=0) %>%  mutate(., Diff = abs(pH1-pH2)) %>% group_by(Z,phi) %>% summarize(meanDiff=mean(Diff)) %>% na.omit

diffTile = ggplot(stoch.tile)+
  geom_tile(aes(x=factor(phi),y=factor(Z),fill=log(meanDiff)))+
  labs(x=expression(phi),y="X-Y",fill=expression("ln|"~p[H1]-p[H2]~"|"))+
  theme_bw()+
  scale_fill_distiller(palette = "RdPu")+
  theme(strip.background = element_blank(),strip.text = element_text(size=14),panel.grid=element_blank(),axis.title=element_text(size=14),axis.text=element_text(size=12),legend.position="right",legend.text=element_text(size=10),plot.title = element_text(size=16,hjust=0.5))+ labs(tag = "c)")

phiDiffTile = ggplot(stoch.tile)+
  geom_tile(aes(x=factor(phi),y=factor(Z),fill=phi*meanDiff))+
  labs(x=expression(phi),y="X-Y",fill=expression(phi~"x"~"|"~p[H1]-p[H2]~"|"))+
  scale_fill_distiller(palette = "YlOrBr")+
  theme_bw()+
  theme(strip.background = element_blank(),strip.text = element_text(size=14),panel.grid=element_blank(),axis.title=element_text(size=14),axis.text=element_text(size=12),legend.position="right",legend.text=element_text(size=10),plot.title = element_text(size=16,hjust=0.5))+ labs(tag = "d)")

phiDiffTile

g = grid.arrange(diffTile,phiDiffTile,nrow=2)

ggsave("~/Desktop/Mideo.lab2/Coevolution/Data sets/FIgures/TwoPatch_Continuous_OGTile.png",plot=g,width=7,height=5,units=c("in"),dpi=300)

ggsave("~/Desktop/Mideo.lab2/Coevolution/Data sets/FIgures/TwoPatch_Continuous_PhiCoevFacet.png",plot=g,width=9.5,height=3.5,units=c("in"),dpi=300)

ggplot()+
  geom_smooth(data=stoch.df,aes(x=time,y=abs(pH1-pH2),col=factor(i)))

stoch.df %>% filter(.,time==400,phi!=0) %>% group_by(Z,phi) %>% summarize(COV=cov(pH1,pP1)) %>% na.omit
stoch.df %>% filter(.,time==400,phi!=0) %>% na.omit %>% group_by(Z,phi) %>% summarize(COVHPhome=mean(c(cov(pH1,pP1),cov(pH2,pP2))),COVHPaway=mean(c(cov(pH1,pP2),cov(pH2,pP1))),Fst=mean((0.5*(pH1+pH2)*(1-0.5*(pH1+pH2))-0.5*pH1*(1-pH1)-0.5*pH2*(1-pH2))/(0.5*(pH1+pH2)*(1-0.5*(pH1+pH2))),na.rm=T)) %>% mutate(.,LAapprox=-COVHPhome+COVHPaway)

ggplot()+
  geom_smooth(data=stoch.df %>% filter(.,time>=400,i%in%c(1,5,9)),aes(x=time,y=0.05*abs(pH1-pH2),col=factor(i)))+
  geom_smooth(data=stoch.df %>% filter(.,time>=400,i%in%c(2,6,10)),aes(x=time,y=0.01*abs(pH1-pH2),col=factor(i)))+
  geom_smooth(data=stoch.df %>% filter(.,time>=400,i%in%c(3,7,11)),aes(x=time,y=0.001*abs(pH1-pH2),col=factor(i)))

ggplot()+
  geom_smooth(data=stoch.df2 %>% filter(.,time>=700,i%in%c(1,5)),aes(x=time,y=0.05*abs(pH1-0.5),col=factor(i)))+
  geom_smooth(data=stoch.df2 %>% filter(.,time>=700,i%in%c(2,6)),aes(x=time,y=0.01*abs(pH1-0.5),col=factor(i)))+
  geom_smooth(data=stoch.df2 %>% filter(.,time>=700,i%in%c(3,7)),aes(x=time,y=0.001*abs(pH1-0.5),col=factor(i)))

ggplot()+geom_smooth(data=stoch.df2,aes(x=time,y=HetH1,col=factor(i)))

