library(tidyverse)
library(gridExtra)

setwd("~/Desktop/ExamplePro/")

colors = c(rgb(0.227, 0.525, 1),rgb(0.224,0,0.6),rgb(1,0.741,0),rgb(1,0,0.33))

paras = read.csv("inputMut_para.csv",header=FALSE)

###### One-patch, continuous ######
x=1
i.list = c(15,16,17)
df = read.csv(paste0("data_final_",x,".csv"),header=FALSE)
names(df) = c("gen","H11","H21","P11","P21","H12","H22","P12","P22")
df = df %>% transmute(.,
                   gen=gen,
                   pH1=(H11)/(H11+H21),
                   pH2=(H12)/(H12+H22),
                   pP1=(P11)/(P11+P21),
                   pP2=(P12)/(P12+P22),
                   pH=(H11)/(H11+H21),
                   pP=(P11)/(P11+P21),
                   HetH = 2*pH1*(1-pH1),
                   HetP = 2*pP1*(1-pP1)) %>% 
  filter(.,gen>295)

dfN = read.csv(paste0("dataN_final_",x,".csv"),header=FALSE)
names(dfN) = c("gen","H11","H21","P11","P21","H12","H22","P12",
              "P22")
dfN = dfN %>% transmute(.,
                   gen=gen,
                   pH1_N=(H11)/(H11+H21),
                   pH2_N=(H12)/(H12+H22),
                   pP1_N=(P11)/(P11+P21),
                   pP2_N=(P12)/(P12+P22),
                   pH_N=(H11)/(H11+H21),
                   pP_N=(P11)/(P11+P21),
                   HetH_N = 2*pH1_N*(1-pH1_N),
                   HetP_N = 2*pP1_N*(1-pP1_N)) %>% 
  filter(.,gen>295) %>% select(.,-gen)

df.joint = bind_cols(df,dfN)
df.joint$b = c(rep(0.1,1983),rep(0.2,1983),rep(0.3,1983))

df.summ = df.joint %>% group_by(b) %>% summarize(HetH = mean(HetH),
                                                 HetH_N=mean(HetH_N))

ggplot(df.joint)+geom_boxplot(aes(x=factor(b),y=HetH-HetH_N))

one.cont.facet=ggplot(df.summ)+geom_tile(aes(x=factor(b),y=1,fill=HetH-HetH_N))+
  labs(y=NULL,x=expression(italic(i)*(lambda[`in`])~"(= cycling rate, "%prop%boundary~ instability*")"),fill=expression(H[coev]-H[neut]))+
  geom_text(data=df.summ %>% filter(.,b!=0.3),aes(x=factor(b),y=1,label="-"),size=10)+
  geom_text(data=df.summ %>% filter(.,b==0.3),aes(x=factor(b),y=1,label="-"),col="white",size=10)+
  scale_fill_gradient(low=colors[1],high="white",guide=guide_colourbar(barheight = 5))+
  theme_bw()+
  theme(strip.background = element_blank(),strip.text = element_text(size=14),panel.grid=element_blank(),axis.title=element_text(size=14),axis.text=element_text(size=12),axis.text.y=element_blank(),axis.ticks.y=element_blank(),legend.position="right",legend.text=element_text(size=12),plot.title = element_text(size=16,hjust=0.5))
one.cont.facet
ggsave("~/Desktop/Mideo.lab2/Coevolution/Data sets/FIgures/OnePatch_Continuous_DegreeFreedomFacet.png",plot=one.cont.facet,width=7,height=2.5,units=c("in"),dpi=300)

 #Distinghishing between how many are 1) polymorphic in both species at the end 2) fixed for pathogen only, 3) fixed for host only, or 4) fixed for both. would be great
df$b = c(rep(0.1,992),rep(0.2,992),rep(0.3,991))
df1 = filter(df,pH1<1,pH1>0,pP1<1,pP1>0)
df2 = filter(df,pH1<1,pH1>0,HetP==0)
df3 = filter(df,HetH==0,pP>0,pP<1)
df4 = filter(df,HetH==0,HetP==0)
direct.sel = as.data.frame(matrix(nrow=4,ncol=4))
direct.sel$V1 = c("both.poly","path.fixed","host.fixed","both.fixed")
direct.sel[1,2:4] = table(df1$b)
direct.sel[2,2:4] = table(df2$b)
direct.sel[3,2:4] = table(df3$b)
direct.sel[4,2:4] = table(df4$b)
names(direct.sel) = c("state","b=0.1","b=0.2","b=0.3")
direct.sel


###### Two-patch, continuous ######
x=4
i.list = c(18:44)
df = read.csv(paste0("data_final_",x,".csv"),header=FALSE)
names(df) = c("gen","H11","H21","P11","P21","H12","H22","P12",
              "P22")
df = df %>% transmute(.,
                      gen=gen,
                      pH1=(H11)/(H11+H21),
                      pH2=(H12)/(H12+H22),
                      pP1=(P11)/(P11+P21),
                      pP2=(P12)/(P12+P22),
                      pH=(H11+H12)/(H11+H21+H12+H22),
                      pP=(P11+P12)/(P11+P21+P12+P22),
                      HetH = 2*pH*(1-pH),
                      HetP = 2*pP*(1-pP))

dfN = read.csv(paste0("dataN_final_",x,".csv"),header=FALSE)
names(dfN) = c("gen","H11","H21","P11","P21","H12","H22","P12",
               "P22")
dfN = dfN %>% transmute(.,
                        gen=gen,
                        pH1_N=(H11)/(H11+H21),
                        pH2_N=(H12)/(H12+H22),
                        pP1_N=(P11)/(P11+P21),
                        pP2_N=(P12)/(P12+P22),
                        pH_N=(H11+H12)/(H11+H21+H12+H22),
                        pP_N=(P11+P12)/(P11+P21+P12+P22),
                        HetH_N = 2*pH1_N*(1-pH1_N),
                        HetP_N = 2*pP1_N*(1-pP1_N)) %>% 
   select(.,-gen)

df.joint = bind_cols(df,dfN)

a.list = c(rep(-0.1,9009),rep(c(-0.2),each=9009),rep(-0.3,9009))
b.list = c(
  rep(rep(c(0.095,0.0975,0.1),each=1001),3),
  rep(rep(c(0.095,0.0975,0.1),each=1001),3),
  rep(rep(c(0.095,0.0975,0.1),each=1001),3)
  )
c.list = c(
  rep(c(0.15),3003),rep(c(0.155),3003),rep(c(0.16),3003),
  rep(c(0.15),3003),rep(c(0.155),3003),rep(c(0.16),3003),
  rep(c(0.15),3003),rep(c(0.155),3003),rep(c(0.16),3003)
)
df.joint$a = a.list
df.joint$b = b.list
df.joint$c = c.list

df.summ = df.joint %>% group_by(a,b,c) %>% summarize(HetH = mean(HetH),
                                                 HetH_N=mean(HetH_N)) %>% mutate(HDiff=HetH-HetH_N)

facetc1=ggplot(df.summ %>% filter(.,c==0.15))+geom_tile(aes(x=a,y=factor(b),fill=HetH-HetH_N))+
  geom_text(data=df.summ %>% filter(.,c==0.15,HDiff<0),aes(x=a,y=factor(b),label="-"),size=10)+
  labs(y=expression(i*(lambda[`in`])~"(= internal cycling)"),
       x="",fill=expression(H[coev]-H[neut]))+
  scale_fill_gradient(low=colors[1],high="white",guide=guide_colourbar(barheight = 10))+
  facet_grid(.~c)+
  ggtitle("")+
  theme_bw()+
  theme(strip.background = element_blank(),strip.text = element_text(size=14),panel.grid=element_blank(),axis.title=element_text(size=14),axis.text=element_text(size=12),legend.position="right",legend.text=element_text(size=12),plot.title = element_text(size=16,hjust=0.5))
facetc1

facetc2=ggplot(df.summ %>% filter(.,c==0.155))+geom_tile(aes(x=a,y=factor(b),fill=HetH-HetH_N))+
  geom_text(data=df.summ %>% filter(.,c==0.155,HDiff<0),aes(x=a,y=factor(b),label="-"),size=10)+
  labs(y="",
       x=expression(R*(lambda[`in`])~"(= internal stability)"),fill=expression(H[coev]-H[neut]))+
  scale_fill_gradient(low=colors[1],high="white",guide=guide_colourbar(barheight = 10))+
  ggtitle(expression(lambda[mb]~"(= matching boundary instability)"))+
  facet_grid(.~c)+
  theme_bw()+
  theme(strip.background = element_blank(),strip.text = element_text(size=14),panel.grid=element_blank(),axis.title=element_text(size=14),axis.text=element_text(size=12),legend.position="right",legend.text=element_text(size=12),plot.title = element_text(size=16,hjust=0.5))
facetc2

facetc3=ggplot(df.summ %>% filter(.,c==0.16))+geom_tile(aes(x=a,y=factor(b),fill=HetH-HetH_N))+
  geom_text(data=df.summ %>% filter(.,c==0.16,HDiff<0),aes(x=a,y=factor(b),label="-"),size=10)+
  geom_text(data=df.summ %>% filter(.,c==0.16,HDiff>0),aes(x=a,y=factor(b),label="+"),size=10)+
  labs(y="",
       x="",fill=expression(H[coev]-H[neut]))+
  scale_fill_gradient2(low=colors[1],mid="white",high=colors[3],midpoint=0,guide=guide_colourbar(barheight = 10))+
  facet_grid(.~c)+
  ggtitle("")+
  theme_bw()+
  theme(strip.background = element_blank(),strip.text = element_text(size=14),panel.grid=element_blank(),axis.title=element_text(size=14),axis.text=element_text(size=12),legend.position="right",legend.text=element_text(size=12),plot.title = element_text(size=16,hjust=0.5))
facetc3
two.cont.facet = grid.arrange(facetc1,facetc2,facetc3,nrow=1)

ggsave("~/Desktop/Mideo.lab2/Coevolution/Data sets/FIgures/TwoPatch_Continuous_DegreeFreedomFacet.png",plot=two.cont.facet,width=11,height=5,units=c("in"),dpi=300)
###### Island-mainland, continuous ######
x=3
i.list = c(45:53)
df = read.csv(paste0("data_final_",x,".csv"),header=FALSE)
names(df) = c("gen","H11","H21","P11","P21","H12","H22","P12",
              "P22")
df = df %>% transmute(.,
                      gen=gen,
                      pH1=(H11)/(H11+H21),
                      pH2=(H12)/(H12+H22),
                      pP1=(P11)/(P11+P21),
                      pP2=(P12)/(P12+P22),
                      pH=(H11+H12)/(H11+H21+H12+H22),
                      pP=(P11+P12)/(P11+P21+P12+P22),
                      HetH = 2*pH*(1-pH),
                      HetP = 2*pP*(1-pP)) #%>% filter(.,gen!=0)

dfN = read.csv(paste0("dataN_final_",x,".csv"),header=FALSE)
names(dfN) = c("gen","H11","H21","P11","P21","H12","H22","P12",
               "P22")
dfN = dfN %>% transmute(.,
                        gen=gen,
                        pH1_N=(H11)/(H11+H21),
                        pH2_N=(H12)/(H12+H22),
                        pP1_N=(P11)/(P11+P21),
                        pP2_N=(P12)/(P12+P22),
                        pH_N=(H11+H12)/(H11+H21+H12+H22),
                        pP_N=(P11+P12)/(P11+P21+P12+P22),
                        HetH_N = 2*pH1_N*(1-pH1_N),
                        HetP_N = 2*pP1_N*(1-pP1_N)) %>% 
  #filter(.,gen!=0) %>% 
  select(.,-gen)

df.joint = bind_cols(df,dfN)

a.list = c(
  rep(-0.07,1001),
  rep(-0.06,1983-1001),
  rep(-0.05,2928-1983),
  rep(-0.07,3881-2928),
  rep(-0.06,4850-3881),
  rep(-0.05,5797-4850),
  rep(-0.07,6739-5797),
  rep(-0.06,7701-6739),
  rep(-0.05,8702-7701)
)
b.list = c(
  rep(0.2,1001),
  rep(0.2,1983-1001),
  rep(0.2,2928-1983),
  rep(0.225,3881-2928),
  rep(0.225,4850-3881),
  rep(0.225,5797-4850),
  rep(0.25,6739-5797),
  rep(0.25,7701-6739),
  rep(0.25,8702-7701)
)

df.joint$a = a.list
df.joint$b = b.list

df.summ = df.joint %>% group_by(a,b) %>% summarize(HetH = mean(HetH),
                                                     HetH_N=mean(HetH_N)) %>% mutate(HDiff=HetH-HetH_N)

island.cont.facet = ggplot(df.summ)+geom_tile(aes(x=factor(a),y=factor(b),fill=HetH-HetH_N))+
  geom_text(data=df.summ %>% filter(.,HDiff>0),aes(x=factor(a),y=factor(b),label="+"),size=10)+
  labs(y=expression(i*(lambda[`in`])~"(= internal cycling, "%prop%boundary~ instability*")"),
       x=expression(R*(lambda[`in`])~"(= internal stability)"),fill=expression(H[coev]-H[neut]))+
  scale_fill_gradient(low=colors[1],high=colors[3],guide=guide_colourbar(barheight = 10))+
  ggtitle("")+
  theme_bw()+
  theme(strip.background = element_blank(),strip.text = element_text(size=14),panel.grid=element_blank(),axis.title=element_text(size=14),axis.text=element_text(size=12),axis.text.y=element_text(size=10),legend.position="right",legend.text=element_text(size=12),plot.title = element_text(size=16,hjust=0.5))
island.cont.facet
ggsave("~/Desktop/Mideo.lab2/Coevolution/Data sets/FIgures/IslandMainland_Continuous_DegreeFreedomFacet.png",width=7,height=5,units=c("in"),dpi=300)
###### One-patch, discrete ######
df.summ = as.data.frame(matrix(
  c(c(0.1,0.1,0.1,0.125,0.125,0.125,0.2,0.2,0.2),
    c(1.6,1.65,1.7,1.6,1.65,1.7,1.6,1.65,1.7),
    c(-0.114302,-0.111278,-0.108821,-0.163369,-0.156595,-0.143837,-0.215259,-0.217714,-0.210231)),byrow=FALSE,ncol=3))
names(df.summ)=c("a","b","HDiff")

ggplot()+geom_tile(data=df.summ,aes(x=factor(1+a),y=factor(b),fill=HDiff))+
  #geom_text(data=df.summ %>% filter(.,HDiff>0),aes(x=factor(1+a),y=factor(b),label="+"),size=10)+
  geom_text(data=df.summ %>% filter(.,HDiff<0),aes(x=factor(1+a),y=factor(b),label="-"),size=10)+
  labs(y=expression(lambda[nm]~"(= non-matching boundary instability)"),
       x=expression(lambda[`in`]~"(= internal instability)"),fill=expression(H[coev]-H[neut]))+
  scale_fill_gradient2(low=colors[1],mid="white",high=colors[3],midpoint=0,guide=guide_colourbar(barheight = 10))+
  ggtitle("")+
  theme_bw()+
  theme(strip.background = element_blank(),strip.text = element_text(size=14),panel.grid=element_blank(),axis.title=element_text(size=14),axis.text=element_text(size=12),legend.position="right",legend.text=element_text(size=12),plot.title = element_text(size=16,hjust=0.5))

ggsave("~/Desktop/Mideo.lab2/Coevolution/Data sets/FIgures/OnePatch_Discrete_DegreeFreedomFacet.png",width=7,height=5,units=c("in"),dpi=300)

###### Two-patch, discrete ######
df.summ = as.data.frame(matrix(
  c(c(0.1,0.1,0.1,0.125,0.125,0.125,0.15,0.15,0.15),
    c(1.6,1.65,1.7,1.6,1.65,1.7,1.6,1.65,1.7),
    c(-0.0290431,-0.0235713,-0.0242811,-0.0743003,-0.0726777,-0.0590879,-0.290422,-0.201128,-0.15208)),byrow=FALSE,ncol=3))

names(df.summ)=c("a","b","HDiff")

ggplot(df.summ)+geom_tile(aes(x=factor(1+a),y=factor(b),fill=HDiff))+
  labs(y=expression(lambda[nm]~"(= non-matching boundary instability)"),
       x=expression(lambda[`in`]~"(= internal instability)"),fill=expression(H[coev]-H[neut]))+
  geom_text(data=df.summ,aes(x=factor(a+1),y=factor(b),label="-"),size=10)+
  scale_fill_gradient2(low=colors[1],mid="white",high=colors[3],midpoint=0,guide=guide_colourbar(barheight = 10))+
  ggtitle("")+
  theme_bw()+
  theme(strip.background = element_blank(),strip.text = element_text(size=14),panel.grid=element_blank(),axis.title=element_text(size=14),axis.text=element_text(size=12),legend.position="right",legend.text=element_text(size=12),plot.title = element_text(size=16,hjust=0.5))

ggsave("~/Desktop/Mideo.lab2/Coevolution/Data sets/FIgures/TwoPatch_Discrete_DegreeFreedomFacet.png",width=7,height=5,units=c("in"),dpi=300)

###### Island-mainland, discrete ######
df.summ = as.data.frame(matrix(
  c(c(0.1,0.1,0.1,0.125,0.125,0.125,0.15,0.15,0.15),
    c(1.6,1.65,1.7,1.6,1.65,1.7,1.6,1.65,1.7),
    c(-0.0279453, -0.024746, -0.0324977, -0.0381983, -0.0315037,-0.0249818, -0.0577349, -0.0603313, -0.0551047)),byrow=FALSE,ncol=3))

names(df.summ)=c("a","b","HDiff")

ggplot(df.summ)+geom_tile(aes(x=factor(1+a),y=factor(b),fill=HDiff))+
  labs(y=expression(lambda[nm]~"(= non-matching boundary instability)"),
       x=expression(lambda[`in`]~"(= internal instability)"),fill=expression(H[coev]-H[neut]))+
  geom_text(data=df.summ,aes(x=factor(a+1),y=factor(b),label="-"),size=10)+
  scale_fill_gradient2(low=colors[1],mid="white",high=colors[3],midpoint=0,guide=guide_colourbar(barheight = 10))+
  ggtitle("")+
  theme_bw()+
  theme(strip.background = element_blank(),strip.text = element_text(size=14),panel.grid=element_blank(),axis.title=element_text(size=14),axis.text=element_text(size=12),legend.position="right",legend.text=element_text(size=12),plot.title = element_text(size=16,hjust=0.5))

ggsave("~/Desktop/Mideo.lab2/Coevolution/Data sets/FIgures/IslandMainland_Discrete_DegreeFreedomFacet.png",width=7,height=5,units=c("in"),dpi=300)
