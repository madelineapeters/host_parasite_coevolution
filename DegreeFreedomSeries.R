library(tidyverse)
library(gridExtra)

setwd("~/Desktop/ExamplePro/")

colors = c(rgb(0.227, 0.525, 1),rgb(0.224,0,0.6),rgb(1,0.741,0),rgb(1,0,0.33))

paras = read.csv("inputMut_para.csv",header=FALSE)
###### One-patch, continuous time ######
i.list = c(15,16,17)

stoch.df = data.frame()
#Interpolated
for (i in i.list){
  
  file.list = Sys.glob(file.path(paste0("~/Desktop/ExamplePro/Data_",i),paste0("data_",i,"_*.csv")))
  file.listN = Sys.glob(file.path(paste0("~/Desktop/ExamplePro/Data_",i),paste0("dataN_",i,"_*.csv")))
  
  df.rep = data.frame()
  df.repN = data.frame()
  for (f in 1:length(file.listN)){
  #for (f in 1:200){
    df = read.csv(file.list[f],header=FALSE) %>% select(.,1:9)
    dfN = read.csv(file.listN[f],header=FALSE) %>% select(.,1:9)
    names(df) = c("time","H11","H21","P11","P21","H12","H22","P12","P22")
    names(dfN) = names(df)
    
    df = df %>% 
      mutate(.,pH = (H11+H12)/(H11+H12+H21+H22)) %>% 
      mutate(.,pP = (P11+P12)/(P11+P12+P21+P22)) %>% 
      mutate(pH1 = H11/(H11+H21)) %>% 
      mutate(.,pH2 = H12/(H12+H22)) %>% 
      mutate(.,pP1 = P11/(P11+P21)) %>% 
      mutate(.,pP2 = P12/(P12+P22)) %>% 
      select(.,time,pH,pP,pH1,pH2,pP1,pP2)
    
    fun.pH = approx(x=df$time,y=df$pH,xout=seq(0,300,5),method="linear",f=0)
    fun.pP = approx(x=df$time,y=df$pP,xout=seq(0,300,5),method="linear",f=0)
    fun.pH1 = approx(x=df$time,y=df$pH1,xout=seq(0,300,5),method="linear",f=0)
    fun.pH2 = approx(x=df$time,y=df$pH2,xout=seq(0,300,5),method="linear",f=0)
    fun.pP1 = approx(x=df$time,y=df$pP1,xout=seq(0,300,5),method="linear",f=0)
    fun.pP2 = approx(x=df$time,y=df$pP2,xout=seq(0,300,5),method="linear",f=0)
    
    df.inter = as.data.frame(matrix(c(fun.pH$x,fun.pH$y,fun.pP$y,fun.pH1$y,fun.pH2$y,fun.pP1$y,fun.pP2$y),ncol=7,byrow=FALSE))
    names(df.inter) = c("time","pH","pP","pH1","pH2","pP1","pP2")
    
    dfN = dfN %>% 
      mutate(.,pH = (H11+H12)/(H11+H12+H21+H22)) %>% 
      mutate(.,pP = (P11+P12)/(P11+P12+P21+P22)) %>% 
      mutate(pH1 = H11/(H11+H21)) %>% 
      mutate(.,pH2 = H12/(H12+H22)) %>% 
      mutate(.,pP1 = P11/(P11+P21)) %>% 
      mutate(.,pP2 = P12/(P12+P22)) %>% 
      select(.,time,pH,pP,pH1,pH2,pP1,pP2)
    
    fun.pH = approx(x=dfN$time,y=dfN$pH,xout=seq(0,300,5),method="linear",f=0)
    fun.pP = approx(x=dfN$time,y=dfN$pP,xout=seq(0,300,5),method="linear",f=0)
    fun.pH1 = approx(x=dfN$time,y=dfN$pH1,xout=seq(0,300,5),method="linear",f=0)
    fun.pH2 = approx(x=dfN$time,y=dfN$pH2,xout=seq(0,300,5),method="linear",f=0)
    fun.pP1 = approx(x=dfN$time,y=dfN$pP1,xout=seq(0,300,5),method="linear",f=0)
    fun.pP2 = approx(x=dfN$time,y=dfN$pP2,xout=seq(0,300,5),method="linear",f=0)
    
    dfN.inter = as.data.frame(matrix(c(fun.pH$x,fun.pH$y,fun.pP$y,fun.pH1$y,fun.pH2$y,fun.pP1$y,fun.pP2$y),ncol=7,byrow=FALSE))
    names(dfN.inter) = c("time","pH","pP","pH1","pH2","pP1","pP2")
    
    df.inter$rep = f
    dfN.inter$rep = f
    df.rep = bind_rows(df.rep,df.inter)
    df.repN = bind_rows(df.repN,dfN.inter)
    
  }
  
  df.rep$type = "Coevolutionary"
  df.repN$type = "Neutral"
  df.tmp = bind_rows(df.rep,df.repN)
  df.tmp = df.tmp %>%
    mutate(.,HetP1 = 2*pP1*(1-pP1)) %>% 
    mutate(.,HetP2 = 2*pP2*(1-pP2)) %>% 
    mutate(.,HetP12 = 1 - (pP1*pP2+(1-pP1)*(1-pP2))) %>% 
    mutate(.,HetP = 2*pP*(1-pP)) %>% 
    mutate(.,HetH1 = 2*pH1*(1-pH1)) %>% 
    mutate(.,HetH2 = 2*pH2*(1-pH2)) %>% 
    mutate(.,HetH12 = 1 - (pH1*pH2+(1-pH1)*(1-pH2))) %>% 
    mutate(.,HetH = 2*pH*(1-pH))
  
  df.tmp$i = as.factor(i)
  stoch.df = bind_rows(stoch.df,df.tmp)
}

stoch.df$i = stoch.df$i %>% as.character %>% as.integer
stoch.df$i = factor(stoch.df$i,levels=c(15,16,17))
stoch.df = na.omit(stoch.df) %>% filter(.,time<300)

stoch.df = transform(stoch.df, b = ifelse(i==15,0.1,ifelse(i==16,0.2,0.3)))
stoch.df = transform(stoch.df, state = ifelse(pH>0&pH<1&pP>0&pP<1,0,ifelse(pH>0.45&pH<0.55&pP>0&pP<1,1,ifelse(pH>0&pH<1&HetP==0,2,ifelse(HetH==0&pP>0&pP<1,3,4)))))

stoch.single = stoch.df %>% group_by(i,time,type) %>% summarize(meanHetH=mean(HetH1),sdHetH=sd(HetH1))

one.cont.series = ggplot()+
  geom_ribbon(data=stoch.single %>% filter(time<=150),aes(x=time,ymin=meanHetH-2*sdHetH/sqrt(1000),ymax=meanHetH+2*sdHetH/sqrt(1000),fill=factor(i),linetype=type),alpha=0.25)+
  geom_line(data=stoch.single %>% filter(time<=150),aes(x=time,y=meanHetH,col=factor(i),linetype=type))+
  scale_color_manual(values=colors,labels=c("0.1","0.2","0.3"))+
  scale_fill_manual(values=colors,labels=c("0.1","0.2","0.3"))+
  labs(x="Host generations",y="Host heterozygosity",linetype=NULL,color=expression(i(lambda[`in`])),fill=expression(i(lambda[`in`])))+
  scale_linetype_manual(values=c("solid","dashed"),breaks=c("Neutral","Coevolutionary"))+
  theme_bw()+
  theme(strip.background = element_blank(),strip.text = element_text(size=14),panel.grid=element_blank(),axis.title=element_text(size=14),axis.text=element_text(size=12),legend.position=c(0.8,0.75),legend.text=element_text(size=12),plot.title = element_text(size=16,hjust=0.5))
one.cont.series

one.cont.grid = grid.arrange(one.cont.series,one.cont.facet,nrow=2,heights=c(2,1))

ggsave("~/Desktop/Mideo.lab2/Coevolution/Data sets/FIgures/OnePatch_Continuous_DegreeFreedomPanels.png",plot=one.cont.grid,width=6,height=8,units=c("in"),dpi=300)

ggsave("~/Desktop/Mideo.lab2/Coevolution/Data sets/FIgures/OnePatch_Continuous_DegreeFreedomSeries.png",width=7,height=5,units=c("in"),dpi=300)


###### Two-patch, continuous time ######
i.list = c(25,38,41,44)

stoch.df = data.frame()
#Interpolated
for (i in i.list){
  
  file.list = Sys.glob(file.path(paste0("~/Desktop/ExamplePro/Data_",i),paste0("data_",i,"_*.csv")))
  file.listN = Sys.glob(file.path(paste0("~/Desktop/ExamplePro/Data_",i),paste0("dataN_",i,"_*.csv")))
  
  df.rep = data.frame()
  df.repN = data.frame()
  for (f in 1:1000){
    #for (f in 1:200){
    df = read.csv(file.list[f],header=FALSE) %>% select(.,1:9)
    dfN = read.csv(file.listN[f],header=FALSE) %>% select(.,1:9)
    names(df) = c("time","H11","H21","P11","P21","H12","H22","P12","P22")
    names(dfN) = names(df)
    
    df = df %>% 
      mutate(.,pH = (H11+H12)/(H11+H12+H21+H22)) %>% 
      mutate(.,pP = (P11+P12)/(P11+P12+P21+P22)) %>% 
      mutate(pH1 = H11/(H11+H21)) %>% 
      mutate(.,pH2 = H12/(H12+H22)) %>% 
      mutate(.,pP1 = P11/(P11+P21)) %>% 
      mutate(.,pP2 = P12/(P12+P22)) %>% 
      select(.,time,pH,pP,pH1,pH2,pP1,pP2)
    
    fun.pH = approx(x=df$time,y=df$pH,xout=seq(0,150,5),method="linear",f=0)
    fun.pP = approx(x=df$time,y=df$pP,xout=seq(0,150,5),method="linear",f=0)
    fun.pH1 = approx(x=df$time,y=df$pH1,xout=seq(0,150,5),method="linear",f=0)
    fun.pH2 = approx(x=df$time,y=df$pH2,xout=seq(0,150,5),method="linear",f=0)
    fun.pP1 = approx(x=df$time,y=df$pP1,xout=seq(0,150,5),method="linear",f=0)
    fun.pP2 = approx(x=df$time,y=df$pP2,xout=seq(0,150,5),method="linear",f=0)
    
    df.inter = as.data.frame(matrix(c(fun.pH$x,fun.pH$y,fun.pP$y,fun.pH1$y,fun.pH2$y,fun.pP1$y,fun.pP2$y),ncol=7,byrow=FALSE))
    names(df.inter) = c("time","pH","pP","pH1","pH2","pP1","pP2")
    
    dfN = dfN %>% 
      mutate(.,pH = (H11+H12)/(H11+H12+H21+H22)) %>% 
      mutate(.,pP = (P11+P12)/(P11+P12+P21+P22)) %>% 
      mutate(pH1 = H11/(H11+H21)) %>% 
      mutate(.,pH2 = H12/(H12+H22)) %>% 
      mutate(.,pP1 = P11/(P11+P21)) %>% 
      mutate(.,pP2 = P12/(P12+P22)) %>% 
      select(.,time,pH,pP,pH1,pH2,pP1,pP2)
    
    fun.pH = approx(x=dfN$time,y=dfN$pH,xout=seq(0,150,5),method="linear",f=0)
    fun.pP = approx(x=dfN$time,y=dfN$pP,xout=seq(0,150,5),method="linear",f=0)
    fun.pH1 = approx(x=dfN$time,y=dfN$pH1,xout=seq(0,150,5),method="linear",f=0)
    fun.pH2 = approx(x=dfN$time,y=dfN$pH2,xout=seq(0,150,5),method="linear",f=0)
    fun.pP1 = approx(x=dfN$time,y=dfN$pP1,xout=seq(0,150,5),method="linear",f=0)
    fun.pP2 = approx(x=dfN$time,y=dfN$pP2,xout=seq(0,150,5),method="linear",f=0)
    
    dfN.inter = as.data.frame(matrix(c(fun.pH$x,fun.pH$y,fun.pP$y,fun.pH1$y,fun.pH2$y,fun.pP1$y,fun.pP2$y),ncol=7,byrow=FALSE))
    names(dfN.inter) = c("time","pH","pP","pH1","pH2","pP1","pP2")
    
    df.inter$rep = f
    dfN.inter$rep = f
    df.rep = bind_rows(df.rep,df.inter)
    df.repN = bind_rows(df.repN,dfN.inter)
    
  }
  
  df.rep$type = "Coevolutionary"
  df.repN$type = "Neutral"
  df.tmp = bind_rows(df.rep,df.repN)
  df.tmp = df.tmp %>%
    mutate(.,HetP1 = 2*pP1*(1-pP1)) %>% 
    mutate(.,HetP2 = 2*pP2*(1-pP2)) %>% 
    mutate(.,HetP12 = 1 - (pP1*pP2+(1-pP1)*(1-pP2))) %>% 
    mutate(.,HetP = 2*pP*(1-pP)) %>% 
    mutate(.,HetH1 = 2*pH1*(1-pH1)) %>% 
    mutate(.,HetH2 = 2*pH2*(1-pH2)) %>% 
    mutate(.,HetH12 = 1 - (pH1*pH2+(1-pH1)*(1-pH2))) %>% 
    mutate(.,HetH = 2*pH*(1-pH))
  
  df.tmp$i = as.factor(i)
  stoch.df = bind_rows(stoch.df,df.tmp)
}
stoch.df$i = factor(stoch.df$i,levels=c(25,38,41,44))
stoch.df = na.omit(stoch.df) %>% filter(.,time<150)

stoch.two = stoch.df %>% group_by(i,time,type) %>% summarize(meanHetH=mean(HetH),sdHetH=sd(HetH))

colors = c(rgb(0.227, 0.525, 1),rgb(0.224,0,0.6),rgb(1,0.741,0),rgb(1,0,0.33),"grey")

two.cont.series = ggplot()+
  geom_ribbon(data=stoch.two,aes(x=time,ymin=meanHetH-2*sdHetH/sqrt(1000),ymax=meanHetH+2*sdHetH/sqrt(1000),fill=factor(i),linetype=type),alpha=0.25)+
  geom_line(data=stoch.two,aes(x=time,y=meanHetH,col=factor(i),linetype=type))+
  scale_color_manual(values=colors,labels=c(
    expression(R*(lambda[`in`])*"= -0.1, "*i*(lambda[`in`])*"= 0.0975, "*lambda[mb]*" = 0.160"),
    expression(R*(lambda[`in`])*"= -0.3, "*i*(lambda[`in`])*"= 0.1, "*lambda[mb]*" = 0.150"),
    expression(R*(lambda[`in`])*"= -0.3, "*i*(lambda[`in`])*"= 0.1, "*lambda[mb]*" = 0.155"),
    expression(R*(lambda[`in`])*"= -0.3, "*i*(lambda[`in`])*"= 0.1, "*lambda[mb]*" = 0.160")
  ))+
  scale_fill_manual(values=colors,labels=c(
    expression(R*(lambda[`in`])*"= -0.1, "*i*(lambda[`in`])*"= 0.0975, "*lambda[mb]*" = 0.160"),
    expression(R*(lambda[`in`])*"= -0.3, "*i*(lambda[`in`])*"= 0.1, "*lambda[mb]*" = 0.150"),
    expression(R*(lambda[`in`])*"= -0.3, "*i*(lambda[`in`])*"= 0.1, "*lambda[mb]*" = 0.155"),
    expression(R*(lambda[`in`])*"= -0.3, "*i*(lambda[`in`])*"= 0.1, "*lambda[mb]*" = 0.160")
  ))+
  labs(x="Host generations",y="Host heterozygosity",linetype=NULL,color=NULL,fill=NULL)+
  scale_linetype_manual(values=c("solid","dashed"),breaks=c("Neutral","Coevolutionary"))+
  theme_bw()+
  theme(strip.background = element_blank(),strip.text = element_text(size=14),panel.grid=element_blank(),axis.title=element_text(size=14),axis.text=element_text(size=12),legend.position=c(0.275,0.22),legend.text=element_text(size=11),plot.title = element_text(size=16,hjust=0.5),legend.background=element_blank())

two.cont.series

ggsave("~/Desktop/Mideo.lab2/Coevolution/Data sets/FIgures/TwoPatch_Continuous_DegreeFreedomSeries.png",plot=two.cont.series,width=7,height=5,units=c("in"),dpi=300)



###### Island-mainland, continuous time ######
i.list = c(47,49,52)

stoch.df2 = data.frame()
#Interpolated
for (i in i.list){
  
  file.list = Sys.glob(file.path(paste0("~/Desktop/ExamplePro/Data_Island_",i),paste0("data_island_",i,"_*.csv")))
  file.listN = Sys.glob(file.path(paste0("~/Desktop/ExamplePro/Data_Island_",i),paste0("dataN_island_",i,"_*.csv")))
  
  df.rep = data.frame()
  df.repN = data.frame()
  for (f in 1:length(file.listN)){
    #for (f in 1:200){
    df = read.csv(file.list[f],header=FALSE) %>% select(.,1:9)
    dfN = read.csv(file.listN[f],header=FALSE) %>% select(.,1:9)
    names(df) = c("time","H11","H21","P11","P21","H12","H22","P12","P22")
    names(dfN) = names(df)
    
    df = df %>% 
      mutate(.,pH = (H11+H12)/(H11+H12+H21+H22)) %>% 
      mutate(.,pP = (P11+P12)/(P11+P12+P21+P22)) %>% 
      mutate(pH1 = H11/(H11+H21)) %>% 
      mutate(.,pH2 = H12/(H12+H22)) %>% 
      mutate(.,pP1 = P11/(P11+P21)) %>% 
      mutate(.,pP2 = P12/(P12+P22)) %>% 
      select(.,time,pH,pP,pH1,pH2,pP1,pP2)
    
    fun.pH = approx(x=df$time,y=df$pH,xout=seq(0,600,5),method="linear",f=0)
    fun.pP = approx(x=df$time,y=df$pP,xout=seq(0,600,5),method="linear",f=0)
    fun.pH1 = approx(x=df$time,y=df$pH1,xout=seq(0,600,5),method="linear",f=0)
    fun.pH2 = approx(x=df$time,y=df$pH2,xout=seq(0,600,5),method="linear",f=0)
    fun.pP1 = approx(x=df$time,y=df$pP1,xout=seq(0,600,5),method="linear",f=0)
    fun.pP2 = approx(x=df$time,y=df$pP2,xout=seq(0,600,5),method="linear",f=0)
    
    df.inter = as.data.frame(matrix(c(fun.pH$x,fun.pH$y,fun.pP$y,fun.pH1$y,fun.pH2$y,fun.pP1$y,fun.pP2$y),ncol=7,byrow=FALSE))
    names(df.inter) = c("time","pH","pP","pH1","pH2","pP1","pP2")
    
    dfN = dfN %>% 
      mutate(.,pH = (H11+H12)/(H11+H12+H21+H22)) %>% 
      mutate(.,pP = (P11+P12)/(P11+P12+P21+P22)) %>% 
      mutate(pH1 = H11/(H11+H21)) %>% 
      mutate(.,pH2 = H12/(H12+H22)) %>% 
      mutate(.,pP1 = P11/(P11+P21)) %>% 
      mutate(.,pP2 = P12/(P12+P22)) %>% 
      select(.,time,pH,pP,pH1,pH2,pP1,pP2)
    
    fun.pH = approx(x=dfN$time,y=dfN$pH,xout=seq(0,600,5),method="linear",f=0)
    fun.pP = approx(x=dfN$time,y=dfN$pP,xout=seq(0,600,5),method="linear",f=0)
    fun.pH1 = approx(x=dfN$time,y=dfN$pH1,xout=seq(0,600,5),method="linear",f=0)
    fun.pH2 = approx(x=dfN$time,y=dfN$pH2,xout=seq(0,600,5),method="linear",f=0)
    fun.pP1 = approx(x=dfN$time,y=dfN$pP1,xout=seq(0,600,5),method="linear",f=0)
    fun.pP2 = approx(x=dfN$time,y=dfN$pP2,xout=seq(0,600,5),method="linear",f=0)
    
    dfN.inter = as.data.frame(matrix(c(fun.pH$x,fun.pH$y,fun.pP$y,fun.pH1$y,fun.pH2$y,fun.pP1$y,fun.pP2$y),ncol=7,byrow=FALSE))
    names(dfN.inter) = c("time","pH","pP","pH1","pH2","pP1","pP2")
    
    df.inter$rep = f
    dfN.inter$rep = f
    df.rep = bind_rows(df.rep,df.inter)
    df.repN = bind_rows(df.repN,dfN.inter)
    
  }
  
  df.rep$type = "Coevolutionary"
  df.repN$type = "Neutral"
  df.tmp = bind_rows(df.rep,df.repN)
  df.tmp = df.tmp %>%
    mutate(.,HetP1 = 2*pP1*(1-pP1)) %>% 
    mutate(.,HetP2 = 2*pP2*(1-pP2)) %>% 
    mutate(.,HetP12 = 1 - (pP1*pP2+(1-pP1)*(1-pP2))) %>% 
    mutate(.,HetP = 2*pP*(1-pP)) %>% 
    mutate(.,HetH1 = 2*pH1*(1-pH1)) %>% 
    mutate(.,HetH2 = 2*pH2*(1-pH2)) %>% 
    mutate(.,HetH12 = 1 - (pH1*pH2+(1-pH1)*(1-pH2))) %>% 
    mutate(.,HetH = 2*pH*(1-pH))
  
  df.tmp$i = as.factor(i)
  stoch.df2 = bind_rows(stoch.df2,df.tmp)
}

stoch.df2$i = stoch.df2$i %>% as.character %>% as.integer
stoch.df2$i = factor(stoch.df2$i,levels=c(47,49,52))
stoch.df2 = na.omit(stoch.df2) %>% filter(.,time<600)

stoch.island = stoch.df2 %>% group_by(i,time,type) %>% summarize(meanHetH=mean(HetH1),sdHetH=sd(HetH1))

island.cont.series = ggplot()+
  #geom_vline(xintercept=145,col="gray")+
  geom_ribbon(data=stoch.island,aes(x=time,ymin=meanHetH-2*sdHetH/sqrt(1000),ymax=meanHetH+2*sdHetH/sqrt(1000),fill=factor(i),linetype=type),alpha=0.25)+
  geom_line(data=stoch.island,aes(x=time,y=meanHetH,col=factor(i),linetype=type))+
  scale_color_manual(values=colors,labels=c("-0.07,0.20","-0.06,0.22","-0.05,0.22"))+
  scale_fill_manual(values=colors,labels=c("-0.07,0.20","-0.06,0.22","-0.05,0.22"))+
  labs(x="Host generations",y="Host heterozygosity",linetype=NULL,color=expression(R(lambda[`in`])*","~i(lambda[`in`])),fill=expression(R(lambda[`in`])*","~i(lambda[`in`])))+
  scale_linetype_manual(values=c("solid","dashed"),breaks=c("Neutral","Coevolutionary"))+
  theme_bw()+
  theme(strip.background = element_blank(),strip.text = element_text(size=14),panel.grid=element_blank(),axis.title=element_text(size=14),axis.text=element_text(size=12),legend.position=c(0.7,0.7),legend.text=element_text(size=12),plot.title = element_text(size=16,hjust=0.5))
island.cont.series

island.cont.grid = grid.arrange(island.cont.series,island.cont.facet,nrow=2,heights=c(1,1))

ggsave("~/Desktop/Mideo.lab2/Coevolution/Data sets/FIgures/IslandMainland_Continuous_DegreeFreedomPanels.png",plot=island.cont.grid,width=6,height=9,units=c("in"),dpi=300)

ggsave("~/Desktop/Mideo.lab2/Coevolution/Data sets/FIgures/IslandMainland_Continuous_DegreeFreedomSeries.png",width=7,height=5,units=c("in"),dpi=300)
