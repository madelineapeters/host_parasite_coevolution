library(tidyverse)

setwd("~/Desktop/Mideo.lab2/Coevolution/Data sets/")

i.list = c(1,2,3,4,5,6,7,8,9,10,11)
i.list = c(9,10,11)
i.list = c(54)
stoch.df = data.frame()
for (i in i.list){
  
  file.list = Sys.glob(file.path(paste0("~/Desktop/ExamplePro/Data_",i),paste0("data_",i,"_*.csv")))
  file.listN = Sys.glob(file.path(paste0("~/Desktop/ExamplePro/Data_",i),paste0("dataN_",i,"_*.csv")))
  
  df.rep = data.frame()
  df.repN = data.frame()
  for (f in 1:length(file.listN)){
    
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
    
    fun.pH = approx(x=df$time,y=df$pH,xout=seq(0,500,5),method="linear",f=0)
    fun.pP = approx(x=df$time,y=df$pP,xout=seq(0,500,5),method="linear",f=0)
    fun.pH1 = approx(x=df$time,y=df$pH1,xout=seq(0,500,5),method="linear",f=0)
    fun.pH2 = approx(x=df$time,y=df$pH2,xout=seq(0,500,5),method="linear",f=0)
    fun.pP1 = approx(x=df$time,y=df$pP1,xout=seq(0,500,5),method="linear",f=0)
    fun.pP2 = approx(x=df$time,y=df$pP2,xout=seq(0,500,5),method="linear",f=0)
    
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
    
    fun.pH = approx(x=dfN$time,y=dfN$pH,xout=seq(0,500,5),method="linear",f=0)
    fun.pP = approx(x=dfN$time,y=dfN$pP,xout=seq(0,500,5),method="linear",f=0)
    fun.pH1 = approx(x=dfN$time,y=dfN$pH1,xout=seq(0,500,5),method="linear",f=0)
    fun.pH2 = approx(x=dfN$time,y=dfN$pH2,xout=seq(0,500,5),method="linear",f=0)
    fun.pP1 = approx(x=dfN$time,y=dfN$pP1,xout=seq(0,500,5),method="linear",f=0)
    fun.pP2 = approx(x=dfN$time,y=dfN$pP2,xout=seq(0,500,5),method="linear",f=0)
    
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
    mutate(.,HetP1 = 1 - (pP1*pP1+(1-pP1)*(1-pP1))) %>% 
    mutate(.,HetP2 = 1 - (pP2*pP2+(1-pP2)*(1-pP2))) %>% 
    mutate(.,HetP12 = 1 - (pP1*pP2+(1-pP1)*(1-pP2))) %>% 
    mutate(.,HetP = 2*pP*(1-pP)) %>% 
    mutate(.,HetH1 = 1 - (pH1*pH1+(1-pH1)*(1-pH1))) %>% 
    mutate(.,HetH2 = 1 - (pH2*pH2+(1-pH2)*(1-pH2))) %>% 
    mutate(.,HetH12 = 1 - (pH1*pH2+(1-pH1)*(1-pH2))) %>% 
    mutate(.,HetH = 2*pH*(1-pH))
  
  df.tmp$i = as.factor(i)
  stoch.df = bind_rows(stoch.df,df.tmp)
}
stoch.df$i = stoch.df$i %>% as.character %>% as.integer
stoch.df$i = factor(stoch.df$i,levels=c(1,2,3,4,5,6,7,8,9,10,11,54))
stoch.df$model = "Two-patch"
stoch.df = na.omit(stoch.df) %>% filter(.,time<500)
write.csv(stoch.df,"TwoPatch_df.csv",row.names=FALSE)

i.list = c(12,13)
stoch.df2 = data.frame()
for (i in i.list){
  
  file.list = Sys.glob(file.path(paste0("~/Desktop/ExamplePro/Data_Island_",i),paste0("data_island_",i,"_*.csv")))
  file.listN = Sys.glob(file.path(paste0("~/Desktop/ExamplePro/Data_Island_",i),paste0("dataN_island_",i,"_*.csv")))
  
  df.rep = data.frame()
  df.repN = data.frame()
  for (f in 1:length(file.listN)){
    
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
    
    fun.pH = approx(x=df$time,y=df$pH,xout=seq(0,800,5),method="linear",f=0)
    fun.pP = approx(x=df$time,y=df$pP,xout=seq(0,800,5),method="linear",f=0)
    fun.pH1 = approx(x=df$time,y=df$pH1,xout=seq(0,800,5),method="linear",f=0)
    fun.pH2 = approx(x=df$time,y=df$pH2,xout=seq(0,800,5),method="linear",f=0)
    fun.pP1 = approx(x=df$time,y=df$pP1,xout=seq(0,800,5),method="linear",f=0)
    fun.pP2 = approx(x=df$time,y=df$pP2,xout=seq(0,800,5),method="linear",f=0)
    
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
    
    fun.pH = approx(x=dfN$time,y=dfN$pH,xout=seq(0,800,5),method="linear",f=0)
    fun.pP = approx(x=dfN$time,y=dfN$pP,xout=seq(0,800,5),method="linear",f=0)
    fun.pH1 = approx(x=dfN$time,y=dfN$pH1,xout=seq(0,800,5),method="linear",f=0)
    fun.pH2 = approx(x=dfN$time,y=dfN$pH2,xout=seq(0,800,5),method="linear",f=0)
    fun.pP1 = approx(x=dfN$time,y=dfN$pP1,xout=seq(0,800,5),method="linear",f=0)
    fun.pP2 = approx(x=dfN$time,y=dfN$pP2,xout=seq(0,800,5),method="linear",f=0)
    
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
    mutate(.,HetP1 = 1 - (pP1*pP1+(1-pP1)*(1-pP1))) %>% 
    mutate(.,HetP2 = 1 - (pP2*pP2+(1-pP2)*(1-pP2))) %>% 
    mutate(.,HetP12 = 1 - (pP1*pP2+(1-pP1)*(1-pP2))) %>% 
    mutate(.,HetP = 2*pP*(1-pP)) %>% 
    mutate(.,HetH1 = 1 - (pH1*pH1+(1-pH1)*(1-pH1))) %>% 
    mutate(.,HetH2 = 1 - (pH2*pH2+(1-pH2)*(1-pH2))) %>% 
    mutate(.,HetH12 = 1 - (pH1*pH2+(1-pH1)*(1-pH2))) %>% 
    mutate(.,HetH = 2*pH*(1-pH))
  
  df.tmp$i = as.factor(i)
  stoch.df2 = bind_rows(stoch.df2,df.tmp)
}
stoch.df2$i = factor(stoch.df2$i,levels=c(1:13))
stoch.df2 = na.omit(stoch.df2) %>% filter(.,time<800)
stoch.df2$model = "Island-mainland"
write.csv(stoch.df2,"IslandMainland_df.csv",row.names=FALSE)

stoch.df = read.csv("TwoPatch_df.csv")
stoch.df2 = read.csv("IslandMainland_df.csv")
