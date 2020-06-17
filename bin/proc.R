setwd("/data/Research/2020nCov/coronatator/")
library(ggplot2)
library(cowplot)

sars2<-read.table("matrix/sars2.c.sg.matrix", header = T, row.names = 1)
x<-grep("PRJNA615032", rownames(sars2))
tmp<-rowSums(sars2[,2:ncol(sars2)])
tmpN<-sars2[,"N"]+1
sars2<-sars2/tmp
#sars2<-sars2/tmpN
sars2$sum<-tmp
sars2$type<-"vivo"
sars2$type[x]<-"vitro"
sars2$type <- as.factor(sars2$type)
sars2<-sars2[sars2$sum>=40,]
semn<-c("S","E","M","N")
acc<-c("ORF1","ORF3a","ORF6","ORF7a","ORF8")
n<-1
p<-list()
for(i in semn){
  p[[i]]<-ggplot(sars2, aes_string(x="type",y=i))+
    geom_boxplot()+
    geom_dotplot(binaxis='y', stackdir='center', dotsize=1)+
    theme_classic()
  n<-n+1 
}
plot_grid(plotlist = p, nrow = 1)
#other normalization methods
gorder<-c("S","ORF3a","E","M","ORF6","ORF7a","ORF8","N")
sars2.gorder <- sars2[,gorder]
tmp <-as.data.frame(t(cumsum(as.data.frame(t(sars2.gorder)))))
sars2.gorder <- sars2.gorder/tmp
#sars2 bp plot
sars2.bp <- read.table("matrix.new/sars2.c.sg.proc.sort.cls",col.names = c("X","Cnt","Ratio","Anno","Overlap","Prj","Srx"))
sars2.bp <- sars2.bp[sars2.bp$Ratio >0.001,]
sars2.bp <- sars2.bp[sars2.bp$Ratio <0.005,]
sars2.bp$Anno <-gsub("ORF","",sars2.bp$Anno)
p<-ggplot()+
    geom_point(data=sars2.bp,aes(x=X, y=Ratio))+
    geom_linerange(data=sars2.bp,aes(x=X,ymin=0,ymax=Ratio))+
    geom_text(data=sars2.bp,aes(x=X,y=Ratio, label=Anno),nudge_y=0.03, size=3,col="red")+
    geom_text(data=sars2.bp,aes(x=X,label=Overlap,y=-0.06, angle=270, hjust=0),size=2.5,col="blue")+
    ylim(-0.2,0.8)+
    theme_classic()

p<-p+theme(axis.title.x = element_blank(),
           axis.line.x = element_blank(),
           axis.ticks.x = element_blank(),
           axis.text.x = element_blank())
p<-p+geom_hline(yintercept = 0)

x.ticks<-c(1:30)*1000
my.df<-data.frame(xt=as.numeric(x.ticks),label=rep("",length(x.ticks)),ymin=-0.005)
my.df$label <- as.character(my.df$label)
tmp<-which(my.df$xt%%5000==0)
my.df$label[tmp] <- as.character(my.df$xt[tmp])
my.df$ymin[tmp]<--0.01
p<-p+geom_text(data=my.df, aes(x=xt, label=label ,y=-0.03),size=3)+
    geom_linerange(data=my.df, aes(x=xt, ymin=ymin, ymax=0))


###add Srx Prj support
tmp <- read.table("matrix.new/sars2.p.s")
sars2.sp <- tmp$V1
names(sars2.sp)<-tmp$V2
sars2.fdot<-read.table("matrix.new/sars2.c.sg.proc.sort.cls.fdot",col.names = c("X","Cnt","Ratio","Srx"))
sars2.fdot$Prj<-sars2.sp[sars2.fdot$Srx]
sars2.fdot <- sars2.fdot[sars2.fdot$Ratio <0.005,]
sars2.fdot <- sars2.fdot[sars2.fdot$Ratio >0.001,]
p+geom_dotplot(data=sars2.fdot, aes(x=X,fill=Prj,y=Cnt),stackdir="up",method = "histodot",
               binwidth = 100,stackgroups = T,binpositions = "all",dotsize = 1.5,
               position=position_nudge(y=sars2.fdot$Ratio+0.04))+
    theme_classic()+theme(axis.title.x = element_blank(),
                          axis.line.x = element_blank(),
                          axis.ticks.x = element_blank(),
                          axis.text.x = element_blank())


#sars
sars<-read.table("matrix/sars.c.sg.matrix", header = T, row.names = 1)
tmp<-rowSums(sars[,2:ncol(sars)])
tmpN<-sars[,"N"]+1
#sars<-sars/tmp
sars<-sars/tmpN
sars$sum<-tmp
sars$type<-"vitro"
sars$type <- as.factor(sars$type)
sars<-sars[sars$sum>=40,]
semn<-c("S","E","M","N")
acc<-c("ORF1","ORF3a","ORF3b","ORF6","ORF7a","ORF8a")
n<-1
p<-list()
for(i in semn){
  p[[i]]<-ggplot(sars, aes_string(x="type",y=i))+
    geom_boxplot()+
    geom_dotplot(binaxis='y', stackdir='center', dotsize=1)+
    theme_classic()
  n<-n+1 
}
plot_grid(plotlist = p, nrow = 1)
#sars bp plot

sars.bp <- read.table("matrix/sars.c.sg.proc.sort.cls",col.names = c("X","Cnt","Ratio","Anno","Overlap","Prj","Srx"))
sars.bp <- sars.bp[sars.bp$Ratio >0.001,]
sars.bp$Anno <-gsub("ORF","",sars.bp$Anno)
p<-ggplot()+
  geom_point(data=sars.bp,aes(x=X, y=Ratio))+
  geom_linerange(data=sars.bp,aes(x=X,ymin=0,ymax=Ratio))+
  geom_text(data=sars.bp,aes(x=X,y=Ratio, label=Anno),nudge_y=0.03, size=3,col="red")+
  geom_text(data=sars.bp,aes(x=X,label=Overlap,y=-0.06, angle=270, hjust=0),size=2.5,col="blue")+
  ylim(-0.2,0.8)+
  theme_classic()

p<-p+theme(axis.title.x = element_blank(),
           axis.line.x = element_blank(),
           axis.ticks.x = element_blank(),
           axis.text.x = element_blank())
p<-p+geom_hline(yintercept = 0)

x.ticks<-c(1:30)*1000
my.df<-data.frame(xt=as.numeric(x.ticks),label=rep("",length(x.ticks)),ymin=-0.005)
my.df$label <- as.character(my.df$label)
tmp<-which(my.df$xt%%5000==0)
my.df$label[tmp] <- as.character(my.df$xt[tmp])
my.df$ymin[tmp]<--0.01
p<-p+geom_text(data=my.df, aes(x=xt, label=label ,y=-0.03),size=3)+
  geom_linerange(data=my.df, aes(x=xt, ymin=ymin, ymax=0))


###add Srx Prj support
tmp <- read.table("matrix/sars.p.s")
sars.sp <- tmp$V1
names(sars.sp)<-tmp$V2
sars.fdot<-read.table("matrix/sars.c.sg.proc.sort.cls.fdot",col.names = c("X","Cnt","Ratio","Srx"))
sars.fdot$Prj<-sars.sp[sars.fdot$Srx]
sars.fdot <- sars.fdot[sars.fdot$Ratio >0.001,]

p+geom_dotplot(data=sars.fdot, aes(x=X,fill=Prj,y=Cnt),stackdir="up",method = "histodot",
               binwidth = 100,stackgroups = T,binpositions = "all",dotsize = 1.5,
               position=position_nudge(y=sars.fdot$Ratio+0.04))+
  theme_classic()+theme(axis.title.x = element_blank(),
                        axis.line.x = element_blank(),
                        axis.ticks.x = element_blank(),
                        axis.text.x = element_blank())


#mers
mers<-read.table("matrix/mers.c.sg.matrix", header = T, row.names = 1)
vivo<-c("PRJNA238265","PRJNA545350")
x<-grep("PRJNA238265|PRJNA545350",rownames(mers))
tmp<-rowSums(mers[,2:ncol(mers)])
tmpN <- mers[,"N"]+1
mers<-mers/tmpN
#mers<-mers/tmp
mers$sum<-tmp
mers$type<-"vitro"
mers$type[x]<-"vivo"
mers$type <- as.factor(mers$type)
mers<-mers[mers$sum>=40,]
semn<-c("S","E","M","N")
acc<-c("ORF1","ORF3","ORF4a","ORF5","ORF8b")
n<-1
p<-list()
for(i in acc){
  p[[i]]<-ggplot(mers, aes_string(x="type",y=i))+
    geom_boxplot()+
    geom_dotplot(binaxis='y', stackdir='center', dotsize=1)+
    theme_classic()
  n<-n+1 
}
plot_grid(plotlist = p, nrow = 1)
#other normalization methods
gorder<-c("S","ORF3a","E","M","ORF6","ORF7a","ORF8","N")
sars2.gorder <- sars2[,gorder]
tmp <-as.data.frame(t(cumsum(as.data.frame(t(sars2.gorder)))))
sars2.gorder <- sars2.gorder/tmp


#mers bp plot

mers.bp <- read.table("matrix/mers.c.sg.proc.sort.cls",col.names = c("X","Cnt","Ratio","Anno","Overlap","Prj","Srx"))
mers.bp <- mers.bp[mers.bp$Ratio >0.005,]
mers.bp$Anno <-gsub("ORF","",mers.bp$Anno)
p<-ggplot()+
  geom_point(data=mers.bp,aes(x=X, y=Ratio))+
  geom_linerange(data=mers.bp,aes(x=X,ymin=0,ymax=Ratio))+
  geom_text(data=mers.bp,aes(x=X,y=Ratio, label=Anno),nudge_y=0.03, size=3,col="red")+
  geom_text(data=mers.bp,aes(x=X,label=Overlap,y=-0.06, angle=270, hjust=0),size=2.5,col="blue")+
  ylim(-0.2,0.8)+
  theme_classic()

p<-p+theme(axis.title.x = element_blank(),
           axis.line.x = element_blank(),
           axis.ticks.x = element_blank(),
           axis.text.x = element_blank())
p<-p+geom_hline(yintercept = 0)

x.ticks<-c(1:30)*1000
my.df<-data.frame(xt=as.numeric(x.ticks),label=rep("",length(x.ticks)),ymin=-0.005)
my.df$label <- as.character(my.df$label)
tmp<-which(my.df$xt%%5000==0)
my.df$label[tmp] <- as.character(my.df$xt[tmp])
my.df$ymin[tmp]<--0.01
p<-p+geom_text(data=my.df, aes(x=xt, label=label ,y=-0.03),size=3)+
  geom_linerange(data=my.df, aes(x=xt, ymin=ymin, ymax=0))


###add Srx Prj support
tmp <- read.table("matrix/mers.p.s")
mers.sp <- tmp$V1
names(mers.sp)<-tmp$V2
mers.fdot<-read.table("matrix/mers.c.sg.proc.sort.cls.fdot",col.names = c("X","Cnt","Ratio","Srx"))
mers.fdot$Prj<-mers.sp[mers.fdot$Srx]
mers.fdot <- mers.fdot[mers.fdot$Ratio >0.005,]

p+geom_dotplot(data=mers.fdot, aes(x=X,fill=Prj,y=Cnt),stackdir="up",method = "histodot",
               binwidth = 100,stackgroups = T,binpositions = "all",dotsize = 1.5,
               position=position_nudge(y=mers.fdot$Ratio+0.04))+
  theme_classic()+theme(axis.title.x = element_blank(),
                        axis.line.x = element_blank(),
                        axis.ticks.x = element_blank(),
                        axis.text.x = element_blank())

