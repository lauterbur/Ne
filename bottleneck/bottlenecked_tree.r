### Show bottlenecked tree
library(phytools)
library(ape)
library(phylobase)
setwd("/home/lauterbur/Documents/Ne/bottleneck")


## Normal tree bottlenecked to 10%
tree<-pbtree(n=50)

plot(tree)
    tree4<-as(tree,"phylo4")
    survive<-sample(tree$tip.label,5)
    survive_nodes<-which(tree$tip.label %in% survive)

    survive_colors<-rep("grey",length(tree$edge.length))    
        survive_colors[c(which(tree$edge[,2] %in% survive_nodes),which(tree$edge[,2] %in% unlist(ancestors(tree4,survive))))]<-"black"    
    survive_width<-rep(5,length(tree$edge.length))
        survive_width[c(which(tree$edge[,2] %in% survive_nodes),which(tree$edge[,2] %in% unlist(ancestors(tree4,survive))))]<-10
    survive_points<-rep(0,length(tree$tip.label))
        survive_points[survive_nodes]<-1
            
png("bottlenecked-tree.png",height=750,width=1500)
    plot(tree,edge.color=survive_colors,show.tip.label=FALSE,edge.width=survive_width,type="cladogram",direction="downwards")
    tiplabels(frame="n",bg=NULL,tip=survive_nodes,pie=survive_points,piecol="red",offset=.05,cex=.5)    
dev.off()  

## Tree with polytomies bottlenecked to 10%
library(castor)
tree<-pbtree(n=50)
    tree<-merge_short_edges(tree,edge_length_epsilon=0.2)$tree
    is.rooted(tree)
    
    tree4<-as(tree,"phylo4")
    survive<-sample(tree$tip.label,5)
    survive_nodes<-which(tree$tip.label %in% survive)

    survive_colors<-rep("grey",length(tree$edge.length))    
        survive_colors[c(which(tree$edge[,2] %in% survive_nodes),which(tree$edge[,2] %in% unlist(ancestors(tree4,survive))))]<-"black"    
    survive_width<-rep(5,length(tree$edge.length))
        survive_width[c(which(tree$edge[,2] %in% survive_nodes),which(tree$edge[,2] %in% unlist(ancestors(tree4,survive))))]<-10
    survive_points<-rep(0,length(tree$tip.label))
        survive_points[survive_nodes]<-1
            
png("multifurcation-bottlenecked-tree.png",height=750,width=1500)
    plot(tree,edge.color=survive_colors,show.tip.label=FALSE,edge.width=survive_width,type="cladogram",direction="downwards")
    tiplabels(frame="n",bg=NULL,tip=survive_nodes,pie=survive_points,piecol="red",offset=.05,cex=.5)    
dev.off()  


## Tree in grid of points
x<-rep(seq(1,10),10)
y<-as.vector(t(sapply(rep(1,10), seq, 10)))
    x<-x[-c(c(1:3,8:10),c(1:3,8:10)+10,c(1:3,8:10)+20,c(1:3,8:10)+30,c(1:3,8:10)+40)]
    y<-y[-c(c(1:3,8:10),c(1:3,8:10)+10,c(1:3,8:10)+20,c(1:3,8:10)+30,c(1:3,8:10)+40)]

    
png("grid_severe_bottleneck.png",height=2000,width=2000)
par(mar=c(2,2,2,2))
plot(x,y,pch=19,cex=8, axes=FALSE, xlab=NULL)
    segments(x[1],y[1],x[6],y[6],col="red",lwd = 25)
    segments(x[2],y[2],x[6],y[6],col="red",lwd = 25)
    segments(x[3],y[3],x[6],y[6],col="red",lwd = 25)
    segments(x[4],y[4],x[8],y[8],col="red",lwd = 25)             
    segments(x[6],y[6],x[10],y[10],col="red",lwd = 25)
    segments(x[8],y[8],x[12],y[12],col="red",lwd = 25)
    segments(x[10],y[10],x[15],y[15],col="red",lwd = 25)    
    segments(x[12],y[12],x[15],y[15],col="red",lwd = 25)    
    segments(x[15],y[15],x[19],y[19],col="red",lwd = 25)    
    segments(x[19],y[19],x[25],y[25],col="red",lwd = 25)    
    segments(x[25],y[25],x[35],y[35],col="red",lwd = 25)    
    segments(x[35],y[35],x[45],y[45],col="red",lwd = 25)    
    segments(x[45],y[45],x[56],y[56],col="red",lwd = 25)    
    segments(x[56],y[56],x[66],y[66],col="red",lwd = 25)    
dev.off()


x<-rep(seq(1,10),10)
y<-as.vector(t(sapply(rep(1,10), seq, 10)))
    x<-x[-c(c(1,10),c(1,10)+10,c(1,10)+20,c(1,10)+30,c(1,10)+40)]
    y<-y[-c(c(1,10),c(1,10)+10,c(1,10)+20,c(1,10)+30,c(1,10)+40)]
    
png("grid_small_bottleneck.png",height=2000,width=2000)
par(mar=c(2,2,2,2))
plot(x,y,pch=19,cex=8, axes=FALSE, xlab=NULL)
    segments(x[1],y[1],x[10],y[10],col="red",lwd = 25)
    segments(x[3],y[3],x[11],y[11],col="red",lwd = 25)
    segments(x[6],y[6],x[14],y[14],col="red",lwd = 25)
    segments(x[7],y[7],x[16],y[16],col="red",lwd = 25)             
    segments(x[10],y[10],x[18],y[18],col="red",lwd = 25)
    segments(x[11],y[11],x[19],y[19],col="red",lwd = 25)
    segments(x[14],y[14],x[23],y[23],col="red",lwd = 25)    
    segments(x[16],y[16],x[24],y[24],col="red",lwd = 25)    
    segments(x[18],y[18],x[26],y[26],col="red",lwd = 25)    
    segments(x[19],y[19],x[26],y[26],col="red",lwd = 25)    
    segments(x[23],y[23],x[31],y[31],col="red",lwd = 25)    
    segments(x[24],y[24],x[32],y[32],col="red",lwd = 25)    
    segments(x[26],y[26],x[35],y[35],col="red",lwd = 25)    
    segments(x[31],y[31],x[39],y[39],col="red",lwd = 25)    
    segments(x[32],y[32],x[40],y[40],col="red",lwd = 25)    
    segments(x[35],y[35],x[45],y[45],col="red",lwd = 25)    
    segments(x[39],y[39],x[48],y[48],col="red",lwd = 25)    
    segments(x[40],y[40],x[48],y[48],col="red",lwd = 25)    
    segments(x[45],y[45],x[55],y[55],col="red",lwd = 25)    
    segments(x[48],y[48],x[57],y[57],col="red",lwd = 25)    
    segments(x[55],y[55],x[65],y[65],col="red",lwd = 25)    
    segments(x[57],y[57],x[67],y[67],col="red",lwd = 25)    
    segments(x[65],y[65],x[76],y[76],col="red",lwd = 25)    
    segments(x[67],y[67],x[76],y[76],col="red",lwd = 25)    
    segments(x[76],y[76],x[85],y[85],col="red",lwd = 25)    
dev.off()


## Palastra and Ruzzante figure redone
Conservation<-c(rep("<50",3),rep("<100",7),rep(">100",9),rep(">200",8),rep(">300",4),rep(">400",0),rep(">500",1),rep(">600",1),rep(">700",0),rep(">800",0),rep(">900",0))
Stable<-c(rep("<50",9),rep("<100",8),rep(">100",18),rep(">200",17),rep(">300",9),rep(">400",7),rep(">500",3),rep(">600",4),rep(">700",3),rep(">800",1),rep(">900",5),rep(">1000",4),rep(">1500",0),rep(">2000",1),rep(">2500",2),rep(">5000",3))
Exploited<-c(rep(">500",1),rep(">600",0),rep(">700",0),rep(">800",1),rep(">900",0),rep(">1000",3),rep(">1500",1),rep(">2000",0),rep(">2500",1),rep(">5000",1),rep(">10000",1))

data<-as.data.frame(cbind(Ne=Conservation,Category=rep("Conservation",length(Conservation))))
    data<-as.data.frame(rbind(data,cbind(Ne=Stable,Category=rep("Stable",length(Stable)))))
    data<-as.data.frame(rbind(data,cbind(Ne=Exploited,Category=rep("Exploited",length(Exploited)))))    
library(RColorBrewer)
library(dplyr)
mycolors<-brewer.pal(n = 3, name = "Dark2")
png("Palstra_Ruzzante.png",height=1000,width=2000)
ggplot(data %>% count(Ne,Category), aes(Ne,n,fill=Category)) + 
    geom_bar(stat="identity") +
    labs(y="Count",x=expression(N[e])) +
    theme(legend.title=element_blank(),legend.position = c(0.8, 0.8),
          text=element_text(size=90),axis.text=element_text(size=70),
          axis.text.x = element_text(angle = 45, margin=margin(t=100)),
          axis.title.x = element_text(margin = margin(t = -10, r = 0, b = 0, l = 0)),
          plot.margin = unit(c(10,220,10,10),"pt")) +
    scale_x_discrete(limits=c("<50","<100",">100",">200",">300",">400",">500",">600",">700",">800",">900",">1000",">1500",">2000",">2500",">5000",">10000")) +
    scale_fill_manual(values=mycolors)
dev.off()
                     