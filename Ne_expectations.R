### Calculating and graphing analytical expectations of probability of a single coalescent event in the previous generation

library(copula)
library(fields)
setwd("/home/lauterbur/Documents/Ne/")

## Probability that i lineages have i-1 ancestors
    ## Do for sample size (i) of 2, 5, 10, 20, 50
    ## Population size 10,20,50,100,200,500,1000,2000,5000,10000

fullGi<-vector("list",100)
approxGi<-vector("list",100)
error<-vector("list",100)
for (i in c(2,5,10,20,50,100)){
    for (N in seq(1,10000,1)){
        fullGi[[i]][N]<-(Stirling2(i,i-1)*fields.pochdown(N,i-1))/(N^i)
        approxGi[[i]][N]<-i*(i-1)/(2*N)
        error[[i]][N]<-abs(100-approxGi[[i]][N]*100/fullGi[[i]][N]) # approxGi percentage of fullGi
    }
}

## multiple merger probabilities
Gi<-list()
for (N in seq(1,1000,1)){
        count<-1
        Gi[[N]]<-list()
    for (i in c(5,10,20,50,100)){
        Gi[[N]][[count]]<-0
        for (j in seq(1,(i-2),1)){
            Gi_j<-(Stirling2(i,j)*fields.pochdown(N,j))/(N^i)
            Gi[[N]][[count]]<-sum(Gi[[N]][[count]]+Gi_j)
        }
        count<-count+1
    }
}

library(reshape2)
library(ggplot2)
Gi_df<-as.data.frame(do.call(rbind,lapply(Gi,unlist)))
    Gi_df<-cbind(seq(1,1000,1),Gi_df)
    colnames(Gi_df)<-c("N","5","10","20","50","100")
Gi_df_long<-melt(Gi_df,id.vars = "N")
    colnames(Gi_df_long)<-c("N","n","mulmerg_prob")

Gi_plot<-ggplot(Gi_df_long, aes(x=N,y=mulmerg_prob,color=n)) + 
    geom_line(aes(linetype = n), size=3) +
    scale_color_grey(start=.6,end=0) +
    theme_classic() +
    theme(text = element_text(size=25),
        axis.text.x = element_text(size=25),
        legend.key.width = unit(2.5, "cm"),
        legend.title.align=0.25) +
    xlab(expression(italic("N"["e"]))) +
    ylab("Probability of multiple mergers") +
    labs(name="sample size")
ggsave("multiple_merger_prob.eps",Gi_plot,width=13.3,height=10)
png("multiple_merger_prob.png",width=1000,height=750)
Gi_plot
dev.off()
    
    
for (i in c(2,5,10,20,50,100)){
    plot(seq(1,10000,1),Gi[[i]],type="l",xlab="N",ylab="multiple merger probability",xlim=c(0,1000),ylim=c(0,100))
        #abline(h=3.6,col="red")
        title(main=paste("i = ",i," genetic lineages",sep=""))
}

## sample sizes 5-50
png("n5-50-compareplots.png",height=1200,width=1200)
    par(mfrow=c(2,2),mar=c(5,8,4,2))
    for (i in c(5,10,20,50)){
    plot(seq(1,10000,1),error[[i]],type="l",xlab=expression("N"["e"]),ylab="% error",xlim=c(0,1000),ylim=c(0,1000),cex.lab=2,cex.axis=2,lwd=3)
        title(main=paste("n = ",i," genetic lineages",sep=""),cex.main=2)
    }
dev.off()


## back to % error
for (i in c(2,5,10,20,50,100)){
    plot(seq(1,10000,1),error[[i]],type="l",xlab="N",ylab="% error",xlim=c(0,1000),ylim=c(0,100))
        abline(h=3.6,col="red")
        title(main=paste("i = ",i," genetic lineages",sep=""))
}

for (i in c(2,5,10,20,50)){    
    plot(seq(1,10000,1),error[[i]],type="l",xlab="N",ylab="% error",ylim=c(0,4),xlim=c(1000,10000))
        abline(h=3.6,col="red")
        abline(v=1000,col="blue")
        title(main=paste("i = ",i," genetic lineages",sep=""))
}

## i = 10, N<1000 on left, N>1000 on right
png("i10-compareplots.png",height=600,width=1200)
    par(mfrow=c(1,2),mar=c(5,8,4,2))
    plot(seq(1,10000,1),error[[10]],type="l",xlab="N",ylab="% error",ylim=c(0,100),xlim=c(0,10000),cex.lab=2,cex.axis=2,lwd=3)
        abline(h=3.6,col="red",lwd=3)
        abline(v=1000,col="blue",lwd=3)
        title(main=paste("i = ","10"," genetic lineages",sep=""),cex.main=2)
        
    plot(seq(1,10000,1),error[[10]],type="l",xlab="N",ylab="% error",xlim=c(0,1000),ylim=c(0,100),cex.lab=2,cex.axis=2,lwd=3)
        abline(h=3.6,col="red",lwd=3)
        title(main=paste("i = ","10"," genetic lineages",sep=""),cex.main=2)
dev.off()

## four plots, N<1000, different i values
png("vary-i-compareplots.png",height=1200,width=1200)
    par(mfrow=c(2,2),mar=c(5,8,4,2))
    for (i in c(2,5,20,50)){
    plot(seq(1,10000,1),error[[i]],type="l",xlab="N",ylab="% error",xlim=c(0,1000),ylim=c(0,100),cex.lab=2,cex.axis=2,lwd=3)
        abline(h=3.6,col="red",lwd=3)
        title(main=paste("i = ",i," genetic lineages",sep=""),cex.main=2)
    }
dev.off()


png("vary-i-compareplots-scale.png",height=1200,width=1200)
    par(mfrow=c(2,2),mar=c(5,8,4,2))
    for (i in c(2,5,20,50)){
    plot(seq(1,10000,1),error[[i]],type="l",xlab="N",ylab="% error",xlim=c(0,1000),ylim=c(0,1000),cex.lab=2,cex.axis=2,lwd=3)
        abline(h=3.6,col="red",lwd=3)
        title(main=paste("i = ",i," genetic lineages",sep=""),cex.main=2)
    }
dev.off()

