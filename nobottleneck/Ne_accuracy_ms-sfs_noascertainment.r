## Compare the mean and variance of ms to SFS_code - look for if they are different, and how different
### RMSE and CV
### 10000 (polymorphic) iterations, 5000 bp loci (no ascertainment bias)

setwd("/home/lauterbur/Documents/Ne/nobottleneck/")
library(ggplot2)
#library(cowplot)
library(plyr)
library(dplyr)
#library(gridExtra)
library(reshape)
library(FSA)

# Root Mean Square Error
RMSE_function=function(obs,exp){
    sqrt(mean((exp-obs) ** 2))
}

# Multiplot function ------------------------------------------------------
        # Multiple plot function
            # http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)
            
            # ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
            # - cols:   Number of columns in layout
            # - layout: A matrix specifying the layout. If present, 'cols' is ignored.
            #
            # If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
            # then plot 1 will go in the upper left, 2 will go in the upper right, and
            # 3 will go all the way across the bottom.
            #
            multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
              library(grid)
            
              # Make a list from the ... arguments and plotlist
              plots <- c(list(...), plotlist)
            
              numPlots = length(plots)
            
              # If layout is NULL, then use 'cols' to determine layout
              if (is.null(layout)) {
                # Make the panel
                # ncol: Number of columns of plots
                # nrow: Number of rows needed, calculated from # of cols
                layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                                ncol = cols, nrow = ceiling(numPlots/cols))
              }
            
             if (numPlots==1) {
                print(plots[[1]])
            
              } else {
                # Set up the page
                grid.newpage()
                pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
            
                # Make each plot, in the correct location
                for (i in 1:numPlots) {
                  # Get the i,j matrix positions of the regions that contain this subplot
                  matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
            
                  print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                                  layout.pos.col = matchidx$col))
                }
              }
            }


# ## Import data ----------------------------------------------------------
ss_pi_ms_temp<-read.csv("collate-all-samplestats-ms.csv")
ss_pi_ms_temp<-na.omit(ss_pi_ms_temp)
    #colnames(ss_pi_ms_temp)<-c("Ne","n","locus","ss","pi")
ss_pi_ms<-data.frame()
    #ss_pi_ms<-ss_pi_ms[c(1:1000,10001:11000,20001:21000,30001:31000,40001:41000,50001:51000,60001:61000,70001:71000,80001:81000,90001:91000,100001:101000,110001:111000,120001:121000,103001:131000,140001:141000,150001:151000,160001:161000,170001:171000,180001:181000,190001:191000,200001:201000,210001:211000,220001:221000,230001:231000,240001:241000,250001:251000,260001:261000,270001:271000,280001:281000,290001:291000,300001:301000,310001:311000,320001:321000,330001:331000,340001:341000),]
    for (i in c(10,20,50,100,200,500,1000)){
        for (j in c(2,5,10,20,50,100,200,500,1000)){
                #numbersample<-min(10000,length(ss_pi_ms_temp[which(ss_pi_ms_temp$N==i & ss_pi_ms_temp$n==j & ss_pi_ms_temp$ss != 0),1])) # in case there are fewer than 10000 rows for that combination
                numbersample<-10000
                ss_pi_ms<-rbind(ss_pi_ms,ss_pi_ms_temp[sample(which(ss_pi_ms_temp$Ne==i & ss_pi_ms_temp$n==j),numbersample),])
        }
    }

ss_pi_sfs_temp<-read.csv("collate-all-samplestats-sfs.csv")
ss_pi_sfs_temp<-na.omit(ss_pi_sfs_temp)
#ss_pi_sfs_10000<-data.frame(Ne=rep(10000,60000),n=c(rep(2,10000),rep(5,10000),rep(10,10000),rep(20,10000),rep(50,10000),rep(100,10000)),ss=rep(0,60000),pi=rep(0,60000))
ss_pi_sfs<-data.frame()
    #ss_pi_sfs<-ss_pi_sfs[c(1:1000,10001:11000,20001:21000,30001:31000,40001:41000,50001:51000,60001:61000,70001:71000,80001:81000,90001:91000,100001:101000,110001:111000,120001:121000,103001:131000,140001:141000,150001:151000,160001:161000,170001:171000,180001:181000,190001:191000,200001:201000,210001:211000,220001:221000,230001:231000,240001:241000,250001:251000,260001:261000,270001:271000,280001:281000,290001:291000,300001:301000,310001:311000,320001:321000,330001:331000,340001:341000),]
    for (i in c(10,20,50,100,200,500,1000)){
        for (j in c(2,5,10,20,50,100,200,500,1000)){
            #numbersample<-min(10000,length(ss_pi_sfs_temp[which(ss_pi_sfs_temp$N==i & ss_pi_sfs_temp$n==j & ss_pi_sfs_temp$ss != 0),1])) # in case there are fewer than 10000 rows for that combination
            numbersample<-10000
            ss_pi_sfs<-rbind(ss_pi_sfs,ss_pi_sfs_temp[sample(which(ss_pi_sfs_temp$N==i & ss_pi_sfs_temp$n==j),numbersample),])
        }
    }

ss_pi_slim_temp<-read.csv("collate-all-samplestats-slim.csv")
ss_pi_slim_temp<-na.omit(ss_pi_slim_temp)
ss_pi_slim<-data.frame()
    for (i in c(10,20,50,100,200,500,1000)){
        for (j in c(2,5,10,20,50,100,200,500,1000)){
            numbersample<-10000
            ss_pi_slim<-rbind(ss_pi_slim,ss_pi_slim_temp[sample(which(ss_pi_slim_temp$N==i & ss_pi_slim_temp$n==j),numbersample,replace=TRUE),])
        }
    }

sample_sizes<-unique(ss_pi_ms$n)
pop_sizes<-unique(ss_pi_ms$Ne)

# Expected values ---------------------------------------------------------
## Expected values - these don't take into account the ascertainment bias
        ss_exp<-vector("list",length(pop_sizes))
        ss_WT<-vector("list",length(pop_sizes))
        pi_exp<-vector("list",length(pop_sizes))
        watt_exp<-vector("list",length(pop_sizes))
        ss_exp_df<-vector("list",length(pop_sizes))
        ss_WT_df<-vector("list",length(pop_sizes))
        pi_exp_df<-vector("list",length(pop_sizes))
        watt_exp_df<-vector("list",length(pop_sizes))
        mu<-0.000001 # highu
        pop_size_iter<-0
        for (pop_size in pop_sizes){
            pop_size_iter<-pop_size_iter+1
            sample_size_iter<-0
            theta<-4*pop_size*mu
            for (sample_size in sample_sizes){
                sample_size_iter<-sample_size_iter+1
                x<-(sample_size*2)/(pop_size*2) # for diploid
                sub<-1-exp(-x)
                g<-(1/6)*sub+(1/36)*sub^2+(37/(36*180))*sub^3+(205/(24*10800))*sub^4-(21625/(1134000*120))*sub^5-(1.54402*10^-4)*sub^6-(2.463*10^-5)*sub^7+(3.1981*10^-5)*sub^8+(1.9962*10^-5)*sub^9
                tau<-log(sample_size*2)+0.57721566+.5*g # for diploid
                h<-sample_size*2 # for diploid
                    multiplier<-0
                    for (k in 1:(h-1)){
                        multiplier<-multiplier+1/k
                    }
                ss_exp[[pop_size_iter]][sample_size_iter]<-theta*multiplier*5000 # multiply by number of sites and loci
                ss_WT[[pop_size_iter]][sample_size_iter]<-theta*5000*tau
                pi_exp[[pop_size_iter]][sample_size_iter]<-theta*5000 # multiply by number of sites
                watt_exp[[pop_size_iter]][sample_size_iter]<-theta*5000
            }
            ss_exp_df[[pop_size_iter]]<-cbind.data.frame(sample_sizes,ss_exp[[pop_size_iter]])
                colnames(ss_exp_df[[pop_size_iter]])<-c("n","ss_exp")
            ss_WT_df[[pop_size_iter]]<-cbind.data.frame(sample_sizes,ss_WT[[pop_size_iter]])
                colnames(ss_WT_df[[pop_size_iter]])<-c("n","ss_WT")
            pi_exp_df[[pop_size_iter]]<-cbind.data.frame(sample_sizes,pi_exp[[pop_size_iter]])
                colnames(pi_exp_df[[pop_size_iter]])<-c("n","pi_exp")
            watt_exp_df[[pop_size_iter]]<-cbind.data.frame(sample_sizes,watt_exp[[pop_size_iter]])
                colnames(watt_exp_df[[pop_size_iter]])<-c("n","watt_exp")
        }

# Segregating sites -------------------------------------------------------
       ## make sure to calculate total ss, avg pi over each locus, avg theta over each locus!
    ## T-tests and distribution plots comparing sfs and ms - segregating sites
        ss_data<-rbind(ss_pi_ms[,c(1:3)],ss_pi_sfs[,c(1:3)],ss_pi_slim[,c(1:3)]) # column 1=N, 2=n, 3=ss, 4=pi
        ss_data<-cbind(ss_data,c(rep("ms",length(ss_pi_ms[,1])),rep("sfs",length(ss_pi_sfs[,1])),rep("slim",length(ss_pi_sfs[,1]))))
        colnames(ss_data)<-c("N","n","seg_sites","simulation")
        ss_data$n_over_Ne<-ss_data$n/ss_data$N
        ss_data$N<-as.factor(ss_data$N)

    # par(mfrow=c(2,3))
        ss_means_sfs<-list()
        ss_means_ms<-list()
        ss_anovas<-list()
        ss_cv_sfs<-list()
        ss_cv_ms<-list()
        ss_tukey_ms_sfs<-list()
        ks_result<-list()
        ks_result_p<-list()
        RMSE<-list()
        
        pop_size_iter<-0
        for (pop_size in pop_sizes){
            pop_size_iter<-pop_size_iter+1
            sample_size_iter<-0
                ss_means_sfs[[pop_size_iter]]<-list()
                ss_means_ms[[pop_size_iter]]<-list()
                ss_anovas[[pop_size_iter]]<-list()
                ss_cv_sfs[[pop_size_iter]]<-list()
                ss_cv_ms[[pop_size_iter]]<-list()
                ss_tukey_ms_sfs[[pop_size_iter]]<-list()
                ks_result[[pop_size_iter]]<-list()
                ks_result_p[[pop_size_iter]]<-list()
                RMSE[[pop_size_iter]]<-list()
                for (sample_size in sample_sizes){
                    sample_size_iter<-sample_size_iter+1
                            ss_means_sfs[[pop_size_iter]][[sample_size_iter]]<-mean(ss_data[which(ss_data$n==sample_size & ss_data$N==pop_size & ss_data$simulation=="sfs"),3])
                            ss_means_sfs[[pop_size_iter]][[sample_size_iter]][2]<-sd(ss_data[which(ss_data$n==sample_size & ss_data$N==pop_size & ss_data$simulation=="sfs"),3])
                            ss_means_ms[[pop_size_iter]][[sample_size_iter]]<-mean(ss_data[which(ss_data$n==sample_size & ss_data$N==pop_size & ss_data$simulation=="ms"),3])
                            ss_means_ms[[pop_size_iter]][[sample_size_iter]][2]<-sd(ss_data[which(ss_data$n==sample_size & ss_data$N==pop_size & ss_data$simulation=="ms"),3])                        
                            ss_cv_sfs[[pop_size_iter]][[sample_size_iter]]<-(ss_means_sfs[[pop_size_iter]][[sample_size_iter]][2]/ss_means_sfs[[pop_size_iter]][[sample_size_iter]][1])*100
                            ss_cv_ms[[pop_size_iter]][[sample_size_iter]]<-(ss_means_ms[[pop_size_iter]][[sample_size_iter]][2]/ss_means_ms[[pop_size_iter]][[sample_size_iter]][1])*100
                            ss_anovas[[pop_size_iter]][[sample_size_iter]]<-summary(aov(ss_data[which(ss_data$n==sample_size & ss_data$N==pop_size),3]~ss_data[which(ss_data$n==sample_size & ss_data$N==pop_size),4]))
                            ss_tukey_ms_sfs[[pop_size_iter]][[sample_size_iter]]<-TukeyHSD(aov(ss_data[which(ss_data$n==sample_size & ss_data$N==pop_size),3]~ss_data[which(ss_data$n==sample_size & ss_data$N==pop_size),4]))
                            ks_result[[pop_size_iter]][[sample_size_iter]]<-ks.test(ss_data[which(ss_data$N==pop_size & ss_data$n==sample_size & ss_data$simulation=="ms"),3],ss_data[which(ss_data$N==pop_size & ss_data$n==sample_size & ss_data$simulation=="sfs"),3])
                            ks_result_p[[pop_size_iter]][[sample_size_iter]]<-ks.test(ss_data[which(ss_data$N==pop_size & ss_data$n==sample_size & ss_data$simulation=="ms"),3],ss_data[which(ss_data$N==pop_size & ss_data$n==sample_size & ss_data$simulation=="sfs"),3])$p.value
                            RMSE[[pop_size_iter]][[sample_size_iter]]<-RMSE_function(ss_data[which(ss_data$N==pop_size & ss_data$n==sample_size & ss_data$simulation=="sfs"),3],ss_data[which(ss_data$N==pop_size & ss_data$n==sample_size & ss_data$simulation=="ms"),3])
                    }
                }

            ## adjusted p-values
                pvalues<-unlist(unlist(unlist(ss_anovas,recursive=F),recursive=F),recursive=F)
                ## fix below for added time
                    p.adj<-p.adjust(unlist(pvalues[which(names(pvalues)=="Pr(>F)")])[which(names(unlist(pvalues[which(names(pvalues)=="Pr(>F)")]))=="Pr(>F)1")],method="fdr")
                    p.adj.list<-relist(unname(p.adj),skeleton=rep(list(rep("a",9)),7))
                pvaluescoalescent<-unlist(lapply(unname(unlist(unlist(unlist(ss_tukey_ms_sfs,recursive=F),recursive=F),recursive=F)),unname))[seq(10,9*7*12,12)]
                    p.adj.coalescent<-p.adjust(pvaluescoalescent,method="fdr")
                    p.adj.coalescent.list<-relist(unname(p.adj.coalescent),list(rep(list(rep("a",9)),7)))
                pvaluesforward<-unlist(lapply(unname(unlist(unlist(unlist(ss_tukey_ms_sfs,recursive=F),recursive=F),recursive=F)),unname))[seq(12,9*7*12,12)]
                    p.adj.forward<-p.adjust(pvaluesforward,method="fdr")
                    p.adj.forward.list<-relist(unname(p.adj.forward),list(rep(list(rep("a",9)),7)))
                kspvalues<-unlist(ks_result_p)
                ks.p.adj<-p.adjust(kspvalues,method="fdr")
                ks.p.adj.list<-relist(unname(ks.p.adj),skeleton=list(rep(list(rep("a",9)),7)))
        ## Density plots of each distribution

        palette<-c("orange","blue","green","turquoise")
        dp<-list()
        vp<-list()
        n_over_Ne_diff<-data.frame(n_over_Ne=unique(ss_data$n_over_Ne)[sort.list(unique(ss_data$n_over_Ne))],check.names=FALSE)
        namevector<-rep("ss.diff",7)
        n_over_Ne_cv_ms<-data.frame(n_over_Ne=unique(ss_data$n_over_Ne)[sort.list(unique(ss_data$n_over_Ne))],check.names=FALSE)
        n_over_Ne_cv_sfs<-data.frame(n_over_Ne=unique(ss_data$n_over_Ne)[sort.list(unique(ss_data$n_over_Ne))],check.names=FALSE)
        namevector_cv<-rep("ss.cv",7)
        n_over_Ne_diff[,namevector]<-NA
        n_over_Ne_cv_ms[,namevector_cv]<-NA
        n_over_Ne_cv_sfs[,namevector_cv]<-NA
       # WT_diff_sfs<-c()
        #exp_diff_sfs<-c()
        pop_size_iter<-0
        for (pop_size in pop_sizes){
            pop_size_iter<-pop_size_iter+1
            sample_size_iter<-0
            ## density plots
                dp[[pop_size_iter]]<-list()
                vp[[pop_size_iter]]<-list()
                for (sample_size in sample_sizes){
                    sample_size_iter<-sample_size_iter+1
                        means<-ddply(ss_data[which(ss_data$n==sample_size & ss_data$N==pop_size),],"simulation",summarise,seg_sites.mean=mean(seg_sites))
                       # cvs<-data.frame(simulation=c("ms","sfs","slim","exp"),seg_sites.cv="NA")
                        # levels(means$simulation)<-c(levels(means$simulation),"expected")
                            # means<-rbind(means,c("expected",ss_exp_df[[pop_size_iter]][l,2]))
                            # means$seg_sites.mean<-as.numeric(means$seg_sites.mean)
                            # levels(ss_data$simulation)<-c(levels(ss_data$simulation),"expected")
                            # ss_data<-rbind(ss_data,c("10",2,0,"expected"))
                            #     ss_data$seg_sites<-as.numeric(ss_data$seg_sites)
                       # exp<-data.frame("exp",ss_exp_df[[pop_size_iter]][sample_size_iter,2])
                            #names(exp)<-c("simulation","seg_sites.mean")
                       # WT<-data.frame("exp_WT",ss_WT_df[[pop_size_iter]][sample_size_iter,2])
                            #names(WT)<-c("simulation","seg_sites.mean")
                        #means<-rbind(means,exp,WT)
                      #  WT_diff_sfs[sample_size_iter]<-((means[2,2]-means[5,2])/means[2,2])*100
                      #  exp_diff_sfs[sample_size_iter]<-((means[2,2]-means[4,2])/means[2,2])*100
                        nNe<-sample_size/as.numeric(as.character(pop_size))
                        perc_diff<-((means[which(means$simulation=="ms"),2]-means[which(means$simulation=="sfs"),2])/means[which(means$simulation=="sfs"),2])*100
                        n_over_Ne_diff[which(n_over_Ne_diff$n_over_Ne==nNe),pop_size_iter+1]<-perc_diff
                        n_over_Ne_cv_ms[which(n_over_Ne_cv_ms$n_over_Ne==nNe),pop_size_iter+1]<-ss_cv_ms[[pop_size_iter]][[sample_size_iter]]
                        n_over_Ne_cv_sfs[which(n_over_Ne_cv_sfs$n_over_Ne==nNe),pop_size_iter+1]<-ss_cv_sfs[[pop_size_iter]][[sample_size_iter]]
                        dp[[pop_size_iter]][[sample_size_iter]]<-ggplot(ss_data[which(ss_data$n==sample_size & ss_data$N==pop_size),], aes(x=seg_sites, fill=simulation)) + geom_density(alpha=.3) +
                            #geom_vline(data=means,aes(xintercept=seg_sites.mean,color=simulation),linetype="dashed",size=1) +
                            #geom_vline(aes(xintercept=as.numeric(pop_size)*4*0.000001*10), color="#BB0000") +
                            #geom_vline(aes(xintercept=ss_means_sfs[[pop_size_iter]][[sample_size_iter]][1],color="mean sfs"),linetype="dashed",size=1) +
                           # geom_point(data=means,aes(x=seg_sites.mean,y=0,color=simulation),alpha=1,size=2) +
                            labs(x="Segregating Sites",title=paste("Ne = ",pop_size,", n = ",sample_size,", T = ",time,sep=""),subtitle=paste("Difference between means, FDR = ",signif(p.adj.list[[pop_size_iter]][[sample_size_iter]],3),"\n",signif(perc_diff,3),"% difference between sfs and ms means\n",sep="")) +
                            scale_fill_manual(values=palette[c(1,2)]) +
                            scale_colour_manual(values=palette[c(1,2)])
                       # vp[[pop_size_iter]][[sample_size_iter]]<-ggplot(ss_data[which(ss_data$n==sample_size & ss_data$N==pop_size),], aes(x=simulation, y=seg_sites)) +
                        #    geom_violin(aes(fill = simulation)) + 
                         #   geom_hline(data=means[c(1:2,4:5),],aes(yintercept=seg_sites.mean,color=simulation),linetype="dashed") +
                          #  labs(x="Simulation",y="Segregating Sites",title=paste("Ne = ",pop_size,", n = ",sample_size,sep=""),subtitle=paste("Difference between means, FDR = ",signif(p.adj.list[[pop_size_iter]][[sample_size_iter]],3),"\n",signif(perc_diff,3),"% difference between means",sep="")) +
                           # scale_fill_manual(values=palette) +
                            #scale_colour_manual(values=palette)
                        #dp[[pop_size_iter]][[sample_size_iter]]
                    #    vp[[pop_size_iter]][[sample_size_iter]]
                }
               # plot(sample_sizes/pop_size,WT_diff_sfs,main=paste("Ne=",pop_size,sep=""))
               # points(sample_sizes/pop_size,exp_diff_sfs,col="red")
            }


            # png("distribution-comparisions-violin-ss-msvssfs.png",height=3000,width=3000)
            #     multiplot(plotlist=unlist(vp,recursive=FALSE),cols=7)
            # dev.off()

            ## n/Ne
            # new_nNe_diff<-melt(n_over_Ne_diff,id="n_over_Ne",na.rm=TRUE)
            #     new_nNe_diff<-new_nNe_diff[,-2]
            #     colnames(new_nNe_diff)<-c("n_over_Ne","perc_diff")
            # png("n_over_Ne_ss-msvssfs.png")
            #     ggplot(data=new_nNe_diff,aes(x=as.factor(n_over_Ne),y=perc_diff,fill=perc_diff)) +
            #         geom_point(shape=21,size=3) +
            #         scale_fill_gradient(low="orange", high="blue") +
            #         scale_color_gradient(low="orange",high="blue") +
            #         labs(x="n/Ne",y="Percent difference in number of segregating sites")
            # dev.off()


# Watterson's theta -------------------------------------------------------
    ## T-tests and distribution plots comparing sfs and ms - Watterson's theta
    # Watterson's theta = S/a_n, a_n = sum(1/i), i = {1..n-1}
        a_n<-rep(0,length(unique(ss_data[,2])))
            for (i in 1:length(unique(ss_data[,2]))){
                k<-unique(ss_data[,2])[i]*2
                for (j in 1:(k-1)){
                    a_n[i]<-a_n[i]+(1/j)
                }
            }

        # theta_data_10<-ss_data_10
        # colnames(theta_data_10)[3]<-"wtheta"
        # j<-1
        # for (i in unique(ss_data_10[,2])){
        #     theta_data_10$wtheta[which(theta_data_10$n==i)]<-((ss_data_10$seg_sites[which(theta_data_10$n==i)])/a_n[j])/10 # divide by number of loci
        # j<-j+1
        # }
        theta_data<-ss_data
        colnames(theta_data)[3]<-"theta"
        j<-1
        for (i in unique(ss_data[,2])){
            theta_data$theta[which(theta_data$n==i)]<-((ss_data$seg_sites[which(theta_data$n==i)])/a_n[j]) # divide by number of loci
        j<-j+1
        }
            
    # par(mfrow=c(2,3))   
    # par(mfrow=c(2,3))
        theta_means_sfs<-list()
        theta_means_ms<-list()
        theta_means_slim<-list()
        theta_anovas<-list()
        theta_cv_sfs<-list()
        theta_cv_ms<-list()
        theta_tukey_ms_sfs<-list()
        theta_tukey_slim_sfs<-list()
        ks_result<-list()
        ks_result_slim<-list()
        ks_result_p<-list()
        ks_result_p_slim<-list()
        kruskalwallis<-list()
        kruskalwallis_p<-list()
        dunn<-list()
        dunn_p<-list()
        dunn_p_adj<-list()
        RMSE<-list()
        
        pop_size_iter<-0
        for (pop_size in pop_sizes){
            pop_size_iter<-pop_size_iter+1
            sample_size_iter<-0
                theta_means_sfs[[pop_size_iter]]<-list()
                theta_means_ms[[pop_size_iter]]<-list()
                theta_means_slim[[pop_size_iter]]<-list()
                theta_anovas[[pop_size_iter]]<-list()
                theta_cv_sfs[[pop_size_iter]]<-list()
                theta_cv_ms[[pop_size_iter]]<-list()
                theta_tukey_ms_sfs[[pop_size_iter]]<-list()
                theta_tukey_slim_sfs[[pop_size_iter]]<-list()
                ks_result[[pop_size_iter]]<-list()
                ks_result_p[[pop_size_iter]]<-list()
                ks_result_slim[[pop_size_iter]]<-list()
                ks_result_p_slim[[pop_size_iter]]<-list()
                kruskalwallis[[pop_size_iter]]<-list()
                kruskalwallis_p[[pop_size_iter]]<-list()
                dunn[[pop_size_iter]]<-list()
                dunn_p[[pop_size_iter]]<-list()
                dunn_p_adj[[pop_size_iter]]<-list()
                RMSE[[pop_size_iter]]<-list()
                for (sample_size in sample_sizes){
                    sample_size_iter<-sample_size_iter+1
                            theta_means_sfs[[pop_size_iter]][[sample_size_iter]]<-mean(theta_data[which(theta_data$n==sample_size & theta_data$N==pop_size & theta_data$simulation=="sfs"),3])
                            theta_means_sfs[[pop_size_iter]][[sample_size_iter]][2]<-sd(theta_data[which(theta_data$n==sample_size & theta_data$N==pop_size & theta_data$simulation=="sfs"),3])
                            theta_means_slim[[pop_size_iter]][[sample_size_iter]]<-mean(theta_data[which(theta_data$n==sample_size & theta_data$N==pop_size & theta_data$simulation=="slim"),3])
                            theta_means_slim[[pop_size_iter]][[sample_size_iter]][2]<-sd(theta_data[which(theta_data$n==sample_size & theta_data$N==pop_size & theta_data$simulation=="slim"),3])                        
                            theta_means_ms[[pop_size_iter]][[sample_size_iter]]<-mean(theta_data[which(theta_data$n==sample_size & theta_data$N==pop_size & theta_data$simulation=="ms"),3])
                            theta_means_ms[[pop_size_iter]][[sample_size_iter]][2]<-sd(theta_data[which(theta_data$n==sample_size & theta_data$N==pop_size & theta_data$simulation=="ms"),3])                        
                            theta_cv_sfs[[pop_size_iter]][[sample_size_iter]]<-(theta_means_sfs[[pop_size_iter]][[sample_size_iter]][2]/theta_means_sfs[[pop_size_iter]][[sample_size_iter]][1])*100
                            theta_cv_ms[[pop_size_iter]][[sample_size_iter]]<-(theta_means_ms[[pop_size_iter]][[sample_size_iter]][2]/theta_means_ms[[pop_size_iter]][[sample_size_iter]][1])*100
                            theta_anovas[[pop_size_iter]][[sample_size_iter]]<-summary(aov(theta_data[which(theta_data$n==sample_size & theta_data$N==pop_size),3]~theta_data[which(theta_data$n==sample_size & theta_data$N==pop_size),4]))
                            theta_tukey_ms_sfs[[pop_size_iter]][[sample_size_iter]]<-TukeyHSD(aov(theta_data[which(theta_data$n==sample_size & theta_data$N==pop_size),3]~theta_data[which(theta_data$n==sample_size & theta_data$N==pop_size),4]))
                           # theta_tukey_slim_sfs[[pop_size_iter]][[sample_size_iter]]<-TukeyHSD(aov(theta_data[which(theta_data$n==sample_size & theta_data$N==pop_size & !theta_data$simulation=="coalescent"),3]~theta_data[which(theta_data$n==sample_size & theta_data$N==pop_size & !theta_data$simulation=="coalescent"),4]))
                            ks_result[[pop_size_iter]][[sample_size_iter]]<-ks.test(theta_data[which(theta_data$N==pop_size & theta_data$n==sample_size & theta_data$simulation=="ms"),3],theta_data[which(theta_data$N==pop_size & theta_data$n==sample_size & theta_data$simulation=="sfs"),3])
                            ks_result_p[[pop_size_iter]][[sample_size_iter]]<-ks.test(theta_data[which(theta_data$N==pop_size & theta_data$n==sample_size & theta_data$simulation=="ms"),3],theta_data[which(theta_data$N==pop_size & theta_data$n==sample_size & theta_data$simulation=="sfs"),3])$p.value
                            ks_result_slim[[pop_size_iter]][[sample_size_iter]]<-ks.test(theta_data[which(theta_data$N==pop_size & theta_data$n==sample_size & theta_data$simulation=="slim"),3],theta_data[which(theta_data$N==pop_size & theta_data$n==sample_size & theta_data$simulation=="sfs"),3])
                            ks_result_p_slim[[pop_size_iter]][[sample_size_iter]]<-ks.test(theta_data[which(theta_data$N==pop_size & theta_data$n==sample_size & theta_data$simulation=="slim"),3],theta_data[which(theta_data$N==pop_size & theta_data$n==sample_size & theta_data$simulation=="sfs"),3])$p.value
                            kruskalwallis[[pop_size_iter]][[sample_size_iter]]<-kruskal.test(theta_data[which(theta_data$n==sample_size & theta_data$N==pop_size),3]~theta_data[which(theta_data$n==sample_size & theta_data$N==pop_size),4])
                            kruskalwallis_p[[pop_size_iter]][[sample_size_iter]]<-kruskal.test(theta_data[which(theta_data$n==sample_size & theta_data$N==pop_size),3]~theta_data[which(theta_data$n==sample_size & theta_data$N==pop_size),4])$p.value
                            dunn[[pop_size_iter]][[sample_size_iter]]<-dunnTest(theta_data[which(theta_data$n==sample_size & theta_data$N==pop_size),3]~theta_data[which(theta_data$n==sample_size & theta_data$N==pop_size),4],method="bh")
                            dunn_p[[pop_size_iter]][[sample_size_iter]]<-dunnTest(theta_data[which(theta_data$n==sample_size & theta_data$N==pop_size),3]~theta_data[which(theta_data$n==sample_size & theta_data$N==pop_size),4])$res$P.unadj
                            dunn_p_adj[[pop_size_iter]][[sample_size_iter]]<-dunnTest(theta_data[which(theta_data$n==sample_size & theta_data$N==pop_size),3]~theta_data[which(theta_data$n==sample_size & theta_data$N==pop_size),4])$res$P.adj
                            RMSE[[pop_size_iter]][[sample_size_iter]]<-RMSE_function(theta_data[which(theta_data$N==pop_size & theta_data$n==sample_size & theta_data$simulation=="sfs"),3],theta_data[which(theta_data$N==pop_size & theta_data$n==sample_size & theta_data$simulation=="ms"),3])
                    
                }
            }

            ## adjusted p-values
                pvalues<-unlist(unlist(unlist(unlist(theta_anovas,recursive=F),recursive=F),recursive=F),recursive=F)
                pvalues_kw<-unlist(unlist(unlist(unlist(kruskalwallis_p,recursive=F),recursive=F),recursive=F),recursive=F)
                ## fix below for added time
                    p.adj<-p.adjust(unlist(pvalues[which(names(pvalues)=="Pr(>F)1")])[which(names(unlist(pvalues[which(names(pvalues)=="Pr(>F)1")]))=="Pr(>F)1")],method="fdr")
                    p.adj.list<-relist(round(unname(p.adj),5),skeleton=rep(list(rep("a",9)),7))
                    p.adj.kw<-p.adjust(pvalues_kw, method="fdr")
                    p.adj.kw.list<-relist(round(unname(p.adj.kw),5),skeleton=rep(list(rep("a",9)),7))
                pvaluescoalescent<-unlist(lapply(unname(unlist(unlist(unlist(unlist(theta_tukey_ms_sfs,recursive=F),recursive=F),recursive=F),recursive=F)),unname))[seq(10,9*7*12,12)]
                    p.adj.coalescent<-p.adjust(pvaluescoalescent,method="fdr")
                    p.adj.coalescent.list<-relist(round(unname(p.adj.coalescent),5),rep(list(rep("a",9)),7))
                pvaluesforward<-unlist(lapply(unname(unlist(unlist(unlist(unlist(theta_tukey_ms_sfs,recursive=F),recursive=F),recursive=F),recursive=F)),unname))[seq(12,9*7*12,12)]
                    p.adj.forward<-p.adjust(pvaluesforward,method="fdr")
                    p.adj.forward.list<-relist(round(unname(p.adj.forward),5),rep(list(rep("a",9)),7))
                pvaluesslim<-unlist(lapply(unname(unlist(unlist(unlist(unlist(theta_tukey_ms_sfs,recursive=F),recursive=F),recursive=F),recursive=F)),unname))[seq(11,9*7*12,12)]
                    p.adj.slim<-p.adjust(pvaluesslim,method="fdr")
                    p.adj.slim.list<-relist(round(unname(p.adj.slim),5),rep(list(rep("a",9)),7))
                kspvalues<-unlist(ks_result_p)
                ks.p.adj<-p.adjust(kspvalues,method="fdr")
                ks.p.adj.list<-relist(round(unname(ks.p.adj),5),skeleton=rep(list(rep("a",9)),7))
                kspvaluesslim<-unlist(ks_result_p_slim)
                ks.p.adj.slim<-p.adjust(kspvaluesslim,method="fdr")
                ks.p.adj.list.slim<-relist(round(unname(ks.p.adj.slim),5),skeleton=rep(list(rep("a",9)),7))
                
        ## Density plots of each distribution

        palette<-c("orange","blue","green","turquoise")
        dp<-list()
        vp<-list()
        n_over_Ne_diff<-list()
        #n_over_Ne_cv_ms<-list()
        #n_over_Ne_cv_sfs<-list()
       # WT_diff_sfs<-c()
        #exp_diff_sfs<-c()
        n_over_Ne_diff<-data.frame(n_over_Ne=unique(theta_data$n_over_Ne)[sort.list(unique(theta_data$n_over_Ne))],check.names=FALSE)
        namevector<-rep("theta.diff",7)
        #n_over_Ne_cv_ms<-data.frame(n_over_Ne=unique(theta_data$n_over_Ne)[sort.list(unique(theta_data$n_over_Ne))],check.names=FALSE)
        #n_over_Ne_cv_sfs<-data.frame(n_over_Ne=unique(theta_data$n_over_Ne)[sort.list(unique(theta_data$n_over_Ne))],check.names=FALSE)
        #namevector_cv<-rep("theta.cv",7)
        n_over_Ne_diff[,namevector]<-NA
        #n_over_Ne_cv_ms[,namevector_cv]<-NA
        #n_over_Ne_cv_sfs[,namevector_cv]<-NA
        pop_size_iter<-0
        dp<-list()
        vp<-list()
        for (pop_size in pop_sizes){
            pop_size_iter<-pop_size_iter+1
            sample_size_iter<-0
            ## density plots
                dp[[pop_size_iter]]<-list()
                vp[[pop_size_iter]]<-list()
                for (sample_size in sample_sizes){
                    sample_size_iter<-sample_size_iter+1
                    dp[[pop_size_iter]][[sample_size_iter]]<-list()
                    vp[[pop_size_iter]][[sample_size_iter]]<-list()
                    means<-ddply(theta_data[which(theta_data$n==sample_size & theta_data$N==pop_size),],"simulation",summarise,theta.mean=mean(theta))
                       # cvs<-data.frame(simulation=c("ms","sfs","slim","exp"),theta.cv="NA")
                        # levels(means$simulation)<-c(levels(means$simulation),"expected")
                            # means<-rbind(means,c("expected",theta_exp_df[[pop_size_iter]][l,2]))
                            # means$theta.mean<-as.numeric(means$theta.mean)
                            # levels(theta_data$simulation)<-c(levels(theta_data$simulation),"expected")
                            # theta_data<-rbind(theta_data,c("10",2,0,"expected"))
                            #     theta_data$theta<-as.numeric(theta_data$theta)
                       # exp<-data.frame("exp",theta_exp_df[[pop_size_iter]][sample_size_iter,2])
                           # names(exp)<-c("simulation","theta.mean")
                       # WT<-data.frame("exp_WT",theta_WT_df[[pop_size_iter]][sample_size_iter,2])
                           # names(WT)<-c("simulation","theta.mean")
                        #means<-rbind(means,exp,WT)
                      #  WT_diff_sfs[sample_size_iter]<-((means[2,2]-means[5,2])/means[2,2])*100
                      #  exp_diff_sfs[sample_size_iter]<-((means[2,2]-means[4,2])/means[2,2])*100
                        nNe<-sample_size/as.numeric(as.character(pop_size))
                        perc_diff<-((means[which(means$simulation=="ms"),2]-means[which(means$simulation=="sfs"),2])/means[which(means$simulation=="sfs"),2])*100
                        n_over_Ne_diff[which(n_over_Ne_diff$n_over_Ne==nNe),pop_size_iter+1]<-perc_diff
                        #n_over_Ne_cv_ms[which(n_over_Ne_cv_ms$n_over_Ne==nNe),pop_size_iter+1]<-theta_cv_ms[[pop_size_iter]][[sample_size_iter]]
                        #n_over_Ne_cv_sfs[which(n_over_Ne_cv_sfs$n_over_Ne==nNe),pop_size_iter+1]<-theta_cv_sfs[[pop_size_iter]][[sample_size_iter]]
                        dp[[pop_size_iter]][[sample_size_iter]]<-ggplot(theta_data[which(theta_data$n==sample_size & theta_data$N==pop_size),], aes(x=theta, fill=simulation)) + geom_density(alpha=.3) +
                            geom_vline(data=means,aes(xintercept=theta.mean,color=simulation),linetype="dashed",size=1) +
                            #geom_vline(aes(xintercept=as.numeric(pop_size)*4*0.000001*10), color="#BB0000") +
                            #geom_vline(aes(xintercept=theta_means_sfs[[pop_size_iter]][[sample_size_iter]][1],color="mean sfs"),linetype="dashed",size=1) +
                           # geom_point(data=means,aes(x=theta.mean,y=0,color=simulation),alpha=1,size=2) +
                            labs(x="Watterson's Theta",title=paste("Ne = ",pop_size,", n = ",sample_size,sep=""),subtitle=paste("Difference between means, FDR = ",signif(p.adj.list[[pop_size_iter]][[sample_size_iter]],3),"\n",signif(perc_diff,3),"% difference between sfs and ms means\n",sep="")) +
                            scale_fill_manual(values=palette[c(1,2,3)]) +
                            scale_colour_manual(values=palette[c(1,2,3)])
                       # vp[[pop_size_iter]][[sample_size_iter]]<-ggplot(theta_data[which(theta_data$n==sample_size & theta_data$N==pop_size),], aes(x=simulation, y=theta)) +
                        #    geom_violin(aes(fill = simulation)) + 
                         #   geom_hline(data=means[c(1:2,4:5),],aes(yintercept=theta.mean,color=simulation),linetype="dashed") +
                          #  labs(x="Simulation",y="Segregating Sites",title=paste("Ne = ",pop_size,", n = ",sample_size,sep=""),subtitle=paste("Difference between means, FDR = ",signif(p.adj.list[[pop_size_iter]][[sample_size_iter]],3),"\n",signif(perc_diff,3),"% difference between means",sep="")) +
                           # scale_fill_manual(values=palette) +
                            #scale_colour_manual(values=palette)
                        #dp[[pop_size_iter]][[sample_size_iter]]
                    #    vp[[pop_size_iter]][[sample_size_iter]]
                    }
                }

            
    ## Boxplot with Ne on x, theta on y
        bp<-list()
        sample_size_iter<-0
        for (s in sample_sizes){
            sample_size_iter<-sample_size_iter+1
            bp[[sample_size_iter]]<-ggplot(theta_data[which(theta_data$n==s),], aes(x=N,y=theta,fill=simulation)) + geom_boxplot() +
                    labs(y="Watterson's theta",title=paste("n =",s))
        }
            png("n-distribution-comparisons-theta-msvssfs-nobottleneck.png",height=1000,width=2000)
                multiplot(plotlist=bp,cols=3)
            dev.off()
            png("n10-distribution-comparisons-theta-msvssfs-nobottleneck.png")
                bp[[3]]
            dev.off()
            
            theta_data<-cbind(theta_data,T=rep("Constant Ne",length(theta_data$N)))
            bp_n<-ggplot(theta_data[which(theta_data$n==20 & theta_data$simulation!="slim"),], aes(x=N,y=theta,fill=simulation,color=simulation)) + geom_boxplot() +
                        labs(y="Watterson's theta", color="Simulation",fill="Simulation") +
                        facet_wrap(~T, ncol=1, labeller=label_both,scales="free") +
                        scale_color_manual(values=c("#5b0000", "deepskyblue4"),labels=c("coalescent","forward")) +
                        scale_fill_manual(values=c("#8e1a00","deepskyblue3"),labels=c("coalescent","forward")) +
                        theme(text=element_text(size=55),axis.text=element_text(size=30),axis.text.x=element_text(angle=45,margin=margin(t=20))) +
                        labs(title="n = 20")
            png("n20-distribution-comparisons-theta-nobottleneck-msvssfs.png",height=3500,width=3450,res=300)
                bp_n
            dev.off()
            
library(ggforce)
        ## boxplot with n on x, theta on y
            bp<-list()
            pop_size_iter<-0
            theta_data_factors<-theta_data
            theta_data_factors$n<-as.factor(theta_data_factors$n)
            for (p in pop_sizes){
                pop_size_iter<-pop_size_iter+1
                bp[[pop_size_iter]]<-ggplot(theta_data_factors[which(theta_data_factors$N==p),], 
                                            aes(x=n,y=theta,fill=simulation)) + 
                    geom_violin() +
                    labs(y="Watterson's theta",title=paste("N =",p))
            }
            
           png("N-distribution-comparisons-theta-msvssfs-nobottleneck.png",height=1000,width=2000)
                multiplot(plotlist=bp,cols=3)
            dev.off()
            
library(scales)            
## heat map of differences
    # sfs vs slim
        sfs_slim<-matrix(100*((unlist(theta_means_slim)-unlist(theta_means_sfs))/unlist(theta_means_sfs))[seq(1,126,2)],nrow=9,ncol=7)
            colnames(sfs_slim)<-c(10,20,50,100,200,500,1000)
            rownames(sfs_slim)<-c(2,5,10,20,50,100,200,500,1000)
            sfs_slim_melt<-melt(sfs_slim)
                sfs_slim_melt<-cbind(sfs_slim_melt,unlist(ks.p.adj.list.slim))
                colnames(sfs_slim_melt)<-c("n","Ne","difference","fdr")
                sfs_slim_melt$n<-as.factor(sfs_slim_melt$n)
                sfs_slim_melt$Ne<-as.factor(sfs_slim_melt$Ne)
                sfs_slim_melt$difference[which(sfs_slim_melt$fdr>=0.001)]<-NA
            sfs_slim_hm<-ggplot(sfs_slim_melt,aes(x=Ne,y=n,fill=difference)) +
                geom_tile() +
                xlab(label=expression("N"["e"])) +
                ylab(label="n") +
                scale_fill_gradientn(name="percent\ndifference",
                                     limits=c(min(na.omit(sfs_slim_melt$difference)),-min(na.omit(sfs_slim_melt$difference))),
                                     colors=c("blue","dodgerblue","white","tan1","darkorange1"),
                                     values=rescale(c(min(na.omit(sfs_slim_melt$difference)),-2,0,2,-min(na.omit(sfs_slim_melt$difference))))) +

                theme_bw()
            sfs_slim_hm
            ggsave("heatmap_sfs_slim.png",sfs_slim_hm)
            ggsave("heatmap_sfs_slim.eps",sfs_slim_hm)
    # sfs vs ms
        sfs_ms<-matrix(100*((unlist(theta_means_ms)-unlist(theta_means_sfs))/unlist(theta_means_sfs))[seq(1,126,2)],nrow=9,ncol=7)
            colnames(sfs_ms)<-c(10,20,50,100,200,500,1000)
            rownames(sfs_ms)<-c(2,5,10,20,50,100,200,500,1000)
            sfs_ms_melt<-melt(sfs_ms)
                sfs_ms_melt<-cbind(sfs_ms_melt,unlist(ks.p.adj.list))
                colnames(sfs_ms_melt)<-c("n","Ne","difference","fdr")
                sfs_ms_melt$n<-as.factor(sfs_ms_melt$n)
                sfs_ms_melt$Ne<-as.factor(sfs_ms_melt$Ne)
                sfs_ms_melt$difference[which(sfs_ms_melt$fdr>=0.001)]<-NA
            range<-max(na.omit(sfs_ms_melt$difference))-min(na.omit(sfs_ms_melt$difference))
            sfs_ms_hm<-ggplot(sfs_ms_melt,aes(x=Ne,y=n,fill=difference)) +
                geom_tile() +
                xlab(label=expression("N"["e"])) +
                ylab(label="n") +
                scale_fill_gradientn(name="percent\ndifference",
                                     limits=c(-max(na.omit(sfs_ms_melt$difference)),max(na.omit(sfs_ms_melt$difference))),
                                     colors=c("blue","dodgerblue","white","tan1","darkorange1"),
                                     values=rescale(c(-max(na.omit(sfs_ms_melt$difference)),-25,0,25,max(na.omit(sfs_ms_melt$difference))))) +
                theme_bw()
            sfs_ms_hm
            ggsave("heatmap_sfs_ms.png",sfs_ms_hm)
            ggsave("heatmap_sfs_ms.eps",sfs_ms_hm)

        ## n/Ne
            # time_iter<-0
            # new_nNe_diff_theta<-list()
            # for (time in times){
            #     time_iter<-time_iter+1
            #     new_nNe_diff<-melt(n_over_Ne_diff,id="n_over_Ne",na.rm=TRUE)
            #      new_nNe_diff<-new_nNe_diff[,-2]
            #      colnames(new_nNe_diff)<-c("n_over_Ne","perc_diff")
            #      new_nNe_diff$Ne<-as.factor(c(rep(10,6),rep(20,6),rep(50,6),rep(100,6),rep(200,6),rep(500,6),rep(1000,6)))
            #      new_nNe_diff_theta<-new_nNe_diff
            # }
                new_nNe_diff<-melt(n_over_Ne_diff,id="n_over_Ne",na.rm=TRUE)
                 new_nNe_diff<-new_nNe_diff[,-2]
                 colnames(new_nNe_diff)<-c("n_over_Ne","perc_diff")
                 new_nNe_diff$Ne<-as.factor(c(rep(10,9),rep(20,9),rep(50,9),rep(100,9),rep(200,9),rep(500,9),rep(1000,9)))
                 new_nNe_diff_theta<-new_nNe_diff
            new_nNe_diff_theta<-cbind(new_nNe_diff_theta,T=rep("Constant Ne",length(new_nNe_diff_theta$Ne)))
            library(RColorBrewer)
                mycolours<-brewer.pal(n = 7, name = "Dark2")
            # nNe<-ggplot(data=new_nNe_diff_theta %>% arrange(desc(Ne)),aes(x=as.factor(n_over_Ne),y=perc_diff,fill=Ne,shape=Ne)) +
            #         geom_point(aes(color=Ne),size=6) +
            #         scale_shape_manual(values=rev(c(3,4,16,17,18,15,25)),name="Ne") +
            #         geom_hline(yintercept=0, linetype="dashed") +
            #         labs(x=expression(n/N[e]),y="Percent difference in theta") +
            #         scale_fill_manual(values = mycolours) +
            #         scale_color_manual(values=mycolours) +
            #         theme(panel.background = NULL,text=element_text(size=30),axis.text=element_text(size=18),axis.text.x = element_text(angle = 45, margin=margin(t=20))) +
            #         facet_wrap(~T, ncol=1,scales="free",labeller="label_both")
                nNe<-ggplot(data=new_nNe_diff_theta %>% arrange(desc(Ne)),
                            aes(x=as.factor(n_over_Ne),y=perc_diff,fill=Ne,shape=Ne)) +
                    geom_point(aes(color=Ne),size=7, alpha=.9) +
                    scale_shape_manual(values=rev(c(3,4,16,17,18,15,25)),name=expression(italic(N[e]))) +
                    geom_hline(yintercept=0, linetype="dashed") +
                    labs(x=expression(frac(italic(n),italic(N[e]))),y=expression(Percent~difference~"in"~theta[W])) +
                    #scale_color_grey(start=.6, end=0,name=expression(italic(N[e]))) +
                    #scale_fill_grey(start=.8, end=.2,name=expression(italic(N[e]))) +
                    scale_color_manual(values=mycolours,name=expression(italic(N[e]))) +
                    scale_fill_manual(values=mycolours,name=expression(italic(N[e]))) +
                    theme(panel.background = NULL,text=element_text(size=25),
                          axis.text=element_text(size=18),
                          axis.text.x = element_text(angle = 45, margin=margin(t=20)),
                          legend.title.align=0.2, legend.key = element_blank(),
                          legend.key.size = unit(2, 'lines')) 
                ggsave("n-Ne-shapes-nobottleneck.eps",nNe,height=12,width=13,device=cairo_ps)
            png("n-Ne-shapes-nobottleneck.png",height=2400,width=2600,res=300)
                nNe
            dev.off()
            
            
            # new_nNe_cv_sfs<-melt(n_over_Ne_cv_sfs,id="n_over_Ne",na.rm=TRUE)
            #     new_nNe_cv_sfs<-new_nNe_cv_sfs[,-2]
            # new_nNe_cv_ms<-melt(n_over_Ne_cv_ms,id="n_over_Ne",na.rm=TRUE)
            #     new_nNe_cv_ms<-new_nNe_cv_ms[,-2]
            #     new_nNe_cv<-rbind(new_nNe_cv_ms,new_nNe_cv_sfs)
            #     new_nNe_cv<-cbind(new_nNe_cv,c(rep("ms",length(new_nNe_cv_ms[,1])),rep("sfs",length(new_nNe_cv_sfs[,1]))))
            #     colnames(new_nNe_cv)<-c("n_over_Ne","CV","simulation")
                mycolours<-c("highlight" = "red", "normal" = "grey50")
                time_iter<-0
        # for (time in times){
        #     time_iter<-time_iter+1
        #     png(paste("T",time,"-n_over_Ne_theta-msvssfs.png",sep=""))
        #         nam<-paste("p",time_iter,sep="")
        #        p<-ggplot(data=new_nNe_diff_theta,aes(x=as.factor(n_over_Ne),y=perc_diff,fill=Ne)) + 
        #             geom_point(shape=21,size=3) + 
        #             geom_hline(yintercept=0, linetype="dashed") +
        #             #scale_fill_gradient(low="orange", high="blue") + 
        #             #scale_color_gradient(low="orange",high="blue") +
        #             labs(x="n/Ne",y="Percent difference in theta") +
        #             theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
        #             scale_color_manual("Ne", values = mycolours) +
        #             ggtitle(paste("T = ",time,sep=""))
        #        assign(nam,p)
        #        print(p)
        #     dev.off()
        # }
        # png("n_over_Ne_theta.png",width=2000,height=500)
        #         multiplot(p1,p2,p3,p4,cols=4)
        # dev.off()  
            # png("n_over_Ne_theta-cv.png")
            #     ggplot(data=new_nNe_cv,aes(x=as.factor(n_over_Ne),y=CV,fill=simulation,alpha=0.3)) + 
            #         geom_point(shape=21,size=3) + 
            #        # geom_hline(yintercept=0, linetype="dashed") +
            #        # scale_fill_gradient(low="orange", high="blue") + 
            #        # scale_color_gradient(low="orange",high="blue") +
            #         labs(x="n/Ne",y="Coefficient of Variation in theta")
            # dev.off()

# 
# Pairwise differences ----------------------------------------------------
   ## T-tests and distribution plots comparing sfs and ms - segregating sites
        pi_data<-rbind(ss_pi_ms[,c(1:2,4)],ss_pi_sfs[,c(1:2,4)],ss_pi_slim[,c(1:2,4)]) # column 1=N, 2=n, 3=ss, 4=pi
        pi_data<-cbind(pi_data,c(rep("ms",length(ss_pi_ms[,1])),rep("sfs",length(ss_pi_sfs[,1])),rep("slim",length(ss_pi_slim[,1]))))
        colnames(pi_data)<-c("N","n","pi","simulation")
        pi_data$N<-as.factor(pi_data$N)
        pi_data$n_over_Ne<-pi_data$n/as.numeric(as.character(pi_data$N))

    # ## Aggregate by l = 10 loci
    #     l <- 10
    #     pi_data_10<-data.frame()
    #  for (i in c(10,20,50,100,200,500,1000)){
    #     for (j in c(2,5,10,20,50,100)){
    #         for (s in c("ms","sfs")){
    #             temp_data<-as.numeric(pi_data[which(pi_data$N==i & pi_data$n==j & pi_data$simulation==s),3]) # pi counts for each set of parameters at each locus
    #             sums<-tapply(temp_data,rep(1:ceiling(length(temp_data)/l),each=l,length.out=length(temp_data)),sum) # pi counts for every 10 loci
    #             #sums<-sapply(temp_list,"[[",2) # unlist sums
    #             remainder<-if (length(temp_data)%%l !=0){1} else {0} # is there a remainder? if so it means one
    #             temp_data<-data.frame(N=rep(i,length(sums)-remainder),n=rep(j,length(sums)-remainder),seg_sites=sums[1:(length(sums)-remainder)],simulation=rep(s,length(sums)-remainder)) # subtract last set if there was a remainder (ie not a full 10 loci)
    #             pi_data_10<-rbind(pi_data_10,temp_data)
    #         }
    #     }
    #  }
    #     pi_data_10$N<-as.factor(pi_data_10$N)
    #     pi_data_10$n_over_Ne<-pi_data_10$n/as.numeric(as.character(pi_data_10$N))

    # par(mfrow=c(2,3))
        pi_means_sfs<-list()
        pi_means_ms<-list()
        pi_means_slim<-list()
        pi_anovas<-list()
        pi_cv_sfs<-list()
        pi_cv_ms<-list()
        pi_cv_slim<-list()
        pop_size_iter<-0
        for (pop_size in pop_sizes){
            pop_size_iter<-pop_size_iter+1
            sample_size_iter<-0
                pi_means_sfs[[pop_size_iter]]<-list()
                pi_means_ms[[pop_size_iter]]<-list()
                pi_means_slim[[pop_size_iter]]<-list()
                pi_anovas[[pop_size_iter]]<-list()
                pi_cv_sfs[[pop_size_iter]]<-list()
                pi_cv_ms[[pop_size_iter]]<-list()
                pi_cv_slim[[pop_size_iter]]<-list()
                for (sample_size in sample_sizes){
                    sample_size_iter<-sample_size_iter+1
                    pi_means_sfs[[pop_size_iter]][[sample_size_iter]]<-mean(pi_data[which(pi_data$n==sample_size & pi_data$N==pop_size & pi_data$simulation=="sfs"),3])
                    pi_means_sfs[[pop_size_iter]][[sample_size_iter]][2]<-sd(pi_data[which(pi_data$n==sample_size & pi_data$N==pop_size & pi_data$simulation=="sfs"),3])
                    pi_means_ms[[pop_size_iter]][[sample_size_iter]]<-mean(pi_data[which(pi_data$n==sample_size & pi_data$N==pop_size & pi_data$simulation=="ms"),3])
                    pi_means_ms[[pop_size_iter]][[sample_size_iter]][2]<-sd(pi_data[which(pi_data$n==sample_size & pi_data$N==pop_size & pi_data$simulation=="ms"),3])
                    pi_means_slim[[pop_size_iter]][[sample_size_iter]]<-mean(pi_data[which(pi_data$n==sample_size & pi_data$N==pop_size & pi_data$simulation=="slim"),3])
                    pi_means_slim[[pop_size_iter]][[sample_size_iter]][2]<-sd(pi_data[which(pi_data$n==sample_size & pi_data$N==pop_size & pi_data$simulation=="slim"),3])
                    pi_cv_sfs[[pop_size_iter]][[sample_size_iter]]<-(pi_means_sfs[[pop_size_iter]][[sample_size_iter]][2]/pi_means_sfs[[pop_size_iter]][[sample_size_iter]][1])*100
                    pi_cv_ms[[pop_size_iter]][[sample_size_iter]]<-(pi_means_ms[[pop_size_iter]][[sample_size_iter]][2]/pi_means_ms[[pop_size_iter]][[sample_size_iter]][1])*100
                    pi_cv_slim[[pop_size_iter]][[sample_size_iter]]<-(pi_means_slim[[pop_size_iter]][[sample_size_iter]][2]/pi_means_slim[[pop_size_iter]][[sample_size_iter]][1])*100
                    pi_anovas[[pop_size_iter]][[sample_size_iter]]<-tryCatch(t.test(pi_data[which(pi_data$n==sample_size & pi_data$N==pop_size),3]~pi_data[which(pi_data$n==sample_size & pi_data$N==pop_size),4]), error=function(e) NA)
                }
            }

            ## adjusted p-values
                pvalues<-unlist(unlist(unlist(pi_anovas,recursive=F),recursive=F),recursive=F)
                    p.adj<-p.adjust(unlist(pvalues[which(names(pvalues)=="Pr(>F)")])[which(names(unlist(pvalues[which(names(pvalues)=="Pr(>F)")]))=="Pr(>F)1")],method="fdr")
                    p.adj.list<-relist(unname(p.adj),skeleton=rep(list(rep("a",6)),7))
                pvaluesforward<-unlist(lapply(unname(unlist(unlist(pi_tukey_forward,recursive=F),recursive=F)),unname))[seq(4,6*7*4,4)]
                    p.adj.forward<-p.adjust(pvaluesforward,method="fdr")
                    p.adj.forward.list<-relist(unname(p.adj.forward),skeleton=rep(list(rep("a",6)),7))
                pvaluescoalescent<-unlist(lapply(unname(unlist(unlist(pi_tukey_ms_sfs,recursive=F),recursive=F)),unname))[seq(4,6*7*4,4)]
                    p.adj.coalescent<-p.adjust(pvaluescoalescent,method="fdr")
                    p.adj.coalescent.list<-relist(unname(p.adj.coalescent),skeleton=rep(list(rep("a",6)),7))
                kspvalues<-unlist(ks_result_p)
                ks.p.adj<-p.adjust(kspvalues,method="fdr")
                ks.p.adj.list<-relist(unname(ks.p.adj,skeleton=rep(list(rep("a",6)),7)))
        ## Density plots of each distribution

        palette<-c("grey","orange","blue","turquoise")
        dp<-list()
        vp<-list()
        n_over_Ne_diff<-data.frame(n_over_Ne=unique(pi_data$n_over_Ne)[sort.list(unique(pi_data$n_over_Ne))],check.names=FALSE)
        namevector<-rep("pi.diff",7)
        n_over_Ne_cv_ms<-data.frame(n_over_Ne=unique(pi_data$n_over_Ne)[sort.list(unique(pi_data$n_over_Ne))],check.names=FALSE)
        n_over_Ne_cv_sfs<-data.frame(n_over_Ne=unique(pi_data$n_over_Ne)[sort.list(unique(pi_data$n_over_Ne))],check.names=FALSE)
        namevector_cv<-rep("pi.cv",7)
        n_over_Ne_diff[,namevector]<-NA
        n_over_Ne_cv_ms[,namevector_cv]<-NA
        n_over_Ne_cv_sfs[,namevector_cv]<-NA
        pop_size_iter<-0
        for (pop_size in pop_sizes){
            pop_size_iter<-pop_size_iter+1
            sample_size_iter<-0
            ## density plots
                dp[[pop_size_iter]]<-list()
                vp[[pop_size_iter]]<-list()
                for (sample_size in sample_sizes){
                    sample_size_iter<-sample_size_iter+1
                    means<-ddply(pi_data[which(pi_data$n==sample_size & pi_data$N==pop_size),],"simulation",summarise,seg_sites.mean=mean(seg_sites))
                   # cvs<-data.frame(simulation=c("ms","sfs","slim","exp"),seg_sites.cv="NA")
                    # levels(means$simulation)<-c(levels(means$simulation),"expected")
                        # means<-rbind(means,c("expected",pi_exp_df[[pop_size_iter]][l,2]))
                        # means$seg_sites.mean<-as.numeric(means$seg_sites.mean)
                        # levels(pi_data$simulation)<-c(levels(pi_data$simulation),"expected")
                        # pi_data<-rbind(pi_data,c("10",2,0,"expected"))
                        #     pi_data$seg_sites<-as.numeric(pi_data$seg_sites)
                    exp<-data.frame("exp",watt_exp_df[[pop_size_iter]][sample_size_iter,2])
                    names(exp)<-c("simulation","seg_sites.mean")
                    means<-rbind(means,exp)
                    nNe<-sample_size/as.numeric(as.character(pop_size))
                    perc_diff<-((means[which(means$simulation=="ms"),2]-means[which(means$simulation=="sfs"),2])/means[which(means$simulation=="sfs"),2])*100
                    perc_diff_forward<-((means[which(means$simulation=="sfs"),2]-means[which(means$simulation=="slim"),2])/means[which(means$simulation=="sfs"),2])*100
                    n_over_Ne_diff[which(n_over_Ne_diff$n_over_Ne==nNe),pop_size_iter+1]<-perc_diff
                    n_over_Ne_cv_ms[which(n_over_Ne_cv_ms$n_over_Ne==nNe),pop_size_iter+1]<-pi_cv_ms[[pop_size_iter]][[sample_size_iter]]
                    n_over_Ne_cv_sfs[which(n_over_Ne_cv_sfs$n_over_Ne==nNe),pop_size_iter+1]<-pi_cv_sfs[[pop_size_iter]][[sample_size_iter]]
                    dp[[pop_size_iter]][[sample_size_iter]]<-ggplot(pi_data[which(pi_data$n==sample_size & pi_data$N==pop_size),], aes(x=seg_sites, fill=simulation)) + geom_density(alpha=.4) +
                        geom_vline(data=means[1:2,],aes(xintercept=seg_sites.mean,color=simulation),linetype="dashed",size=1) +
                        #geom_vline(aes(xintercept=as.numeric(pop_size)*4*0.000001*10), color="#BB0000") +
                        #geom_vline(aes(xintercept=pi_means_sfs[[pop_size_iter]][[sample_size_iter]][1],color="mean sfs"),linetype="dashed",size=1) +
                        geom_point(data=means,aes(x=seg_sites.mean,y=0,color=simulation),alpha=1,size=2) +
                        labs(x="Segregating Sites",title=paste("Ne = ",pop_size,", n = ",sample_size,sep=""),subtitle=paste("Difference between means, FDR = ",signif(p.adj.list[[pop_size_iter]][[sample_size_iter]],3),"\n",signif(perc_diff,3),"% difference between sfs and ms means\n",signif(perc_diff_forward,3),"% difference between forward means",sep="")) +
                        scale_fill_manual(values=palette) +
                        scale_colour_manual(values=palette)
                    vp[[pop_size_iter]][[sample_size_iter]]<-ggplot(pi_data[which(pi_data$n==sample_size & pi_data$N==pop_size),], aes(x=simulation, y=seg_sites)) +
                        geom_violin(aes(fill = simulation)) + 
                        geom_hline(data=means[1:2],aes(yintercept=seg_sites.mean,color=simulation),linetype="dashed") +
                        labs(x="Simulation",y="Segregating Sites",title=paste("Ne = ",pop_size,", n = ",sample_size,sep=""),subtitle=paste("Difference between means, FDR = ",signif(p.adj.list[[pop_size_iter]][[sample_size_iter]],3),"\n",signif(perc_diff,3),"% difference between means",sep="")) +
                        scale_fill_manual(values=palette) +
                        scale_colour_manual(values=palette)

                    dp[[pop_size_iter]][[sample_size_iter]]
                    vp[[pop_size_iter]][[sample_size_iter]]
                }
            }
# 
#             png("distribution-comparisons-pi-msvssfs.png",height=3000,width=3000)
#                 multiplot(plotlist=unlist(dp,recursive=FALSE),cols=6)
#             dev.off()
# 
#             png("distribution-comparisions-violin-pi-msvssfs.png",height=3000,width=3000)
#                 multiplot(plotlist=unlist(vp,recursive=FALSE),cols=6)
#             dev.off()
# 
#             ## n/Ne
#             new_nNe_diff<-melt(n_over_Ne_diff,id="n_over_Ne",na.rm=TRUE)
#                 new_nNe_diff<-new_nNe_diff[,-2]
#                 colnames(new_nNe_diff)<-c("n_over_Ne","perc_diff")
#             png("n_over_Ne_pi-msvssfs.png")
#                 ggplot(data=new_nNe_diff,aes(x=as.factor(n_over_Ne),y=perc_diff,fill=perc_diff)) +
#                     geom_point(shape=21,size=3) +
#                     scale_fill_gradient(low="orange", high="blue") +
#                     scale_color_gradient(low="orange",high="blue") +
#                     labs(x="n/Ne",y="Percent difference in mean pairwise differences")
#             dev.off()
# 
#  
# 
# 
# save.image()
# 
# ## Individual plots - CAREFUL with which version of dp you use (will be the most recent set of plots run, ss, theta, or pi)
png("cv-msvssfs.png")
    cv_frame<-data.frame(N=rep(c(10,20,50,100,200,500,1000),6),n=c(rep(2,7),rep(5,7),rep(10,7),rep(20,7),rep(50,7),rep(100,7)),simulation=c(rep("sfs",42),rep("ms",42)),sumstat=rep(rep("wtheta",42),2))
        cv_frame_temp<-c(unlist(unlist(theta_cv_sfs,recursive=F),recursive=F),unlist(unlist(theta_cv_ms,recursive=F),recursive=F))
        cv_frame<-cbind(cv_frame,cv_frame_temp)
            colnames(cv_frame)<-c(colnames(cv_frame)[-5],"cv")
    ggplot(data=cv_frame[which(cv_frame$n==10),],aes(x=N,y=cv,group=simulation)) +
        geom_point(aes(color=simulation)) + geom_line(aes(color=simulation), linetype="dashed")
dev.off()
#             
#     # Ne=1000,n=50
#         png("plots/theta-1000-50.png",height=600,width=800)
#             dp[[5]][[7]] + theme_classic(base_size = 30, base_family = "")
#         dev.off()
#     # Ne=1000,n=2
#         png("plots/theta-1000-2.png",height=600,width=800)
#             dp[[1]][[7]] + theme_classic(base_size = 30, base_family = "")
#         dev.off()
#     # Ne=100,n=50
#         png("plots/theta-100-50.png",height=600,width=800)
#             dp[[5]][[4]] + theme_classic(base_size = 30, base_family = "")
#         dev.off()
#     # Ne=100,n=2
#         png("plots/theta-100-2.png",height=600,width=800)
#             dp[[1]][[4]] + theme_classic(base_size = 30, base_family = "")
#         dev.off()
#     # Ne=10,n=50
#         png("plots/theta-10-50.png",height=600,width=800)
#             dp[[5]][[1]] + theme_classic(base_size = 30, base_family = "")
#         dev.off()
#     # Ne=10,n=2
#         png("plots/theta-10-2.png",height=600,width=800)
#             dp[[1]][[1]] + theme_classic(base_size = 30, base_family = "")
#         dev.off()
#     # n_over_Ne_theta
#             new_nNe_diff_theta$n_over_Ne[which(new_nNe_diff_theta$n_over_Ne<.1)]<-"n<<Ne"
#             new_nNe_diff_theta$n_over_Ne[which(new_nNe_diff_theta$n_over_Ne<1 & new_nNe_diff_theta$n_over_Ne>=.1)]<-"n<Ne"
#             new_nNe_diff_theta$n_over_Ne[which(new_nNe_diff_theta$n_over_Ne==1)]<-"n=Ne"
#             new_nNe_diff_theta$n_over_Ne[which(new_nNe_diff_theta$n_over_Ne>1)]<-"n>Ne"
#             new_nNe_diff_theta$n_over_Ne<-as.factor(new_nNe_diff_theta$n_over_Ne)
#         png("plots/n_over_Ne_theta.png",height=600,width=800)
#             ggplot(data=new_nNe_diff,aes(x=as.factor(n_over_Ne),y=perc_diff,fill=perc_diff)) + 
#                     geom_point(shape=21,size=6) + 
#                     geom_hline(yintercept=0, linetype="dashed") +
#                     scale_fill_gradient(low="orange", high="blue", name="Percent\nDifference") + 
#                     scale_color_gradient(low="orange",high="blue") +
#                     labs(x="n/Ne",y="Percent difference in mean theta") + theme_classic(base_size = 30, base_family = "") +
#                     theme(text = element_text(size=30), axis.text.x = element_text(angle=90, hjust=1))
#             dev.off()
#     # CV n=2
#         png("plots/cv-n2.png",height=600,width=800)
#             cv_frame<-data.frame(N=rep(c(10,20,50,100,200,500,1000),5),n=c(rep(2,7),rep(5,7),rep(10,7),rep(20,7),rep(50,7)),simulation=c(rep("sfs",35),rep("ms",35)),sumstat=rep(rep("wtheta",35),2))
#             cv_frame_temp<-c(unlist(unlist(theta_cv_sfs,recursive=F),recursive=F),unlist(unlist(theta_cv_ms,recursive=F),recursive=F))
#             cv_frame<-cbind(cv_frame,cv_frame_temp)
#                 colnames(cv_frame)<-c(colnames(cv_frame)[-5],"cv")
#             ggplot(data=cv_frame[which(cv_frame$n==2),],aes(x=N,y=cv,group=simulation)) +
#                 geom_point(aes(color=simulation),size=6) + geom_line(aes(color=simulation), linetype="dashed") +
#                 theme_classic(base_size = 30, base_family = "")
#         dev.off()
#     # CV n=10
#         png("plots/cv-n10.png",height=600,width=800)
#             cv_frame<-data.frame(N=rep(c(10,20,50,100,200,500,1000),5),n=c(rep(2,7),rep(5,7),rep(10,7),rep(20,7),rep(50,7)),simulation=c(rep("sfs",35),rep("ms",35)),sumstat=rep(rep("wtheta",35),2))
#             cv_frame_temp<-c(unlist(unlist(theta_cv_sfs,recursive=F),recursive=F),unlist(unlist(theta_cv_ms,recursive=F),recursive=F))
#             cv_frame<-cbind(cv_frame,cv_frame_temp)
#                 colnames(cv_frame)<-c(colnames(cv_frame)[-5],"cv")
#             ggplot(data=cv_frame[which(cv_frame$n==10),],aes(x=N,y=cv,group=simulation)) +
#                 geom_point(aes(color=simulation),size=6) + geom_line(aes(color=simulation), linetype="dashed") +
#                 theme_classic(base_size = 30, base_family = "")
#         dev.off()
#     # CV n=50
#         png("plots/cv-n50.png",height=600,width=800)
#             cv_frame<-data.frame(N=rep(c(10,20,50,100,200,500,1000),5),n=c(rep(2,7),rep(5,7),rep(10,7),rep(20,7),rep(50,7)),simulation=c(rep("sfs",35),rep("ms",35)),sumstat=rep(rep("wtheta",35),2))
#             cv_frame_temp<-c(unlist(unlist(theta_cv_sfs,recursive=F),recursive=F),unlist(unlist(theta_cv_ms,recursive=F),recursive=F))
#             cv_frame<-cbind(cv_frame,cv_frame_temp)
#                 colnames(cv_frame)<-c(colnames(cv_frame)[-5],"cv")
#             ggplot(data=cv_frame[which(cv_frame$n==50),],aes(x=N,y=cv,group=simulation)) +
#                 geom_point(aes(color=simulation),size=6) + geom_line(aes(color=simulation), linetype="dashed") +
#                 theme_classic(base_size = 30, base_family = "")
#         dev.off()