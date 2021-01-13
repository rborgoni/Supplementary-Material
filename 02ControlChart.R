rm(list=ls()); ls()
### libraries #####
require(spatstat); require(boot); require(ggplot2)
set.seed(123)
##############  IMPORT DATA  ####

dd <- read.table("data.csv", sep=";",header=TRUE)

##############  CONTROL CHART ####
# pvalue and confidence interval 
##############  sub-function  ####
control_chart <- function(lag, numero_sim=500, input_data) {
  pval <- function(data,indices) { #function for p-value bootstrapping
    d <- data[indices]
    pvalue <- round(mean(d>D_oss),3)
    return(pvalue)
  }
  set.seed(123)
  w <- disc(1,c(0,0))
  X <- as.ppp(input_data[which(input_data$lag<=lag),c("X","Y")], W=w) 
  GGb     <- Gest(X, r=NULL, breaks = NULL,correction = "none")
  envpp   <- spatstat::envelope(X, fun=Gest, nsim = numero_sim, verbose = F, savefuns =T, savepatterns = T)
  simfuns <- attr(envpp, "simfuns", exact = T)
  G_star  <- as.data.frame(simfuns)  
  F_theo  <- GGb$theo
  F_oss   <- GGb$raw
  D_oss   <- sum((F_theo-F_oss)^2)
  D_star  <- numeric()
  for (i in 1:(ncol(G_star)-1)) {
    D_star[i] <- sum((F_theo-G_star[,i+1])^2)
    }
  output <- list()
  pvalue <- round(mean(D_star>D_oss),3) #pvalue MC
  bootobject    <- boot(D_star, statistic=pval,R=1000) #ic bootstrap
  bootCI <- boot.ci(bootobject, type = "perc")
  output[[1]]   <- pvalue; output[[2]] <- bootCI
  names(output) <- c("pvalue", "bootCI")
  output
}
### End sub-function 
#
startt <- Sys.time() 
r=1000   # number of MC replicates
c1 <- control_chart(1,r, dd) 
c2 <- control_chart(2,r, dd) 
c3 <- control_chart(3,r, dd) 
c4 <- control_chart(4,r, dd) 
c5 <- control_chart(5,r, dd) 

endd <- Sys.time()
print(paste("execution time: ",endd -startt))
#
#prepare plot
punti <- data.frame(x=c(1,2,3,4,5), y=c(c1$pvalue, c2$pvalue, c3$pvalue, c4$pvalue, c5$pvalue)); 
punti <- cbind(punti,rbind(c1$bootCI$percent[4:5],
                           c2$bootCI$percent[4:5],
                           c3$bootCI$percent[4:5],
                           c4$bootCI$percent[4:5],
                           c5$bootCI$percent[4:5])); 
names(punti) <- c("lag","pvalue","lower 0.95","upper 0.95")


punti=punti[-nrow(punti),]
plot <- ggplot(punti, aes(x=lag, y=pvalue, group=1))+geom_point(col="red",size=3)+ylim(0,0.25)+
  geom_text(aes(label=pvalue, hjust=-.5))+
  geom_line(size=1.1)+ scale_linetype_identity()+
  geom_segment(aes(x=1,xend=5,y=0,yend=0)) + 
  geom_segment(aes(x=1,xend=5,y=0.1, yend=0.1,linetype="dotdash")) + 
  geom_segment(aes(x=1,xend=5,y=0.05,yend=0.05,linetype="dashed")) + 
  geom_segment(aes(x=1,xend=5,y=0.01,yend=0.01,linetype="twodash"))+
  annotate("text",x=1,y=0.11,label="0.1",size=3)+
  annotate("text",x=1,y=0.06,label="0.05",size=3)+
  annotate("text",x=1,y=0.02,label="0.01",size=3)+
  geom_segment(aes(x=1,   xend=1,y= `lower 0.95`[1],yend= `upper 0.95`[1]))+
  geom_segment(aes(x=2,   xend=2,y= `lower 0.95`[2],yend= `upper 0.95`[2]))+
  geom_segment(aes(x=3,   xend=3,y= `lower 0.95`[3],yend= `upper 0.95`[3]))+
  geom_segment(aes(x=4,   xend=4,y= `lower 0.95`[4],yend= `upper 0.95`[4]))+
  geom_segment(aes(x=0.95,xend=1.05,y= `lower 0.95`[1],yend= `lower 0.95`[1]))+
  geom_segment(aes(x=0.95,xend=1.05,y= `upper 0.95`[1],yend= `upper 0.95`[1]))+
  geom_segment(aes(x=1.95,xend=2.05,y= `lower 0.95`[2],yend= `lower 0.95`[2]))+
  geom_segment(aes(x=1.95,xend=2.05,y= `upper 0.95`[2],yend= `upper 0.95`[2]))+
  geom_segment(aes(x=2.95,xend=3.05,y= `lower 0.95`[3],yend= `lower 0.95`[3]))+
  geom_segment(aes(x=2.95,xend=3.05,y= `upper 0.95`[3],yend= `upper 0.95`[3]))
plot
