rm(list=ls()); ls()
### libraries #####
require(spatstat); require(igraph); require(alphahull); 
#
##############  IMPORT DATA  ####
ok=read.table("clusters.csv", header = T, sep=";")
#
##############  sub-function    ####
ashapeToOwin <- function (ashape){ 
  # alpha-shape to matrix - to be use in owin
  alpha <- ashape  
  n10g_21 = graph.edgelist(cbind(as.character(alpha$edges[, "ind1"]), as.character(alpha$edges[,"ind2"])), directed = FALSE)
  if (!is.connected(n10g_21)) return(F)
  if (any(degree(n10g_21) != 2)) return(F)
  if (clusters(n10g_21)$no > 1) return(F)
  cutg_21 = n10g_21 - E(n10g_21)[1]
  ends_21 = names(which(degree(cutg_21) == 1))
  path_21 = get.shortest.paths(cutg_21, ends_21[1], ends_21[2])[[1]]
  pathX_21 = as.numeric(V(n10g_21)[path_21[[1]]]$name)
  pathX_21 = c(pathX_21, pathX_21[1])
  obj_21 <- alpha$x[pathX_21, ]
  return(obj_21)
}
### End sub-function 
##############  MAIN START  ####
mindefect=11
# select largest clusters
ok$id <- as.character(ok$id); 
group <- names(table(ok[ok$cluster!=0 & ok$Freq>mindefect, "id"]))
ncp=length(group); ncp
#
### Alphashape construction for each cluster #####
win <- disc(1,c(0,0))
result=data.frame()
lalpa=l=NULL
  plot(win, lwd=3,main="")
  for (i in 1:ncp)  {
    m <- as.matrix(ok[ok$id==group[i],c("X", "Y")]) 
      points(m, pch=i, col=i)
    w<-convexhull.xy(x=m[,"X"],y=m[,"Y"])
    k=max(dist(m))
    alpha0 = ashape(m[,"X"], m[,"Y"], alpha = k^2); #plot(alpha0)
    v=sort(unique(as.vector(alpha0$edges[, c("ind1","ind2")])))
    check=T
    a1=aa=k;     AA=area(w);     count=conta=1;     eps=0.00001;     a=k/2;     convergence=F
    while(check){
      count=count+1
      alpha = ashape(m[,"X"], m[,"Y"], alpha = a) 
      v1=sort(unique(as.vector(alpha$edges[, c("ind1","ind2")])))
      CK.CNT=length(setdiff(v,v1))==0
      if(CK.CNT) { 
        conta=conta+1
        objCL1 = ashapeToOwin(alpha)
        if(is.matrix(objCL1)){
          aa=c(aa,a)
          a1=a
          a=a/2
          winCL1 <- try(owin(poly=data.frame(x=(objCL1[,1]),  y=(objCL1[,2]))),silent=TRUE)
          if(class(winCL1) == "try-error") 
            winCL1 <- owin(poly=data.frame(x=rev(objCL1[,1]),y=rev(objCL1[,2])))
          A=area(winCL1)
          AA=c(AA,A)
         if( (abs(A-AA[(length(AA)-1)])/AA[(length(AA)-1)]) <eps) {check=F; convergence=T
            }  else {a=a+(a1-a)/2 }
        }  else {check=F }
      } else {check=F }
    }
    print(paste("alfa = ", aa[length(aa)]))
    result=rbind(result,c(count=count,conta=conta,convergence=convergence,eps=eps,alfa=aa[length(aa)]))
    names(result)=c("count","conta","convergence","eps","alfa")
    alpha <- ashape(m[,"X"], m[,"Y"], alpha = aa[length(aa)])    
     plot(alpha,lwd=2, add=T);     points(ok[!ok$id %in% group,c("X", "Y")], pch=1, col="lightgray",cex=0.5)
    aaa1 = ashapeToOwin(alpha) #  area 
    aaa <- try(owin(poly=data.frame(x=(aaa1[,1]),y=(aaa1[,2]))),silent=TRUE)
    if(class(aaa) == "try-error") 
      aaa <- owin(poly=data.frame(x=rev(aaa1[,1]),y=rev(aaa1[,2])))
    l=c(l,  area(aaa) )
    lalpa[[i]]=aaa
  }
  # label the clusters
  text(apply(ok[ok$id==group[1],c("X", "Y")],2,"mean")[1]*3,   apply(ok[ok$id==group[1],c("X", "Y")],2,"mean")[2], "A",cex=1.5)
  text(apply(ok[ok$id==group[4],c("X", "Y")],2,"mean")[1]*0.9, apply(ok[ok$id==group[4],c("X", "Y")],2,"mean")[2], "B",cex=1.5,col="blue")
  text(apply(ok[ok$id==group[2],c("X", "Y")],2,"mean")[1]*0.9, apply(ok[ok$id==group[2],c("X", "Y")],2,"mean")[2], "C",cex=1.5,col="red")
  text(apply(ok[ok$id==group[3],c("X", "Y")],2,"mean")[1]*0.6, apply(ok[ok$id==group[3],c("X", "Y")],2,"mean")[2], "D",cex=1.5,col="green")
  