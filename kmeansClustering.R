# Initialize: remove files, include libraries
options(warn=-1)
invisible(if (file.exists("BTplot.jpg")) file.remove("BTplot.jpg"))
invisible(if (file.exists("heatmap.jpg")) file.remove("heatmap.jpg"))
invisible(if (file.exists("beeswarm.jpg")) file.remove("beeswarm.jpg"))
invisible(if (file.exists("ecdf.jpg")) file.remove("ecdf.jpg"))
invisible(if (file.exists("stats.txt")) file.remove("stats.txt"))
invisible(if (file.exists("clusters.txt")) file.remove("clusters.txt"))
invisible(install.packages("ggplot2",repos = "http://cran.us.r-project.org"))
invisible(install.packages("beeswarm",repos = "http://cran.us.r-project.org"))
invisible(install.packages("phangorn",repos = "http://cran.us.r-project.org"))
library("ggplot2")
library('beeswarm')
library("phangorn")

args = commandArgs(trailingOnly=TRUE)
if (length(args)!=5) {
  stop("Wrong parameters. Run as: Rscript kmeansClustering.R <distance matrix> <k1> <k2> <k*> <cut>", call.=FALSE)
}

# Read input file
mx <- as.matrix(read.table(args[1],header=T,row.names=1))

k1 <- as.numeric(args[2])
k2 <- as.numeric(args[3])
kx <- as.numeric(args[4])
cut <- as.numeric(args[5])

b <- c()
t <- c()
bt <- c()

message("Calculating k-means...")
for (i in k1:k2) {
  res <- kmeans(mx,i)
  b <- c(b,res$betweenss)
  t <- c(t,res$totss)
}

bt <- b/t
btm <- max(bt)
i_max_bt <- which.max(bt)
res_max <- kmeans(mx,i_max_bt)

# If k* is set
if (kx > 0) {
  i_max_bt <- kx
  res_max <- kmeans(mx,i_max_bt)
}

######################## Plots ###########################

message("Drawing plots...")

# Betweenness and withinness plots
jpeg("BTplot.jpg")
plot(bt,type="line",col="green",xlab="k",ylab="",axes=F, ylim=c(0,1))
axis(side = 1, at = c(1:k2), cex.axis=0.65)
axis(side = 2)
legend("topleft",c("betw./tot.wi"),col=c("green"),lty=c(1,1,1))
invisible(dev.off())

# Heatmap
species <- row.names(mx)
myBreaks <- c(seq(0,1,by=0.01))
cexx = 1-length(species)/200
ceyy = cexx

# normalize mx to 0-1
mx_orig <- mx
mx2 <- (mx - min(mx))/(max(mx) - min(mx))
mx <- mx2

jpeg(filename = "heatmap.jpg", height = 5000, width = 5000, units = "px", res = 600)
heatmap(mx, symkey =F, symbreaks=F, scale="none", dendrogram = F, Rowv=F, Colv=F,col = topo.colors(100), breaks = myBreaks, na.color="white", margin = c(12,16), 
        cexRow=cexx,cexCol=ceyy, key=T, trace="none", lmat=rbind( c(4, 3), c(2,1), c(0,0) ), lhei=c(1.8,6.5, 1 ), dendrogram="none")
invisible(dev.off())

# Beeswarm plot
dvals <- mx[upper.tri(mx)]
jpeg(filename = "beeswarm.jpg")
beeswarm(dvals,pch=16,cex=0.5,ylim=c(0,1),xlab="",ylab="Similarity")
invisible(dev.off())

# ECDF plot
jpeg(filename = "ecdf.jpg")
dvals_ecdf <- ecdf(dvals)
plot(dvals_ecdf,xlim=c(0,1),cex=0.5,xlab="Similarity",ylab="CDF",main="CDF plot as a function of similarity")
invisible(dev.off())

######################### Stats files #########################

# stats and clusters files

write.table(res_max$cluster, file="clusters.txt", col.names=F, quote=F, sep="\t")

header = "baramin\tspecies\tmean\tstdev\tmin\tmax\tp-value\tneglog"
write(header, file="stats.txt", sep="\t", append=T)

mx <- mx_orig
cluster_sizes = res_max$size
for (n_cluster in 1:i_max_bt) {
  csize = cluster_sizes[n_cluster]
  if (csize >= 3) {
    m1 = as.matrix(mx[res_max$cluster == n_cluster,res_max$cluster == n_cluster])
    m1_2 = as.matrix(mx2[res_max$cluster == n_cluster,res_max$cluster == n_cluster])
    
    # draw cladogram
    tree_upgma <- upgma(m1_2)
    treename <- paste("cluster",n_cluster,".jpg",sep="")
    jpeg(filename = treename)
    plot(tree_upgma, main=paste("cluster",n_cluster,sep=" "))
    dev.off()
    
    x = m1[upper.tri(m1)]
    ll = dim(m1)[1]
    
    m2 = as.matrix(cbind(mx[res_max$cluster != n_cluster,res_max$cluster == n_cluster],t(mx[res_max$cluster == n_cluster,res_max$cluster != n_cluster])))
    m2b = m2[!duplicated(colnames(m2))]
    
    t = t.test(x,m2b)
    pval = t$p.value
    nglog = -log10(pval)
    min = min(x)
    max = max(x)
    
    mean2 = sprintf("%.3f", mean(x))
    sd2 = sprintf("%.3f", sd(x))
    min2 = sprintf("%.3f", min)
    max2 = sprintf("%.3f", max)
    pval2 = sprintf("%.3f", pval)
    nglog2 = sprintf("%.3f", nglog)
    
    stats = paste(n_cluster, ll, mean2, sd2, min2, max2, pval, nglog2, sep="\t")
    stats2 = gsub("\n","\t",stats)
    write(stats, file="stats.txt", sep="\t", append = T)
  }
}

# Estimation of similarity cutoff value
sdvals <- sort(dvals)
mdvals <- c()

w <- 50
inc <- 0.001
for(i in seq(from=0, to=1, by=inc)) {
  lw <- length(sdvals[sdvals>i & sdvals<=i+w*inc])
  mdvals <- c(mdvals, lw)
}

# Generate Cytoscape attr file
header <- paste("Species1","interaction","Species2","=","Score",sep="\t")
write(header, file="matrix.attr", sep="", append = T)
for (nr in row.names(mx)) {
  for (nc in colnames(mx)) {
    if ((!(nr == nc)) && (which(row.names(mx) == nr) < which(colnames(mx) == nc))) {
      lin <- paste(nr,"(pp)",nc,"=",mx[nr,nc],sep=" ")
      if (mx[nr,nc] >= cut) {
        write(lin, file="matrix.attr", sep="", append = T)
      }
    }
  }
}

message("Complete!")