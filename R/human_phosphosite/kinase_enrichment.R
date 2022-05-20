library(directPA)
library(calibrate)
library(pheatmap)

## The input is the processed phospho data
Tc <- as.matrix(grand.tab[, c(1:6)])
rownames(Tc) <- sapply(strsplit(rownames(Tc), ";"), function(x)paste(x[1], x[2], "", sep=";"))
##

## KinasePA function (edited)
perturbPlot2d_edit <- function(Tc, annotation, minSize=5, ...) {
  
  # step 1. convert statistics into z-scores
  Tc.zscores <- apply(Tc, 2, function(x){qnorm(rank(x)/(nrow(Tc)+1))})
  
  # step 2. filter the groups that are smaller than the minimun cutoff
  DE = lapply(annotation, function(x){
    if(sum(rownames(Tc.zscores) %in% x) >= minSize) {
      X <- Tc.zscores[rownames(Tc.zscores)%in%x,]
      n = nrow(X)
      Z1 = sum(X[,1])/sqrt(n)
      Z2 = sum(X[,2])/sqrt(n)
      list(Z1=Z1, Z2=Z2)
    }
  })
  
  # step3. filter DE that has 0 element
  DE <- DE[which(sapply(DE, length) != 0)]
  Z1 <- unlist(sapply(DE, function(x){x[1]}))
  Z2 <- unlist(sapply(DE, function(x){x[2]}))
  
  # visualization
  s <- sqrt(abs(Z1)+abs(Z2))
  plot(Z1, Z2, col="darkblue", pch=16, cex=1, xlab=colnames(Tc)[1], ylab=colnames(Tc)[2], ...)
  
  textxy(Z1, Z2, names(DE), col = "black", cex = 0.5)
  abline(v=0, h=0, col="gold", lty=2)
  abline(a=0, b=1, col="darkgreen", lty=2)
  abline(a=0, b=-1, col="darkgreen", lty=2)
  
  r <- ceiling(max(sqrt(Z1^2 + Z2^2)))
  for(i in seq(0, 8, 8/5)) {
    theta = seq(-3.14,3.14,0.05)
    lines(i*cos(theta),i*sin(theta),col = 'gray', type="l")
  }
  
  
  ## return the results
  result <- list()
  result$Z1 <- Z1
  result$Z2 <- Z2
  return(result)
}

z1 <- perturbPlot2d_edit(Tc=Tc[,c(1,2)], annotation=PhosphoSite.human,
                         xlim=c(-5, 6), ylim=c(-4,4), main="Kinase perturbation analysis")

z2 <- perturbPlot2d_edit(Tc=Tc[,c(3,4)], annotation=PhosphoSite.human,
                         xlim=c(-5, 6), ylim=c(-4,4), main="Kinase perturbation analysis")

z3 <- perturbPlot2d_edit(Tc=Tc[,c(5,6)], annotation=PhosphoSite.human,
                         xlim=c(-5, 6), ylim=c(-4,4), main="Kinase perturbation analysis")

Zs1 <- perturbPlot3d(Tc=Tc[,c(1,3,5)], annotation=PhosphoSite.human,
                     size=10, main="Kinase perturbation analysis")

Zs2 <- perturbPlot3d(Tc=Tc[,c(2,4,6)], annotation=PhosphoSite.human,
                     size=10, main="Kinase perturbation analysis")

# create heatmap
mat1 <- do.call(rbind, Zs1)
colnames(mat1) <- gsub(".Z1", "", names(z1$Z1))
rownames(mat1) <- colnames(Tc[,c(1,3,5)])
pheatmap(mat1, cluster_rows = FALSE, cluster_cols = TRUE, fontsize=7)

mat2 <- do.call(rbind, Zs2)
colnames(mat2) <- gsub(".Z1", "", names(z1$Z1))
rownames(mat2) <- colnames(Tc[,c(2,4,6)])
pheatmap(mat2, cluster_rows = FALSE, cluster_cols = TRUE, fontsize=7)



# Export Figures
svglite::svglite(filename = "../../data/export/human_phosphosite/kinase_enrichment_endurance.svg", width = 4.2, height = 4)
perturbPlot2d_edit(Tc=Tc[,c(1,2)], annotation=PhosphoSite.human, xlim=c(-5, 6), ylim=c(-4,4), main="Endurance")
dev.off()

svglite::svglite(filename = "../../data/export/human_phosphosite/kinase_enrichment_sprint.svg", width = 4.2, height = 4)
perturbPlot2d_edit(Tc=Tc[,c(3,4)], annotation=PhosphoSite.human, xlim=c(-5, 6), ylim=c(-4,4), main="Sprint")
dev.off()

svglite::svglite(filename = "../../data/export/human_phosphosite/kinase_enrichment_strength.svg", width = 4.2, height = 4)
perturbPlot2d_edit(Tc=Tc[,c(5,6)], annotation=PhosphoSite.human, xlim=c(-5, 6), ylim=c(-4,4), main="Resistance")
dev.off()

svglite::svglite(filename = "../../data/export/human_phosphosite/kinase_heatmap_post_pre.svg", width = 6, height = 2.5)
pheatmap(mat1, cluster_rows = FALSE, cluster_cols = TRUE, fontsize=7)
dev.off()


# Export heatmap data
write.table(mat1, file="../../data/export/human_phosphosite/heatmap_post_pre_v2.txt", sep="\t", quote=F)
write.table(mat2, file="../../data/export/human_phosphosite/heatmap_reco_pre_v2.txt", sep="\t", quote=F)
