library(directPA)
library(calibrate)
library(pheatmap)

## The input is the processed phospho data
Tc <- as.matrix(grand.tab[, c(1:3)])
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

z1 <- perturbPlot2d_edit(Tc=Tc[,c(1,2)], annotation=PhosR:::PhosphoSite.mouse, xlim=c(-5, 6), ylim=c(-4,4), main="Kinase perturbation analysis")


# Export Figure
svglite::svglite(filename = "../../data/export/c18orf25_phosphosite/kinase_enrichment_WTstim_con_vs_KOstim_con.svg", width = 4.2, height = 4)
perturbPlot2d_edit(Tc=Tc[,c(1,2)], annotation=PhosR:::PhosphoSite.mouse, xlim=c(-5, 6), ylim=c(-4,4), main="Kinase perturbation analysis")
dev.off()
