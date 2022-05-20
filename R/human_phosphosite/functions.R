# Modified function for PCA, from PhosR
# Modification includes an additional 'shape' parameter
pcaPlot_shapes <- function(mat, col, shape, labels){
  result <- pcaMethods::pca(t(mat), method = "ppca", nPcs = 2, seed = 123, main = "PCA")
  
  dat = data.frame(PC1 = pcaMethods::scores(result)[, 1], PC2 = pcaMethods::scores(result)[,2], col = col, shape = shape, Samples = labels)
  
  ggplot(dat, ggplot2::aes(x = .data$PC1, y = .data$PC2, color = .data$col, shape =  .data$shape)) + geom_point(size = 3, alpha = 0.8) + 
    labs(x = paste("PC1", round(result@R2[1] * 100), "%"), y = paste("PC2", round(result@R2[2] * 100), "%")) + theme_bw()
}

