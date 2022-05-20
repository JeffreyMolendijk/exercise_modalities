# Human protein analysis

# Description:

# Data generated:
# Figure XX, XX, XX, XX
# Table xx, xx


##############
### SETUP  ###
##############

# Load libraries and functions
library(tidyverse)
library(PhosR)
library(limma)
library(funscoR)
library(WGCNA)
library(ggpubr)
library(forecast)
library(sva)

source("./functions.R")

dir.create(path = "../../data/export/human_protein", recursive = TRUE, mode = 744)


#################
### LOAD DATA ###
#################

prot.raw <- read.delim(file = "../../data/input/human_protein.txt", na.strings = c("NaN", "NA")) %>% as.data.frame()
prot.raw <- prot.raw[!grepl("REV|CON_", prot.raw$Protein.IDs), ]

# Fix minor error in column names, where Sprnt should be Sprint
colnames(prot.raw) <- sub("Sprnt", "Sprint", colnames(prot.raw)) %>% gsub(" ", ".", .)

# Generate an ID column
Ids <- prot.raw$Protein.IDs

# Assign Ids to rownames of prot.raw
rownames(prot.raw) <- Ids



########################################
### TMT INTERNAL STANDARD CORRECTION ###
########################################

# Subtract internal standard values from other samples
prot.exp <- list()

for(i in 1:8){
  
  df <- prot.raw %>% as.data.frame() %>% select(contains(paste("Subject", i, "_", sep = "")))
  df[df == 0] <- NA
  df <- log2(df)
  
  prot.exp[[i]] <- (df %>% select(-contains("internal.standard"))) - ((df %>% select(contains("internal.standard"))) %>% as.matrix() %>% .[,1])
}

# Combine internal standard subtracted data back into a matrix
prot.exp <- do.call(cbind, prot.exp)



############################
### PHOSR PRE-PROCESSING ###
############################

# Define experimental groups
grps <- gsub(".+_", "", colnames(prot.exp))

# Perform filtering to retain sites measured in >3/8 Subjects in at least 1 treatment group
prot.filtered <- selectGrps(prot.exp, grps, percent = (3/8), n=1)

# Median scaling
set.seed(1)
prot.scale <- medianScaling(prot.filtered)


######################
### LIMMA ANALYSIS ###
######################

# Define the model matrix for lmFit
Treat <- factor(grps)
design <- model.matrix(~0+Treat)
colnames(design) <- levels(Treat)

# Estimate intra-block correlation where block refers to subjects
corfit3 <- duplicateCorrelation(prot.scale, design, block=gsub("_.+", "", colnames(prot.scale)))

lf <- lmFit(prot.scale, design, block=gsub("_.+", "", colnames(prot.scale)), correlation=corfit3$consensus)
cm <- makeContrasts(PostEndurance-PreEndurance, RecoveryEndurance-PreEndurance,
                    PostSprint-PreSprint, RecoverySprint-PreSprint,
                    PostStrength-PreStrength, RecoveryStrength-PreStrength, levels=design)
fit <- contrasts.fit(lf, cm)
fit <- eBayes(fit)
tab <- topTableF(fit, number = Inf, adjust.method="BH")

# Export the q-value distribution
svglite::svglite(filename = "../../data/export/human_protein/padjhist.svg", width = 5, height = 4)
hist(tab$adj.P.Val, main="adj p-value distribution")
dev.off()

# Combine DE results into a table (grand.tab)
c1 <- topTable(fit, number = Inf, coef="PostEndurance - PreEndurance")
c2 <- topTable(fit, number = Inf, coef="RecoveryEndurance - PreEndurance")
c3 <- topTable(fit, number = Inf, coef="PostSprint - PreSprint")
c4 <- topTable(fit, number = Inf, coef="RecoverySprint - PreSprint")
c5 <- topTable(fit, number = Inf, coef="PostStrength - PreStrength")
c6 <- topTable(fit, number = Inf, coef="RecoveryStrength - PreStrength")

o <- rownames(tab)
c1.tab <- c1[o, c(3:5)]
c2.tab <- c2[o, c(3:5)]
c3.tab <- c3[o, c(3:5)]
c4.tab <- c4[o, c(3:5)]
c5.tab <- c5[o, c(3:5)]
c6.tab <- c6[o, c(3:5)]
colnames(c1.tab) <- paste("PostEndurance-PreEndurance", colnames(c1[o, c(3:5)]))
colnames(c2.tab) <- paste("RecoveryEndurance-PreEndurance", colnames(c2[o, c(3:5)]))
colnames(c3.tab) <- paste("PostSprint-PreSprint", colnames(c3[o, c(3:5)]))
colnames(c4.tab) <- paste("RecoverySprint-PreSprint", colnames(c4[o, c(3:5)]))
colnames(c5.tab) <- paste("PostStrength-PreStrength", colnames(c5[o, c(3:5)]))
colnames(c6.tab) <- paste("RecoveryStrength-PreStrength", colnames(c6[o, c(3:5)]))

grand.tab <- cbind(tab, c1.tab, c2.tab, c3.tab, c4.tab, c5.tab, c6.tab)
write.table(grand.tab, file="../../data/export/human_protein/prot_grand_DE.txt", quote=F, sep="\t")



#################################
### Perform removeBatcheffect ###
#################################

# Perform removeBatcheffect on data
prot.rba <- prot.scale %>% removeBatchEffect(batch = gsub("_.+", "", colnames(prot.scale)), design = design, correlation=corfit3$consensus)

# Export data and PCA plots (before and after processing)
write.csv(prot.rba, file = "../../data/export/human_protein/human_protein_rba.csv")

svglite::svglite(filename = "../../data/export/human_protein/pca_before.svg", width = 5.5, height = 4)
pcaPlot_shapes(prot.raw %>% select(contains("Subject")) %>% select(-contains("internal")) %>% tImpute(), col=sub("_.*", "", colnames(prot.rba)), shape=rep(c("Pre", "Post", "Recovery"), 24), labels=gsub("Subject\\d_", "", colnames(prot.scale)))
dev.off()

svglite::svglite(filename = "../../data/export/human_protein/pca_after.svg", width = 5.5, height = 4)
pcaPlot_shapes(prot.rba, col=rep(c(rep(c("Endurance"), 3), rep(c("Sprint"), 3), rep(c("Resistance"), 3)), 8), shape=rep(c("Pre", "Post", "Recovery"), 24), labels=gsub("Subject\\d_", "", colnames(prot.scale)))
dev.off()

# Remove c* objects
rm(list = ls(pattern = '^c[0-9].*'))
rm(lf, cm, fit, corfit3, tab, design, o)

