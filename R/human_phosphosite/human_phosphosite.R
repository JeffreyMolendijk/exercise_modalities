# Human phosphosite analysis


##############
### SETUP  ###
##############

# Load libraries and functions
library(tidyverse)
library(PhosR)
library(limma)
library(WGCNA)
library(ggpubr)
library(forecast)
library(sva)

source("./functions.R")

# Load PSP data required for kinase enrichment analysis
data(PhosphoSitePlus)

# Create folders to save exported plots
dir.create(path = "../../data/export/human_phosphosite", recursive = TRUE, mode = 744)
dir.create(path = "../../data/export/human_phosphosite/scatterplot", recursive = TRUE, mode = 744)
dir.create(path = "../../data/export/human_phosphosite/spearman", recursive = TRUE, mode = 744)


#################
### LOAD DATA ###
#################

phospho.raw <- readxl::read_excel(path = "../../data/input/human_phospho.xlsx", na = c("NaN", "NA")) %>% as.data.frame()

# Fix minor error in column names, where Sprnt should be Sprint and replace spaces with dots
colnames(phospho.raw) <- sub("Sprnt", "Sprint", colnames(phospho.raw)) %>% gsub(" ", ".", .)

# Generate an ID column
Ids <- paste(paste(sapply(strsplit(as.character(phospho.raw[,"Protein"]), ";"), function(x){x[1]}),
                   sapply(strsplit(as.character(phospho.raw[,"Gene.names"]), ";"), function(x){x[1]}),
                   paste(as.character(phospho.raw[,"Amino.acid"]), as.character(phospho.raw[,"Position"]), sep=""),
                   sapply(strsplit(as.character(phospho.raw[,"Sequence.window"]), ";"), function(x){x[1]}), sep="~"),
             gsub("___", "", as.character(phospho.raw[,"Multiplicity"])), as.character(phospho.raw[,"Unique.identifier"]), sep="~")

# Assign Ids to rownames of phospho.raw
rownames(phospho.raw) <- Ids



########################################
### TMT INTERNAL STANDARD CORRECTION ###
########################################

# Subtract internal standard values from other samples
phospho.exp <- list()

for(i in 1:8){
  
  phospho.exp[[i]] <- phospho.raw %>% as.data.frame() %>% select(contains(paste("Subject", i, "_", sep = "")))
  phospho.exp[[i]] <- (phospho.exp[[i]] %>% select(-contains("internal.standard"))) - ((phospho.exp[[i]] %>% select(contains("internal.standard"))) %>% as.matrix() %>% .[,1])

  }

# Combine internal standard subtracted data back into a matrix
phospho.exp <- do.call(cbind, phospho.exp)


############################
### PHOSR PRE-PROCESSING ###
############################

# Define experimental groups
grps <- gsub(".+_", "", colnames(phospho.exp))

# Collapse phosphosites
id <- sapply(strsplit(rownames(phospho.exp), "~"), function(x){paste(toupper(x[2]), x[3], x[4], sep=";")})
phospho.site <- phosCollapse(phospho.exp, id = id, stat = apply(abs(phospho.exp), 1, max, na.rm=TRUE), by = "max")

# Perform filtering to retain sites measured in >3/8 Subjects in at least 1 treatment group
phospho.filtered <- selectGrps(phospho.site, grps, percent = (3/8), n=1)

# Median scaling
set.seed(1)
phospho.scale <- medianScaling(phospho.filtered)


######################
### LIMMA ANALYSIS ###
######################

# Define the model matrix for lmFit
Treat <- factor(grps)
design <- model.matrix(~0+Treat)
colnames(design) <- levels(Treat)

# Estimate intra-block correlation where block refers to subjects
corfit3 <- duplicateCorrelation(phospho.scale, design, block=gsub("_.+", "", colnames(phospho.scale)))

lf <- lmFit(phospho.scale, design, block=gsub("_.+", "", colnames(phospho.scale)), correlation=corfit3$consensus)
cm <- makeContrasts(PostEndurance-PreEndurance, RecoveryEndurance-PreEndurance,
                    PostSprint-PreSprint, RecoverySprint-PreSprint,
                    PostStrength-PreStrength, RecoveryStrength-PreStrength, levels=design)
fit <- contrasts.fit(lf, cm)
fit <- eBayes(fit)
tab <- topTableF(fit, number = Inf, adjust.method="BH")

# Export the q-value distribution
svglite::svglite(filename = "../../data/export/human_phosphosite/padjhist.svg", width = 5, height = 4)
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
write.table(grand.tab, file="../../data/export/human_phosphosite/phospho_grand_DE.txt", quote=F, sep="\t")



#################################
### Perform removeBatcheffect ###
#################################

# Perform removeBatcheffect on data
phospho.rba <- phospho.scale %>% removeBatchEffect(batch = gsub("_.+", "", colnames(phospho.scale)), design = design, correlation=corfit3$consensus)

# Export data and PCA plots (before and after processing)
write.csv(phospho.rba, file = "../../data/export/human_phosphosite/human_phosphosite_rba.csv")

svglite::svglite(filename = "../../data/export/human_phosphosite/pca_before.svg", width = 5.5, height = 4)
pcaPlot_shapes(phospho.raw %>% select(contains("Subject")) %>% select(-contains("internal")), col=sub("_.*", "", colnames(phospho.rba)), shape=rep(c("Pre", "Post", "Recovery"), 24), labels=gsub("Subject\\d_", "", colnames(phospho.scale)))
dev.off()

svglite::svglite(filename = "../../data/export/human_phosphosite/pca_after.svg", width = 5.5, height = 4)
pcaPlot_shapes(phospho.rba, col=rep(c(rep(c("Endurance"), 3), rep(c("Sprint"), 3), rep(c("Resistance"), 3)), 8), shape=rep(c("Pre", "Post", "Recovery"), 24), labels=gsub("Subject\\d_", "", colnames(phospho.scale)))
dev.off()




##############################
### Run FunScor annotation ###
##############################

# Execute workflow in funscor.R
source("./funscor.R")


##############################
### Site-trait correlation ###
##############################

# Execute workflow in phosphosite_trait_correlation.R
source("./phosphosite_trait_correlation.R")


##############################
### RUN KINASE ENRICHMENT  ###
##############################

# Execute workflow in kinase_enrichment.R
source("./kinase_enrichment.R")
