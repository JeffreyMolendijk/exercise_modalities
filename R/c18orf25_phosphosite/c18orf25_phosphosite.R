# C18ORF25 phosphosite analysis


##############
### SETUP  ###
##############

# Load libraries and functions
library(PhosR)
library(limma)


#################
### LOAD DATA ###
#################

phospho.processed <- read.delim("../../data/input/c18orf25_phospho.txt", head=TRUE)

# Concatenate columns into unique phosphosite identifiers
Ids <- paste(paste(sapply(strsplit(as.character(phospho.processed[,"Protein"]), ";"), function(x){x[1]}),
             sapply(strsplit(as.character(phospho.processed[,"Gene.names"]), ";"), function(x){x[1]}),
             paste(as.character(phospho.processed[,"Amino.acid"]), as.character(phospho.processed[,"Position"]), sep=""),
             sapply(strsplit(as.character(phospho.processed[,"Sequence.window"]), ";"), function(x){x[1]}), sep="~"),
             gsub("___", "", as.character(phospho.processed[,"Multiplicity"])), as.character(phospho.processed[,"Unique.identifier"]), sep="~")
rownames(phospho.processed) <- Ids

phospho.seqs <- sapply(strsplit(as.character(phospho.processed[,"Sequence.window"]), ";"), function(x){x[1]})
names(phospho.seqs) <- Ids

# Replace NaN with NA
phospho.exp <- as.matrix(phospho.processed[,1:18])
phospho.exp[is.nan(phospho.exp)] <- NA


#################
### filtering ###
#################

# Define experimental groupings used for filtering and imputation
grps <- gsub("_\\d", "", colnames(phospho.exp))

id <- sapply(strsplit(rownames(phospho.exp), "~"), function(x){paste(toupper(x[2]), x[3], x[4], sep=";")})
phospho.site <- phosCollapse(phospho.exp, id = id, stat = apply(abs(phospho.exp), 1, max, na.rm=TRUE), by = "max")

# filtering by group quantification
phospho.filtered <- selectGrps(phospho.site, grps, percent = 0.75, n=1)
dim(phospho.filtered)


##################
### imputation ###
##################

# Perform scImpute and ptImpute
phospho.impute <- medianScaling(scImpute(mat=phospho.filtered, percent = 0.5, grps = grps), grps=grps, reorder = TRUE)
phospho.impute[,1:5] <- ptImpute(phospho.impute[,6:10], phospho.impute[,1:5], percent1 = 0.5, percent2 = 0, paired = FALSE)
phospho.impute[,11:14] <- ptImpute(phospho.impute[,15:18], phospho.impute[,11:14], percent1 = 0.5, percent2 = 0, paired = FALSE)

grps <- gsub("_\\d", "", colnames(phospho.impute))



######################
### LIMMA ANALYSIS ###
######################

# Define the model matrix for lmFit
Treat <- factor(grps)
design <- model.matrix(~0+Treat)
colnames(design) <- levels(Treat)

# Estimate intra-block correlation
corfit <- duplicateCorrelation(phospho.impute, design, block=sapply(strsplit(colnames(phospho.impute), "_"), function(x)x[3]))

fit <- lmFit(phospho.impute, design, block=sapply(strsplit(colnames(phospho.impute), "_"), function(x)x[3]), correlation=corfit$consensus)
cm <- makeContrasts(WT=WT_stim-WT_ctrl, KO=KO_stim-KO_ctrl, Diff=(KO_stim-KO_ctrl)-(WT_stim-WT_ctrl), levels=design)

fit2 <- contrasts.fit(fit, cm)
fit2 <- eBayes(fit2)
tab <- topTable(fit2, number = Inf, adjust.method="BH")
tab[1:10,]

c1 <- topTable(fit2, number = Inf, coef="WT")
c2 <- topTable(fit2, number = Inf, coef="KO")
c3 <- topTable(fit2, number = Inf, coef="Diff")

# Combine DE results into a table (grand.tab)
o <- rownames(tab)
c1.tab <- c1[o, c(3:5)]
c2.tab <- c2[o, c(3:5)]
c3.tab <- c3[o, c(3:5)]

colnames(c1.tab) <- paste("WT=WT_stim-WT_ctrl", colnames(c1[o, c(3:5)]))
colnames(c2.tab) <- paste("KO=KO_stim-KO_ctrl", colnames(c2[o, c(3:5)]))
colnames(c3.tab) <- paste("Diff=(KO_stim-KO_ctrl)-(WT_stim-WT_ctrl)", colnames(c3[o, c(3:5)]))

sum(c1.tab[,3] < 0.05, na.rm = TRUE) / nrow(c1.tab)
sum(c2.tab[,3] < 0.05, na.rm = TRUE) / nrow(c2.tab)

kinases <- rep(names(PhosR:::PhosphoSite.mouse[1:250]), time=sapply(PhosR:::PhosphoSite.mouse[1:250], length))
names(kinases) <- unlist(PhosR:::PhosphoSite.mouse[1:250])
sapply(strsplit(o, ";"), function(x)paste(x[1], x[2], "", sep=";"))
sites <- sapply(strsplit(o, ";"), function(x)paste(x[1], x[2], "", sep=";"))
idx <- which(sites %in% names(kinases))
kinaseAnnotation <- rep(NA, length(o))
kinaseAnnotation[idx] <- kinases[sites[idx]]

grand.tab <- cbind(tab, c1.tab, c2.tab, c3.tab, kinaseAnnotation)

write.csv(grand.tab, file="../../data/export/c18orf25_phosphosite/phospho_grand_DE.csv", quote=F, row.names = TRUE)



##############################
### RUN KINASE ENRICHMENT  ###
##############################

# Execute workflow in kinase_enrichment.R
source("./kinase_enrichment.R")
