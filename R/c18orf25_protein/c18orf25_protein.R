library(tidyverse)
library(PhosR)
library(limma)
library(msigdb)
library(GSEABase)
library(clusterProfiler)


prot <- read.delim("../../data/input/c18orf25_protein.txt", head=TRUE)
prot <- prot[!grepl("REV_|CON_", prot$Protein.IDs), ]


Ids <- prot %>% select(Protein.IDs, Gene.names)

prot <- prot %>% select(contains("Tibialis"))
prot[prot == 0] <- NA

Ids <- Ids[apply(is.na(prot), 1, sum) == 0, ]
prot <- prot[apply(is.na(prot), 1, sum) == 0, ]

prot <- log2(prot)


#### grouping
grps <- str_extract(colnames(prot), "_.._") %>% gsub("_", "", .)



# Median scale, two options
prot <- scale(prot, center = TRUE, scale = FALSE) %>% `rownames<-`(Ids$Protein.IDs)

plotQC(prot, grps = grps, panel="quantify", labels=colnames(prot))
plotQC(prot, grps = grps, panel="pca", labels=colnames(prot))



#######################
Treat <- factor(grps)
design <- model.matrix(~0+Treat)
colnames(design) <- levels(Treat)

corfit <- duplicateCorrelation(prot, design)
corfit$consensus

fit <- lmFit(prot, design, correlation=corfit$consensus)
cm <- makeContrasts(KO-WT, levels=design)

fit2 <- contrasts.fit(fit, cm)
fit2 <- eBayes(fit2)
tab <- topTable(fit2, number = Inf, adjust.method="BH")

c1 <- topTable(fit2, number = Inf, coef="KO - WT")
colnames(c1) <- paste("KO-WT", colnames(c1))



grand.tab <- c1 %>% mutate(Protein.IDs = rownames(.)) %>% left_join(Ids) %>% select(Protein.IDs, Gene.names, everything()) %>% rename("gene" = "Gene.names")

write.csv(grand.tab, file="../../data/export/c18orf25_protein/prot_grand_DE.csv", quote=F, row.names = TRUE)






comparisons <- "KO-WT logFC"
df <- grand.tab

## Create input data
d <- df %>% dplyr::select(gene, comparisons) %>% filter(!is.nan(comparisons) & !is.na(gene)) %>% dplyr::rename(gene = 1, fc = 2) %>% mutate(gene = sub("\\;.*", "", .$gene)) %>% group_by(gene) %>% summarize(fc = mean(fc))
  
geneList <- d[,2] %>% unlist
names(geneList) <- as.character(d[,1] %>% unlist)
geneList <- base::sort(geneList, decreasing = TRUE)
  
  
  
## -----------------------------------------------------------------------------

  
#load the MSigDB from the msigdb package
msigdb_mm = msigdb.v7.2.mm.SYM()
#append KEGG gene-sets
msigdb_mm = appendKEGG(msigdb_mm)
#dplyr::select h, c2, and c5 collections (recommended)
msigdb_mm = subsetCollection(msigdb_mm, c('h', 'c2', 'c5'))

# Create an input gene set for the GSEA function
genedb <- geneIds(msigdb_mm)
msig_db_gsea <- data.frame(gs_name = rep(names(genedb), lengths(genedb)), gene_symbol = unlist(genedb, use.names=TRUE))
  
# Run GSEA
set.seed(36)
em <- GSEA(geneList, TERM2GENE = msig_db_gsea)

write.csv(em@result, "../../data/export/c18orf25_protein/GSEA_result_KO_WT.csv", row.names = FALSE)