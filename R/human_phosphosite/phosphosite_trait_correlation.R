##############################
### Site-trait correlation ###
##############################

# Read traits data

# List the sheets in physiological data excel
sheets <- readxl::excel_sheets("../../data/input/human_trait.xlsx")

columns <- c(
  paste("Subject", c(1:8), "_Pre", c(rep("Endurance", 8), rep("Sprint", 8), rep("Strength", 8)),  sep = ""),
  paste("Subject", c(1:8), "_Post", c(rep("Endurance", 8), rep("Sprint", 8), rep("Strength", 8)),  sep = ""),
  paste("Subject", c(1:8), "_Recovery", c(rep("Endurance", 8), rep("Sprint", 8), rep("Strength", 8)),  sep = "")
)

traits <- list()

for(i in 1:length(sheets)){
  
  x <- readxl::read_excel("../../data/input/human_trait.xlsx", sheet = sheets[i])
  
  traits[[sheets[i]]] <- c(
    x[2,-1] %>% unlist,
    x[3,-1] %>% unlist,
    x[4,-1] %>% unlist
  ) %>% data.frame() %>% t %>% `colnames<-`(columns)
  
}

traits <- do.call(rbind, traits) %>% as.data.frame() %>% `rownames<-`(sheets)


# Check whether each trait is normal, otherwise transform
shapiro <- list()

for(i in 1:nrow(traits)){
  
  shapiro[[i]] <-  c(rownames(traits)[i], 
                     shapiro.test(traits[i,] %>% unlist)$p.value %>% formatC(digits = 2),
                     shapiro.test(log10(traits[i,] %>% unlist))$p.value %>% formatC(digits = 2),
                     shapiro.test((traits[i,] %>% unlist)^2)$p.value %>% formatC(digits = 2),
                     shapiro.test((traits[i,] %>% unlist)^0.5)$p.value %>% formatC(digits = 2),
                     shapiro.test(BoxCox(traits[i,] %>% unlist, BoxCox.lambda(traits[i,] %>% unlist)))$p.value %>% formatC(digits = 2)
  ) 
  
}

do.call(rbind, shapiro) %>% `colnames<-`(c("trait", "raw", "log", "square", "root", "boxcox"))


x <- WGCNA::corAndPvalue(x = phospho.rba %>% as.data.frame() %>% dplyr::select(colnames(traits)) %>% t %>% `colnames<-`(rownames(phospho.rba)), traits %>% t, method = "spearman")

x <- left_join(x$cor %>% reshape2::melt() %>% `colnames<-`(c("Var1", "Var2", "cor")), x$p %>% reshape2::melt() %>% `colnames<-`(c("Var1", "Var2", "p"))) %>% mutate(q = p.adjust(p = .$p, method = "BH"))

write.csv(x, file = "../../data/export/human_phosphosite/phosphosite_trait_correlation_spearman.csv", row.names = FALSE)



traitnames <- x$Var2 %>% unique

for(i in 1:length(traitnames)){
  
  range <- x %>% filter(Var2 == traitnames[i]) %>% .$cor %>% max(abs(.))
  
  p <-  ggplot(x %>% filter(Var2 == traitnames[i]) %>% arrange(cor) %>% mutate(index = 1:nrow(.)), aes(x = index, y = cor, fill = q < 0.05)) + geom_bar(stat = "identity") + labs(subtitle = traitnames[i]) + theme_bw() + annotate("text", x = 1, y = Inf, hjust = 0, vjust = 1, label = paste("# significant: ", x %>% filter(Var2 == traitnames[i]) %>% filter(q < 0.05) %>% nrow(), "\n(q < 0.05)", sep = "")) + lims(y = c(-range, range)) + theme(legend.position = "null")
  
  pdf(file = paste("../../data/export/human_phosphosite/spearman/traitcor", traitnames[i],".pdf", sep = ""), width = 4, height = 3)
  plot(p)
  dev.off()
  
}

######################################
### Create individual scatterplots ###
######################################

df_ind <- phospho.rba %>% as.data.frame %>% select(intersect(colnames(phospho.rba), colnames(traits))) 
trait_ind <- traits %>% select(intersect(colnames(phospho.rba), colnames(traits)))

colnames(df_ind) == colnames(trait_ind)

df_ind$sites <- rownames(phospho.rba)
trait_ind$trait <- trait_ind %>% rownames()

ab <- c("C18ORF25;S67;ISSMPCLLMELRRDSSESQLASTESDKPTTG", 
        "C18ORF25;S67;ISSMPCLLMELRRDSSESQLASTESDKPTTG", 
        "C18ORF25;S67;ISSMPCLLMELRRDSSESQLASTESDKPTTG", 
        "LDHA;Y39;FGSKSNMATLKDQLIYNLLKEEQTPQNKITV", 
        "AKAP2;S393;GPPEDSGASAAKGQKSPGALETPSAAGSQGN", 
        "MAP2K4;S80;PTGVQNPHIERLRTHSIESSGKLKISPEQHW",
        "RPS6;S205;KEKRQEQIAKRRRLSSLRASTSKSESSQK__", 
        "PPP1R2;S87;STPYHSMMGDDEDACSDTEATEAMAPDILAR",
        "PAK1;S174;LNVKAVSETPAVPPVSEDEDDDDDDATPPPV", 
        "GYS1;S585;YRYPRPASVPPSPSLSRHSSPHQSEDEEDPR", 
        "DHCR7;S14;__MAAKSQPNIPKAKSLDGVTNDRTASQGQW", 
        "PDHA1;S269;QTYRYHGHSMSDPGVSYRTREEIQEVRSKSD")

ba <- c("Adrenalin", "Lactate", "Noradrenalin", "Lactate", "FFA","Glucose", "Glucose", "Insulin", "Insulin", "Glycogen", "FFA", "FFA")


for(i in 1:length(ab)){
  
  x <- cbind(df_ind %>% filter(sites == ab[i]) %>% select(-sites) %>% t, 
             trait_ind %>% filter(trait == ba[i]) %>% select(-trait) %>% t) %>% `colnames<-`(c("site", "trait")) %>% as.data.frame() %>% mutate(exercise = rownames(.) %>% sub(".*_", "", .) %>% sub("Pre|Post|Recovery", "", .)) %>% mutate(condition = rownames(.) %>% sub(".*_", "", .) %>% sub("Endurance|Sprint|Strength", "", .))
  
  pdf(file = paste("../../data/export/human_phosphosite/scatterplot/", ba[i], "_",janitor::make_clean_names(ab[i]),".pdf", sep = ""), width = 6, height = 3)
  print(ggplot(x, aes(x = site, y = log10(trait), col = exercise)) + geom_smooth(data = x, aes(group = 1), alpha=0.3, size=0, method = "lm") + stat_smooth(data = x, aes(group = 1), geom="line", alpha=0.3, size=1, method = "lm") + geom_point(aes(shape = condition)) + labs(x = ab[i], y = paste0('log10 ', ba[i], sep = "")) + theme_bw())
  dev.off()
  
}