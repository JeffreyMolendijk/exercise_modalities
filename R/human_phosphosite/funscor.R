###################################
### Load packages and functions ###
###################################

library(funscoR)


#############################
### Create input data     ###
#############################

# Create an annotation table to map between the 'Ids' and 'id' arrays
phospho.annotation <- phospho.raw %>% select(-id) %>% 
  mutate(Ids = Ids, id = id) %>% select(-contains("Subject")) %>% 
  select(Ids, id, everything()) %>% filter(id %in% rownames(grand.tab))

# To run funScor, generate a table containing protein accession, residue position and residue amino acid
funscor <- phospho.annotation %>% select(Protein, Position, Amino.acid) %>% 
  `colnames<-`(c("acc", "position", "residue")) %>% `rownames<-`(c(1:nrow(.))) %>% 
  filter(!is.na(acc) & !is.na(position)) %>% distinct()

## annotate phosphoproteome with features
funscor <- annotate_sites(funscor)


##########################################
### Preprocess data and generate model ###
##########################################

## preprocess features for training
ST_features <- preprocess_features(funscor, "ST")
Y_features <- preprocess_features(funscor, "Y")

## train new model
ST_model <- train_funscore(ST_features, "ST", psp, ncores = 4)
Y_model <- train_funscore(Y_features, "Y", psp, ncores = 4)

## predict funcscoR for all sites
ST_scores <- predict_funscore(ST_features, ST_model, ncores = 4)
Y_scores <- predict_funscore(Y_features, Y_model, ncores = 4)


############################################
### Generate predictions and annotations ###
############################################

## gather all predictions
all_scores <- bind_rows(ST_scores, Y_scores) %>% mutate(probabilities = log_scaling(probabilities))

# Join predictions onto the annotated sites table
funscor <- left_join(funscor %>% mutate(sites = paste(acc, position, sep = "_")), all_scores)

# Join all annotations with the original sites
phospho.annotation <- phospho.annotation %>% mutate(sites = paste(Protein, Position, sep = "_")) %>% left_join(., funscor, by = "sites")


############################################
### Export results                       ###
############################################

# Export all annotations
write.csv(phospho.annotation, file = "../../data/export/human_phosphosite/site_annotations.csv")
