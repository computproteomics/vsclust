## code to prepare `DATASET` dataset goes here
artificial_clusters <- read.csv("data-raw/ArtData.csv")
protein_expressions <- read.csv("data-raw/ProtExample.csv",row.names=1)
protein_expressions <- protein_expressions[2:nrow(protein_expressions),]
usethis::use_data(protein_expressions, overwrite = TRUE)
usethis::use_data(artificial_clusters, overwrite = TRUE)
