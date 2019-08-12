library(XML)
library(xml2)
# devtools::install_github('dantonnoriega/xmltools')
library(xmltools)
library(dplyr)

setwd("C:/Users/lc8/Downloads/sweat_metabolites")
xmlFileName = "sweat_metabolites.xml"

hmdb <- read_xml(xmlFileName) %>%
  xml2::as_list()

getremove <- function(list){
  list$synthesis_reference <- NULL
  list$general_references <- NULL
  list$version <- NULL
  list$creation_date <- NULL
  list$update_date <- NULL
  list$secondary_accessions <- NULL
  list$description <- NULL
  list$synonyms <- NULL
  list$traditional_iupac <- NULL
  list$spectra <- NULL
  list$predicted_properties <- NULL
  list$taxonomy$alternative_parents <- NULL
  list$taxonomy$substituents <- NULL
  list$taxonomy$description<- NULL
  list$normal_concentrations <- NULL
  list$abnormal_concentrations <- NULL
  list$diseases <- NULL
  list$ontology <- NULL
  list$direct_parent <- list$taxonomy$direct_parent
  list$kingdom <- list$taxonomy$kingdom
  list$super_class <- list$taxonomy$super_class
  list$class <- list$taxonomy$class
  list$sub_class <- list$taxonomy$sub_class
  list$molecular_framework <- list$taxomomy$molecular_framework
  list$taxonomy <- NULL
  list$protein_associations <- NULL
  list$biofluid_locations <- NULL
  list$pathways <- NULL
  list$tissue_locations <- NULL
  list$experimental_properties <- NULL
  return(list)
}
subhmdb <- lapply(hmdb,getremove)
id = unlist(sapply(subhmdb, "[[", "accession"))
name = unlist(sapply(subhmdb, "[[", "name"))
hmdbclass = sapply(subhmdb, "[[", "class")
hmdbclass2 = sapply(hmdbclass, unlist)
class <- as.character(hmdbclass2)
hmdbsubclass = sapply(subhmdb, "[[", "sub_class")
hmdbsubclass2 = sapply(hmdbsubclass, unlist)
subclass <- as.character(hmdbsubclass2)
hmdbMW = sapply(subhmdb, "[[", "monisotopic_molecular_weight")
hmdbMW2 = sapply(hmdbMW, unlist)
MW <- as.numeric(as.character(hmdbMW2))
table <- cbind.data.frame(id = id,name = name,class = class,subclass = subclass, MW = MW)
write.csv(table,file = 'hmdb.csv')