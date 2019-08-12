
# BiocManager::install("ChemmineR")


library(ChemmineR)

setwd(dirname(rstudioapi::getSourceEditorContext()$path))

hmdb_sdf <- read.SDFset("structures.sdf") 
valid <- validSDF(hmdb_sdf) # Identifies invalid SDFs in SDFset objects 
hmdb_sdf <- hmdb_sdf[valid] # Removes invalid SDFs, if there are any 

hmdb_summary <- data.frame(
   HMDB_ID = datablocktag(hmdb_sdf, tag="HMDB_ID"),
   Name = datablocktag(hmdb_sdf, tag="GENERIC_NAME"),
   Charge = datablocktag(hmdb_sdf, tag="JCHEM_FORMAL_CHARGE"),
   MF = datablocktag(hmdb_sdf, tag="FORMULA"),
   Exact_Mass = datablocktag(hmdb_sdf, tag="EXACT_MASS"),
   ALOGPS_LOGP = datablocktag(hmdb_sdf, tag="ALOGPS_LOGP"),
   JCHEM_LOGP = datablocktag(hmdb_sdf, tag="JCHEM_LOGP"),
   JCHEM_NUMBER_OF_RINGS = datablocktag(hmdb_sdf, tag="JCHEM_NUMBER_OF_RINGS"),
   SMILES = datablocktag(hmdb_sdf, tag="SMILES"),
   JCHEM_PHYSIOLOGICAL_CHARGE = datablocktag(hmdb_sdf, tag="JCHEM_PHYSIOLOGICAL_CHARGE"),
   JCHEM_PKA = datablocktag(hmdb_sdf, tag="JCHEM_PKA"),
   JCHEM_PKA_STRONGEST_ACIDIC = datablocktag(hmdb_sdf, tag="JCHEM_PKA_STRONGEST_ACIDIC"),
   JCHEM_PKA_STRONGEST_BASIC = datablocktag(hmdb_sdf, tag="JCHEM_PKA_STRONGEST_BASIC"),
   JCHEM_POLAR_SURFACE_AREA = datablocktag(hmdb_sdf, tag="JCHEM_POLAR_SURFACE_AREA"),
   JCHEM_RULE_OF_FIVE = datablocktag(hmdb_sdf, tag="JCHEM_RULE_OF_FIVE"),
   JCHEM_ROTATABLE_BOND_COUNT = datablocktag(hmdb_sdf, tag="JCHEM_ROTATABLE_BOND_COUNT"),
   JCHEM_GHOSE_FILTER = datablocktag(hmdb_sdf, tag="JCHEM_GHOSE_FILTER"),
   JCHEM_MDDR_LIKE_RULE = datablocktag(hmdb_sdf, tag="JCHEM_MDDR_LIKE_RULE"),
   JCHEM_VEBER_RULE = datablocktag(hmdb_sdf, tag="JCHEM_VEBER_RULE"),
   JCHEM_REFRACTIVITY = datablocktag(hmdb_sdf, tag="JCHEM_REFRACTIVITY")
   )

test = datablocktag(hmdb_sdf)
test2 = cbind(test)

test = as.data.frame(datablocktag(hmdb_sdf[[1]]))

plot(hmdb_sdf[1], print=FALSE)

saveRDS(hmdb_summary, "HMDB_structure_sdf.rds")
write.table(hmdb_summary, file = "HMDB_structure_sdf.tsv", sep = "\t",
            col.names = T, row.names = F, quote = F)

# 
# 
# 
# 
# setwd("C:/Users/lc8/Desktop/hmdb_metabolites")
# 
# library(XML)
# input <- "hmdb_metabolites.xml"
# items <- NULL
# maxItems <- 50
# 
# parseItem = function (parser, node, ...) {
#   children <- xmlChildren(node)
#   items <<- rbind(items, sapply(children, xmlValue))
#   if (nrow(items) == maxItems) {
#     xmlStopParser(parser)
#   }
# }
# 
# # with XMLParserContextFunction, we get the parser as first parameter
# # so we can call xmlStopParser
# class(parseItem) = c("XMLParserContextFunction", "SAXBranchFunction")
# 
# xmlEventParse(input,
#               branches = list(item = parseItem),
#               ignoreBlanks = T
# )
# 
# items <- as.data.frame(items)



