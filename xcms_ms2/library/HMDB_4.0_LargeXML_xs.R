require(XML)
require(rlist)

result <- xmlParse(file = "C:/Users/Li Chen/Desktop/sweat_metabolites_origin.xml")
rootnode <- xmlRoot(result)
MetaboliteNumber <- xmlSize(rootnode)

if.not.null <- function(x) if(!is.null(x)) x else "NA"

print(paste("This version of HMDB contains",MetaboliteNumber,"metabolites."))

#MetaboliteNumber <- 100

Sys.time()

ResultsTable <- data.frame(No.=c(1:MetaboliteNumber)) 
for (i in 1:MetaboliteNumber) {
  NodeXML <- xmlToList(rootnode[[i]])
  ResultsTable[i,2] <- NodeXML$accession
  ResultsTable[i,3] <- NodeXML$name
  ResultsTable[i,4] <- NodeXML$chemical_formula
  ResultsTable[i,5] <- if.not.null(NodeXML$monisotopic_molecular_weight)
  if(is.list(NodeXML$ontology)) {
    Endo <- unlist(list.search(NodeXML$ontology,identical(., 'Endogenous')))
    if(!is.null(Endo)) {ResultsTable[i,6] <-  Endo}
  }
  if(is.list(NodeXML$taxonomy)) {
  ResultsTable[i,7] <- if.not.null(NodeXML$taxonomy$super_class)
  ResultsTable[i,8] <- if.not.null(NodeXML$taxonomy$class)
  }
  ResultsTable[i,9] <- NodeXML$status
  ResultsTable[i,10] <- NodeXML$inchi
  if(i%%100==0) print(paste(i,"metabolites are done. Time is",Sys.time()))
}

Sys.time()

ResultsTable_ls = list()
for (i in 1:MetaboliteNumber) {
  NodeXML <- xmlToList(rootnode[[i]])
  accession <- NodeXML$accession
  name <- NodeXML$name
  chemical_formula <- NodeXML$chemical_formula
  monisotopic_molecular_weight <- if.not.null(NodeXML$monisotopic_molecular_weight)
  if(is.list(NodeXML$ontology)) {
    Endo <- unlist(list.search(NodeXML$ontology,identical(., 'Endogenous')))
    if(!is.null(Endo)) {endo <-  Endo}
  }
  if(is.list(NodeXML$taxonomy)) {
    super_class <- if.not.null(NodeXML$taxonomy$super_class)
    class <- if.not.null(NodeXML$taxonomy$class)
  }
  status <- NodeXML$status
  inchi <- NodeXML$inchi
  
  temp_list = list(No. = i,
                   accession = accession,
                   name = name,
                   chemical_formula = chemical_formula,
                   monisotopic_molecular_weight = monisotopic_molecular_weight,
                   endo = endo,
                   super_class = super_class,
                   class = class,
                   status = status,
                   inchi = inchi
                   )
  ResultsTable_ls[[i]]=temp_list
  
  if(i%%100==0) print(paste(i,"metabolites are done. Time is",Sys.time()))
}

ResultsTable = dplyr::bind_rows(ResultsTable_ls)

Sys.time()

write.csv(ResultsTable,file="HMDB_Extracted_0801.csv")