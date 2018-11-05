#install.packages("enviPat")
#Testedit

library(enviPat)
library(rcdk)

data(adducts)

data(chemforms)
data(isotopes)
checked<-check_chemform(isotopes,chemforms)
checked

chemforms<-c("C900Cl4H49","O82394","C8G500Zn9","Br1","6DBr9889")
data(isotopes)
checked<-check_chemform(isotopes,chemforms)$new_formula
checked

test = "NNNN[15]NNNN[13]CCCHHDDDHH"
test_checked = check_chemform(isotopes, test)

formulas<-c("C8H4Cl2","C10H16O2","C3H10")
deduct<-c("C4H10")
# TRUE if deduct is not contained and FALSE otherwise.
check_ded(formulas, deduct)

data(isotopes)
data(chemforms)
pattern<-isopattern(
  isotopes,
  chemforms,
  threshold=0.1,
  plotit=TRUE,
  charge=FALSE,
  emass=0.00054858,
  algo=1
)
#It will report if any isotope peaks may be confusing with others.
test=check_several(pattern,dmz=0.001,ppm=FALSE)


test["mass"]=NA
for (i in 1:nrow(test)){
  test$mass[i]=get.formula(as.character(test$compound[i]))@mass
  
}

data(isotopes)
data(chemforms)
chemforms<-chemforms[1:5]
pattern<-isopattern(
  isotopes,
  chemforms,
  threshold=1,
  plotit=F,
  charge=FALSE,
  emass=0.00054858,
  algo=2
)
profiles<-envelope(
  pattern,
  ppm=FALSE,
  dmz=0.0001,
  frac=1/4,
  env="Gaussian",
  resolution=1E6,
  plotit=TRUE
)




data(resolution_list)
resmass<-resolution_list[[4]]
data(isotopes)
data(chemforms)
checked<-check_chemform(isotopes,chemforms)
resolution<-getR(checked,resmass,nknots=13,spar=0.1,plotit=TRUE)

checked<-check_chemform(isotopes,chemforms)
checked[,3]<-(checked[,3]/abs(-2))
resolution<-getR(checked,resmass,nknots=13,spar=0.1,plotit=TRUE)




# batch of chemforms #######
data(isotopes)
data(chemforms)
pattern<-isopattern(
  isotopes,
  chemforms,
  threshold=0.1,
  plotit=TRUE,
  charge=FALSE,
  emass=0.00054858,
  algo=1
)
############################
# Single chemical formula ##
data(isotopes)


pattern<-isopattern(
  isotopes,
  "C100H200S2Cl5",
  threshold=0.1,
  plotit=TRUE,
  charge=FALSE,
  emass=0.00054858,
  algo=1
)
############################




data(isotopes)
data(resolution_list)
data(chemforms)
chemforms<-chemforms[1:10]
checked<-check_chemform(
  isotopes,
  chemforms
);
resmass<-resolution_list[[1]]

centro<-isowrap(
  isotopes,
  checked,
  resmass=resolution_list[[4]],
  resolution=FALSE,
  nknots=4,
  spar=0.2,
  threshold=0.1,
  charge=1,
  emass=0.00054858,
  algo=2,
  ppm=FALSE,
  dmz="get", # retrieve dm from R=m/dm
  frac=1/4,
  env="Gaussian",
  detect="valley",
  plotit=TRUE
)



#####
#formula handling
{
  formula_table = function(X){
    Xformula = get.formula(X)
    temp=Xformula@isotopes
    df = data.frame(as.numeric(temp[,2]))
    rownames(df)=temp[,1]
    colnames(df)=Xformula@string
    return(df)
  }
  table_to_formula = function(result_formula){
    r = as.character()
    for(i in 1:nrow(result_formula)){
      if(result_formula$result[i]!=0){
        if(result_formula$result[i]==1){r=paste(r, result_formula$Row.names[i], sep="")}
        else{r=paste(r, result_formula$Row.names[i],  result_formula$result[i], sep="")}
      }
    }
    if(sum(result_formula$result<0)!=0){
      r=paste(r, "Illegal_formula")
    }
    return (r)
  }
  formula_manipulate = function(F1, F2, plusminus) {
    temp_1 = formula_table(F1)
    temp_2 = formula_table(F2)
    result_formula = merge(temp_1,temp_2, all=T,  by="row.names")
    result_formula[is.na(result_formula)]=0
    result_formula["result"]=result_formula[,2]+result_formula[,3]*plusminus
    return (table_to_formula(result_formula))
  }
}

time=Sys.time()
for(i in 1:100){
  #test=get.formula("C2H3O1+CH2O1")@string #5.4s for 1000
  test=check_chemform(isotopes,mergeform("[13]C12C2H3O1","CH2"))$new_formula #0.24s for 1000
  #test=formula_manipulate("C2H3O1","CH2O1",1) #13s for 1000
}
time = Sys.time()-time
time

time=Sys.time()
for(i in 1:100){
  
  test=check_chemform(isotopes,subform("C2H3SAsO1Cl21","C1H2"))$new_formula #12.2s for 1000
  #test=formula_manipulate("C2H3SAsO1Cl21","CH2",-1) #12.5s for 1000
}
time = Sys.time()-time
time


  
formula1<-c("C10[13]C2H10Cl10")
formula2<-c("C2H5Na1D2[18]O1")
mergeform(formula1,formula2)

check_chemform(isotopes,formula1)


####

formula1<-c("C10H8O1Cl0")
formula2<-c("C2H5O2")
data(isotopes)
check_chemform(isotopes, "CH3COOH")


check_chemform(isotopes,mergeform(formula1,formula2),get_sorted = F)

test=subform(formula1,formula2)

test=subform("C2H3O1","C1H2O1Cl2")

#Fastest formula handling.
#Every element has to have a number come after it. Otherwise it will combine two things into one element.


##My_mergefrom&My_subfrom
formula1<-c("C10H8O3Cl0")
formula2<-c("C10H8O3Cl0")
#formula2<-c("C2H5O2")
My_mergefrom=function (formula1, formula2)
{
  formula2 <- gsub("D", "[2]H", formula2)
  ende2 <- nchar(formula2)
  element2 <- c()
  number2 <- c()
  j <- c(1)
  while (j <= ende2) {
    #browser()
    if (substr(formula2, j, j) == c("[")) {
      b <- j
      while (any(substr(formula2, j, j) == c("]")) != 
             TRUE) {
        j <- c(j + 1)
      }
      k <- j
      while (any(substr(formula2, j, j) == c("0", "1", 
                                             "2", "3", "4", "5", "6", "7", "8", "9")) != 
             TRUE) {
        j <- c(j + 1)
      }
      m <- c(j - 1)
      element2 <- c(element2, substr(formula2, b, m))
    }
    if (any(substr(formula2, j, j) == c("0", "1", "2", "3", 
                                        "4", "5", "6", "7", "8", "9")) != TRUE) {
      k <- j
      while (any(substr(formula2, j, j) == c("0", "1", 
                                             "2", "3", "4", "5", "6", "7", "8", "9")) != 
             TRUE) {
        j <- c(j + 1)
      }
      m <- c(j - 1)
      j <- c(j - 1)
      element2 <- c(element2, substr(formula2, k, m))
    }
    if (any(substr(formula2, j, j) == c("0", "1", "2", "3", 
                                        "4", "5", "6", "7", "8", "9")) == TRUE) {
      k <- j
      while (any(substr(formula2, j, j) == c("0", "1", 
                                             "2", "3", "4", "5", "6", "7", "8", "9")) == 
             TRUE) {
        j <- c(j + 1)
      }
      m <- c(j - 1)
      j <- c(j - 1)
      number2 <- c(number2, as.numeric(substr(formula2, 
                                              k, m)))
    }
    j <- j + 1
  }
  formulas <- c()
  for (i in 1:length(formula1)) {
    #i=1
    formula1[i] <- gsub("D", "[2]H", formula1[i])
    ende1 <- nchar(formula1[i])
    element1 <- c()
    number1 <- c()
    j <- c(1)
    while (j <= ende1) {
      if (substr(formula1[i], j, j) == c("[")) {
        b <- j
        while (any(substr(formula1[i], j, j) == c("]")) != 
               TRUE) {
          j <- c(j + 1)
        }
        k <- j
        while (any(substr(formula1[i], j, j) == c("0", 
                                                  "1", "2", "3", "4", "5", "6", "7", "8", "9")) != 
               TRUE) {
          j <- c(j + 1)
        }
        m <- c(j - 1)
        element1 <- c(element1, substr(formula1[i], 
                                       b, m))
      }
      if (any(substr(formula1[i], j, j) == c("0", "1", 
                                             "2", "3", "4", "5", "6", "7", "8", "9")) != 
          TRUE) {
        k <- j
        while (any(substr(formula1[i], j, j) == c("0", 
                                                  "1", "2", "3", "4", "5", "6", "7", "8", "9")) != 
               TRUE) {
          j <- c(j + 1)
        }
        m <- c(j - 1)
        j <- c(j - 1)
        element1 <- c(element1, substr(formula1[i], 
                                       k, m))
      }
      if (any(substr(formula1[i], j, j) == c("0", "1", 
                                             "2", "3", "4", "5", "6", "7", "8", "9")) == 
          TRUE) {
        k <- j
        while (any(substr(formula1[i], j, j) == c("0", 
                                                  "1", "2", "3", "4", "5", "6", "7", "8", "9")) == 
               TRUE) {
          j <- c(j + 1)
        }
        m <- c(j - 1)
        j <- c(j - 1)
        number1 <- c(number1, as.numeric(substr(formula1[i], 
                                                k, m)))
      }
      j <- j + 1
    }
    both <- unique(c(element1, element2))
    both = both[order(both)]
    counts <- c()
    for (i in 1:length(both)) {
      if (any(element1 == both[i])) {
        it1 <- c(number1[element1 == both[i]])
      }
      else {
        it1 <- c(0)
      }
      if (any(element2 == both[i])) {
        it2 <- c(number2[element2 == both[i]])
      }
      else {
        it2 <- c(0)
      }
      counts <- c(counts, it1 + it2)
    }
    formula_all <- ""
    for (i in 1:length(both)) {
      if(counts[i]==0){next}
      formula_all <- paste(formula_all, both[i], counts[i], 
                           sep = "")
    }
    formulas <- c(formulas, formula_all)
  }
  return(formulas)
}
My_mergefrom(formula1,formula2)

My_subfrom =function (formula1, formula2) 
{
  formula2 <- gsub("D", "[2]H", formula2)
  ende2 <- nchar(formula2)
  element2 <- c()
  number2 <- c()
  j <- c(1)
  while (j <= ende2) {
    if (substr(formula2, j, j) == c("[")) {
      b <- j
      while (any(substr(formula2, j, j) == c("]")) != 
             TRUE) {
        j <- c(j + 1)
      }
      k <- j
      while (any(substr(formula2, j, j) == c("0", "1", 
                                             "2", "3", "4", "5", "6", "7", "8", "9")) != 
             TRUE) {
        j <- c(j + 1)
      }
      m <- c(j - 1)
      element2 <- c(element2, substr(formula2, b, m))
    }
    if (any(substr(formula2, j, j) == c("0", "1", "2", "3", 
                                        "4", "5", "6", "7", "8", "9")) != TRUE) {
      k <- j
      while (any(substr(formula2, j, j) == c("0", "1", 
                                             "2", "3", "4", "5", "6", "7", "8", "9")) != 
             TRUE) {
        j <- c(j + 1)
      }
      m <- c(j - 1)
      j <- c(j - 1)
      element2 <- c(element2, substr(formula2, k, m))
    }
    if (any(substr(formula2, j, j) == c("0", "1", "2", "3", 
                                        "4", "5", "6", "7", "8", "9")) == TRUE) {
      k <- j
      while (any(substr(formula2, j, j) == c("0", "1", 
                                             "2", "3", "4", "5", "6", "7", "8", "9")) == 
             TRUE) {
        j <- c(j + 1)
      }
      m <- c(j - 1)
      j <- c(j - 1)
      number2 <- c(number2, as.numeric(substr(formula2, 
                                              k, m)))
    }
    j <- j + 1
  }
  formulas <- c()
  for (i in 1:length(formula1)) {
    i=1
    formula1[i] <- gsub("D", "[2]H", formula1[i])
    ende1 <- nchar(formula1[i])
    element1 <- c()
    number1 <- c()
    j <- c(1)
    while (j <= ende1) {
      if (substr(formula1[i], j, j) == c("[")) {
        b <- j
        while (any(substr(formula1[i], j, j) == c("]")) != 
               TRUE) {
          j <- c(j + 1)
        }
        k <- j
        while (any(substr(formula1[i], j, j) == c("0", 
                                                  "1", "2", "3", "4", "5", "6", "7", "8", "9")) != 
               TRUE) {
          j <- c(j + 1)
        }
        m <- c(j - 1)
        element1 <- c(element1, substr(formula1[i], 
                                       b, m))
      }
      if (any(substr(formula1[i], j, j) == c("0", "1", 
                                             "2", "3", "4", "5", "6", "7", "8", "9")) != 
          TRUE) {
        k <- j
        while (any(substr(formula1[i], j, j) == c("0", 
                                                  "1", "2", "3", "4", "5", "6", "7", "8", "9")) != 
               TRUE) {
          j <- c(j + 1)
        }
        m <- c(j - 1)
        j <- c(j - 1)
        element1 <- c(element1, substr(formula1[i], 
                                       k, m))
      }
      if (any(substr(formula1[i], j, j) == c("0", "1", 
                                             "2", "3", "4", "5", "6", "7", "8", "9")) == 
          TRUE) {
        k <- j
        while (any(substr(formula1[i], j, j) == c("0", 
                                                  "1", "2", "3", "4", "5", "6", "7", "8", "9")) == 
               TRUE) {
          j <- c(j + 1)
        }
        m <- c(j - 1)
        j <- c(j - 1)
        number1 <- c(number1, as.numeric(substr(formula1[i], 
                                                k, m)))
      }
      j <- j + 1
    }
    formula_all <- TRUE
    for (i in 1:length(element2)) {
      if (any(element2[i] == element1) == FALSE) {
        return(F)
        formula_all <- paste(element2[i], " from formula 2 not part of formula1", 
                             sep = "")
      }
      else {
        if (number2[i] > number1[element2[i] == element1]) {
          return(F)
          formula_all <- paste("Atom number of ", element2[i], 
                               " from formula 2 not fully subset of formula1 atom number", 
                               sep = "")
        }
      }
    }
    if (formula_all == TRUE) {
      formula_all <- ""
      counts <- c()
      for (i in 1:length(element1)) {
        if (any(element2 == element1[i])) {
          counts <- c(counts, number1[i] - (number2[element2 == 
                                                      element1[i]]))
        }
        else {
          counts <- c(counts, number1[i])
        }
      }
      element1 <- element1[counts > 0]
      counts <- counts[counts > 0]
      for (i in 1:length(counts)) {
        formula_all <- paste(formula_all, element1[i], 
                             counts[i], sep = "")
      }
    }
    formulas <- c(formulas, formula_all)
  }
  return(formulas)
}
My_subfrom(formula1, formula2)


time = Sys.time()
for(i in 1:10000){
  My_subfrom(formula2, formula1)
  
  
}
time = Sys.time()-time
time



