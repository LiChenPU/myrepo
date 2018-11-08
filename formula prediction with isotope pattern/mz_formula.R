

library(readr)
library(rcdk)
library(dplyr)
library(enviPat)

#Read data
{
  setwd(dirname(rstudioapi::getSourceEditorContext()$path))
  getwd()
  raw_node_list=read_csv("merge_node_list.csv")
  raw_pred_formula=read_csv("All_formula_predict.csv")
}

#Clean up
{
  lib_nodes = raw_node_list[raw_node_list$category==0,]
  lib_nodes_cutoff = nrow(raw_node_list)-nrow(lib_nodes)
  unknown_nodes = raw_node_list[raw_node_list$category!=0,]
  unknown_nodes_mz = unknown_nodes$mz
  
  pred_formula=raw_pred_formula[!grepl("-",raw_pred_formula$formula),]
  
  sf=list()
  for(i in 1:nrow(raw_node_list)){
    sf[[i]]=pred_formula[pred_formula$id==i,]
  }
}

#Parameter
{
mode = -1
H_mass = 1.00782503224
e_mass = 0.00054857990943
mz_input = 1120.4598
mz_neutral = mz_input - (H_mass-e_mass)*mode


C_range = 1:99
H_range = 1:100
O_range = 0:20
N_range = 0:12
Cl_range = 0:3
P_range = 0:5
S_range = 0:5
Na_range = 0:0
K_range = 0:0
F_range = 0:0
Br_range = 0:0
I_range = 0:0
Si_range = 0:0

db_min = -1
db_max = 99
H_min = 0
C_min = 0
C_max = 99



N_rule = 1


H_iso=(1.00782503224-1)
N_iso=(14.0030740048-14)
O_iso=(15.99491462-16)
Cl_iso=(34.96885268-35)
F_iso=(18.99840322-19)
Br_iso=(78.9183371-79)
I_iso=(126.904473-127)
S_iso=(31.97207100-32)
P_iso=(30.97376163-31)
Si_iso=(27.97692653-28)
B_iso=(10.0129370-10)
Na_iso=(22.9897692809-23)
K_iso=(38.96370668-39)
M_ELECTRON=0.00054857990943
}

#my function to process formula 
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
      if(result_formula[i,1]!=0){
        r=paste(r, rownames(result_formula)[i],  result_formula[i,1], sep="")
      }
    }
    if(sum(result_formula$result<0)!=0){
      r=paste(r, "Illegal_formula")
    }
    return (r)
  }
}

time = Sys.time()
iteration=0


mz_neutral = 307.0838
#mz_neutral = 506.9958
mz_neutral = 181.0739
#mz_neutral = 1115.3088
mz_neutral = 745.0911
mz_neutral = 127.9430




##Core code calculate formula given mz
mz_formula = function(mz_neutral, ppm)
{
  
  
  mz_integer = floor(mz_neutral+0.2)
  mz_decimal = mz_neutral - mz_integer
  tolerance = mz_neutral*ppm/10^6
  
  temp_formula = formula_table("CHNOClSPNaKFBrISi")
  temp_formula = temp_formula[order(row.names(temp_formula)), , drop=FALSE]
  temp_formula_list = list()
  temp = data.frame(formula = "formula", differ = as.numeric(0), db_r=as.numeric(0), stringsAsFactors=F)
  
  for(Si in Si_range){
  for(I in I_range){
  for(Br in Br_range){
  for(fluorine in F_range){
      for(K in K_range){
        for(Na in Na_range){
          if((K + Na) > 1) next
          for(P in P_range){
            for(S in S_range){
              for(Cl in Cl_range){
                CHNO = mz_integer - 35*Cl - 32*S - 31*P - 23*Na - 39*K
                - 19*fluorine - 79*Br - 127*I - 28*Si 
                if(CHNO<0) break
                for(N in N_range){
                  #(1) N rule
                  if(N_rule==1){
                    if(N%%2!=mz_integer%%2) next
                  }
                  for(O in O_range){
                    CH = CHNO - 14*N - 16*O 
                    
                    if(CH<0) break
                    else {
                      Heavy_atom_decimal = (Cl_iso*Cl) + (F_iso*fluorine) + (Br_iso*Br) + (I_iso*I) + (S_iso*S) + 
                        (P_iso*P) + (Si_iso*Si) + (Na_iso*Na) + (K_iso*K) + (N_iso*N) + (O_iso*O) 
                      H_min = max(0, floor((mz_decimal-Heavy_atom_decimal-tolerance)/H_iso))
                      C_lower = max(C_min, floor(CH/14)-1)
                      C_upper = min(C_max, floor(CH/12)+1)
                    }
                    for(C in C_lower:C_upper){
                      
                      iteration = iteration+1

                      ##(2) H number
                      H = CH - 12*C 
                      if(H < H_min) next
                      
                      ##(3) Saturation, ring&double bond
                      db_r = (C+1) - ((H+Cl+fluorine+Br+I+Na+K) - (N+P))/2;
                      if(db_r < db_min | db_r >db_max) next
                      
                      ##(4) Decimal place
                      differ = (H_iso*H) + Heavy_atom_decimal - mz_decimal
                      if(differ > tolerance | differ < -tolerance) next
                      
                      ## record

                      temp_formula[,1] = c(Br, C, Cl, fluorine, H, I, K, N, Na, O, P, S, Si)
                      temp[1,] = list(table_to_formula(temp_formula),differ,db_r)
                      temp_formula_list[[length(temp_formula_list)+1]]=temp
}}}}}}}}}}}}
          
formula_df = bind_rows(temp_formula_list)
return(formula_df)
}

ppm = 10
timer = Sys.time()

for(i in 1:nrow(unknown_nodes)){
#for(i in 1:10){
  if(i %% 50==0){
    print(paste("i =", i))
    print(Sys.time()-timer)
  }
  mz_neutral = unknown_nodes_mz[i]
  mz_formula_df = mz_formula(mz_neutral, ppm)
  temp = sf[[i]]
  
  sf[[i]] = temp[sf[[i]]$formula %in% mz_formula_df$formula,]
}
print(Sys.time()-timer)


pred_formula_prune = bind_rows(sf)
write_csv(pred_formula_prune,"pred_formula_prune.csv")




#####
#Calculate isotope pattern
{
data(isotopes)

pattern<-isopattern(
  isotopes,
  formula_df$formula,
  threshold=0.1,
  plotit=F,
  charge=FALSE,
  emass=0.00054858,
  algo=2
)
pattern$C41Cl5H66N9O16

profiles<-envelope(
  pattern,
  ppm=F,
  dmz="get",
  env="Gaussian",
  resolution = 10000,
  plotit=F
)
centro <- vdetect(profiles, detect = "centroid", plotit = F, verbose = F)
centro_origin = centro
report_isotope = 3
iso_mat = matrix(0,nrow = length(centro), ncol = report_isotope+1)
flagi=0

for(i in 1:length(centro)){
  
  for(j in 2:nrow(centro[[i]])){
    if((centro[[i]][j,1]-centro[[i]][j-1,1])<0.8|(centro[[i]][j,1]-centro[[i]][j-1,1])>1.8){
      flagi=i
      print("warning")
    }
    
  }
    
  n = min(nrow(centro[[i]]),report_isotope+1)
  iso_mat[i,1:n]=centro[[i]][1:n,2]
}
iso_df = data.frame(iso_mat)
colnames(iso_df)[1]="Parent"
for(i in 1:report_isotope){
  colnames(iso_df)[i+1]=paste("M+",i,sep="")
}
formula_iso = cbind(formula_df, iso_df)
formula_iso[,which(colnames(formula_iso)=="Parent"):ncol(formula_iso)] = sweep(formula_iso[,which(colnames(formula_iso)=="Parent"):ncol(formula_iso)],
                                                                                1,
                                                                                formula_iso[,which(colnames(formula_iso)=="Parent")],
                                                                                FUN = "/")

plot(formula_iso$`M+1`,formula_iso$`M+2`)
#plot(formula_iso$`M+1`,formula_iso$`M+2`, xlim = range(0,0.22), ylim = range(0,.3))

test = formula_iso[formula_iso$`M+2`<0.3&formula_iso$`M+1`<0.22,]
test2 = formula_iso[formula_iso$`M+2`<0.3&formula_iso$`M+1`>0.22,]


#write_excel_csv(pattern,"pattern.csv")


}


time = Sys.time()-time
time
# 
# time = Sys.time()
# for(i in 1:1000000){
#   min(C_range)
#   #nrow(centro[[1]])
#   #length(centro[[1]][,2])
#   # for(i in 1:length(test)){
#   #
#   #   test[[i]]["formula"]=names(test[i])
#   # }
# }
# time = Sys.time()-time
# time

##
# indx <- sapply(centro,length)
# res <- as.data.frame(do.call(rbind,lapply(centro, `length<-`,max(indx))))