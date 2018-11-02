library(rcdk)
library(dplyr)
library(enviPat)


Em = 307.0838
Em = 506.9958
ppm = 5
tolerance = Em*ppm/10^6

Element_list = list(list("C",1,99),
                    list("H",1,99),
                    list("N",0,20),
                    list("O",0,12),
                    list("P",0,5),
                    list("S",0,5),
                    list("Cl",0,5),
                    list("F",0,0),
                    list("Br",0,0),
                    list("I",0,0),
                    list("K",0,0),
                    list("Na",0,0)
                    )



####Do not change from here###

Elem_ls = data.frame(Elements = as.character(), 
                     Nmin = as.numeric(), 
                     Nmax = as.numeric(),
                     Mass = as.numeric(),
                     stringsAsFactors = F
                     )

for(i in 1:length(Element_list)){
  Elem_ls[i,]=Element_list[[i]]
  Elem_ls$Mass[i]= get.formula(Elem_ls$Elements[i])@mass
  
}
Elem_ls$Mass = as.numeric(Elem_ls$Mass)
Elem_ls = Elem_ls[order(Elem_ls$Mass, decreasing = T),]
Elem_ls["Mass_int"]=floor(Elem_ls$Mass+0.2)
Elem_ls["Mass_deci"]=Elem_ls$Mass-Elem_ls["Mass_int"]
Elem_ls["Ratio_deci/int"]=Elem_ls["Mass_deci"]/Elem_ls["Mass_int"]

for(i in 1:length(Element_list)){
  Element_list[[i]]=as.list(t(Elem_ls[i,]))
}


formula = ""
formula_list = list()
formula_str = data.frame(formula = "", 
                         em_error = 0,
                         stringsAsFactors = F)

decompose_formula = function (Element_list, Em, formula){
  
  Element_info = Element_list[[1]]
  mass = as.numeric(Element_info[4])
  
  end_of_recursion=F
  Element_list_corr = Element_list [-1]
  if(length(Element_list_corr) == 0){end_of_recursion=T}
  if(end_of_recursion){
    Nmin = floor((Em-tolerance)/mass)
    Nmax = floor((Em+tolerance)/mass)
  }
  else{
    Nmin = Element_info[[2]]
    Nmax = min(Element_info[[3]], floor((Em+tolerance)/mass))
  }
  
  formula_in = paste(formula,Element_info[[1]],sep="")
  
  for(i in Nmin:Nmax){
    
    #browser()
    formula_out = paste(formula_in,i,sep="")
    
    Em_corr = Em - mass * i
    if(Em_corr<(-tolerance)){
      break
    }
    if(abs(Em_corr)<tolerance){
      formula_str[1,] = list(formula_out, Em_corr)
      formula_list[[length(formula_list)+1]] <<- formula_str
      return(1)
    }
    if(!end_of_recursion){
      decompose_formula(Element_list_corr,Em_corr,formula_out)
    }
  }
  #browser()
  return(1)
}


time = Sys.time()
iteration=0
test=decompose_formula(Element_list, Em, formula)
result_formula = bind_rows(formula_list)


time = Sys.time()-time
time

# 
# 
# time = Sys.time()
# for(i in 1:1000000){
#   #mass = Elem_ls$Mass[1] #0.7s/100k
#   #Elem_ls_corr = Elem_ls [-1,] #4.7s/100k
#   #min(Elem_ls$Nmax[1], floor((Em+tolerance)/mass))
#   #formula_in = paste(formula,Elem_ls$Elements[1],sep="")
#   #Elem_ls_corr= Element_list[-1]
#   #test2=test$`1`
#   
# }
# 
# time = Sys.time()-time
# time

