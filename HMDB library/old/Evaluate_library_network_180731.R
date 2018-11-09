#显示中文 
#Sys.setlocale(category = "LC_ALL", locale = "Chinese")
#import library
{
  library(readr)
  #install.packages("rcdk")
  library(rcdk)
  library(igraph)
  #install.packages("fitdistrplus")
  library("fitdistrplus")
  library(tidyr)
  library(dplyr)
  
  #install.packages("lpSolve")
  #library(lpSolve)
  #install.packages("Rglpk")
  library(Rglpk)
  
}

##Read basic data information##
time = data.frame(Sys.time())
colnames(time)="read_data"


#Read data
{
  #setwd("C:/Users/lc8/Dropbox/HRMass_ID/")
  setwd("C:/Users/Li Chen/Dropbox/HRMass_ID/180629 Triwizard Tournament/HMDB library")
  data_known <- read_csv("hmdb_unique.csv")
  hmdb_neutralformula= read_csv("hmdb_full.csv")
  hmdb_endogenous = read_csv("HMDB_Extracted_0801.csv")
}

test=merge(hmdb_neutralformula,hmdb_endogenous, by.x = "HMDB_ID", by.y = "Accession Number", all=T)

#Reorganize datatable
{
  raw = data_known[data_known$Exact_Mass>70&data_known$Exact_Mass<1500,]
  raw = raw[with(raw, order(Exact_Mass)),]
  node_list = data.frame(c(1:nrow(raw)),raw$Exact_Mass, raw$MF, raw$Name, raw$category, stringsAsFactors=FALSE)
  colnames(node_list)=c("ID", "mz","MF","Name", "category")
  raw_rnum = nrow(raw)
  merge_node_list=node_list
  merge_nrow=raw_rnum
  merge_library = merge(merge_node_list, hmdb_endogenous, by.x ="MF", by.y="Formula", all=T)
  merge_library = merge_library[!is.na(merge_library$ID),]
}



##find mass difference corresponding to a functional group, such as "CH2"
time["edge_list"]=Sys.time()

#function group information
{
  fun_group_1 = data.frame(c("H2",
                             "C2H2O", #ACETYL
                             "C10H12N5O6P", #AMP
                             "C16H30O", #Palmityl
                             "C6H10O5", #Glycosyl
                             "C11H17NO9", #Sialic acid
                             "HPO3",
                             # "CO", #FORMYL
                             # "NH3",
                             # "NH3-O" ,#transaminase
                             # "OH-NH2", #amino group
                             # "C",
                             # "CH2O",
                             "C9H13N3O10P2", #CDP
                             "S",
                             "SO3",
                             "C2H4",
                             "O",
                             "NH",
                             "CH2",
                             "CO2",
                             "H2O"
                             ),
                           stringsAsFactors=F)
  for (i in 1:nrow(fun_group_1)){
    fun_group_1[i,2]=get.formula(fun_group_1[i,1])@mass
  }
  colnames(fun_group_1) = c("fun_group", "mass")
  fun_group_1[fun_group_1$fun_group=="NH3-O",2]=get.formula("NH3")@mass-get.formula("O")@mass
  fun_group_1[fun_group_1$fun_group=="OH-NH2",2]=get.formula("OH")@mass-get.formula("NH2")@mass
  fun_group_1=rbind(c("Same", 0),fun_group_1)
  fun_group_1$mass=as.numeric(fun_group_1$mass)
}


#Edge list
mass_tol = 0.001
merge_node_list = merge_node_list[with(merge_node_list, order(mz)),]
{
  edge_ls = list()
  temp_mz_list=merge_node_list$mz
  for (k in 1:nrow(fun_group_1)){
    temp_fg=fun_group_1[k,2]
    i=j=1
    temp_edge_list = data.frame(node1=as.numeric(), node2=as.numeric(), linktype=as.numeric(), mass_dif=as.numeric())
    while(i<=merge_nrow){
      temp_ms=0
      while(1){
        j=j+1
        if(j>merge_nrow){break}
        temp_ms = temp_mz_list[j]-temp_mz_list[i]
        if(abs(temp_ms-temp_fg)<mass_tol){
          temp_edge_list[nrow(temp_edge_list)+1,]=c(merge_node_list$ID[i], merge_node_list$ID[j], k, (temp_ms-temp_fg))
        }
        if((temp_ms-temp_fg)>mass_tol){break}
      }
      i=i+1
      
      while(j>i){
        j=j-1
        temp_ms = temp_mz_list[j]-temp_mz_list[i]
        if((temp_ms-temp_fg)<(-mass_tol)){break}
      }
    }
    edge_ls[[k]]=temp_edge_list
  }
}
edge_list = bind_rows(edge_ls)


fun_group_1["summary"]=table(edge_list$linktype)

#Scoring edge based on mass accuracy
{
edge_list = edge_list[edge_list$node1!=edge_list$node2,]
edge_mzdif_FIT <- fitdist(as.numeric(edge_list[,4])*10^6, "norm")    
hist(edge_list$mass_dif)

plot(edge_mzdif_FIT)  
edge_list["edge_massdif_score"]=dnorm(edge_list[,4]*10^6, edge_mzdif_FIT$estimate[1], edge_mzdif_FIT$estimate[2])
edge_list["edge_massdif_score"]=edge_list["edge_massdif_score"]/max(edge_list["edge_massdif_score"])
}


##Predict formula based on known formula and edgelist
time["formula_pred"]=Sys.time()
#Initialize predict_formula from HMDB known formula
{
  data_str = data.frame(id=as.numeric(),
                        formula=as.character(),
                        degree=as.numeric(),
                        parent=as.numeric(),
                        score=as.numeric(), 
                        parent_formula = as.character(), stringsAsFactors = F)
  #sf stands for summary_formula
  sf = lapply(1:nrow(node_list),function(i)(data_str))
  Initial_formula =sf[[1]]
  Initial_formula[1,]= list(1,"C3H2O2",0,1,1,"C3H2O2")
  sf[[1]]=Initial_formula

  merge_node_list["Predict_formula"]=NA
  merge_node_list["flag_head"]=-1
  merge_node_list["flag_tail"]=-1

  merge_node_list[1,6]="C3H2O2"
  merge_node_list[1,7]=1
  merge_node_list[1,8]=1
  mnl = merge_node_list[with(merge_node_list, order(ID)),]
}
##Predict formula based on known formula and edgelist & update edge list

#formula_manipulate is used to process two formula, calling subfunction formula_table and table_to_formula
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
formula_manipulate = function(F1, F2, plusminus){
  temp_1 = formula_table(F1)
  temp_2 = formula_table(F2)
  result_formula = merge(temp_1,temp_2, all=T,  by="row.names")
  result_formula[is.na(result_formula)]=0
  result_formula["result"]=result_formula[,2]+result_formula[,3]*plusminus
  return (table_to_formula(result_formula))
}
}
#while loop to Predict formula based on known formula and edgelist
flag = 1
while(flag==1){
  flag = 0
  for (head in 1:nrow(mnl)){
    if(mnl$flag_head[head]!=1){next}
    mnl$flag_head[head]=0

    temp_edge_list=subset(edge_list, edge_list$node1==head)
    if(nrow(temp_edge_list)==0){next}

    if(nrow(sf[[head]])==0){
      print(paste("Caution. Head is flagged 1 without formula imput. Head", head))
      next
    }
    if(sf[[head]]$score[1]==0){next}

    for(i in 1:nrow(temp_edge_list)){
      tail=temp_edge_list$node2[i]
      temp_fg = fun_group_1$fun_group[temp_edge_list$linktype[i]]
      if(temp_fg=="Same"){temp_formula=sf[[head]]$formula[1]}
      else{temp_formula = get.formula(paste(sf[[head]]$formula[1],"+", temp_fg))@string}
      temp_degree = sf[[head]]$degree[1]+1
      temp_parent = head
      temp_parent_formula = sf[[head]]$formula[1]
      temp_score = sf[[head]]$score[1]*temp_edge_list$edge_massdif_score[i]

      #Criteria to enter new entry into formula list
      #1. new formula
      temp = sf[[tail]]
      temp_subset=subset(temp, temp$formula==temp_formula)
      if(nrow(temp_subset)!=0){
        #2. smaller degree or much higher scores
        if((temp_degree>=min(temp_subset$degree))&(temp_score<=(1.2*max(temp_subset$score)))){
          next
        }
      }

      #Enter new entry
      temp[nrow(temp)+1, ] = list(tail, temp_formula, temp_degree, temp_parent, temp_score, temp_parent_formula)
      temp = temp[with(temp, order(-score)),]
      sf[[tail]]=temp
      mnl$flag_head[tail]=mnl$flag_tail[tail]=1
      flag=1
    }
  }

  for (tail in 1:nrow(mnl)){
    if(mnl$flag_tail[tail]!=1){next}
    mnl$flag_tail[tail]=0

    temp_edge_list=subset(edge_list, edge_list$node2==tail)
    if(nrow(temp_edge_list)==0){next}

    if(nrow(sf[[tail]])==0){
      print(paste("Caution. Tail is flagged 1 without formula imput. Tail"), tail)
      next
    }
    if(sf[[tail]]$score[1]==0){next}
    for(i in 1:nrow(temp_edge_list)){
      head=temp_edge_list$node1[i]
      temp_fg = fun_group_1$fun_group[temp_edge_list$linktype[i]]
      if(temp_fg=="Same"){temp_formula=sf[[tail]]$formula[1]}
      else{temp_formula = formula_manipulate(sf[[tail]]$formula[1], temp_fg, -1)}
      temp_degree = sf[[tail]]$degree[1]+1
      temp_parent_formula = sf[[tail]]$formula[1]
      temp_parent = tail
      if(grepl("Illegal_formula", temp_formula)){temp_score=0}
      else{temp_score = sf[[tail]]$score[1]*temp_edge_list$edge_massdif_score[i]}

      #Criteria to enter new entry into formula list
      #1. new formula
      temp = sf[[head]]
      temp_subset=subset(temp, temp$formula==temp_formula)
      if(nrow(temp_subset)!=0){
        #2. smaller degree or much higher scores
        if((temp_degree>=min(temp_subset$degree))&(temp_score<=(1.2*max(temp_subset$score)))){
          next
        }
      }

      #Output step
      temp[nrow(temp)+1,] = list(head, temp_formula, temp_degree, temp_parent, temp_score, temp_parent_formula)
      temp = temp[with(temp, order(-score)),]
      sf[[head]]=temp
      mnl$flag_head[head]=mnl$flag_tail[head]=1
      flag=1
    }
  }
}


#Update the mnl from [List]
for(i in 1:nrow(mnl)){
  mnl$Predict_formula[[i]]=sf[[i]]$formula[1]
}



merge_node_list = mnl


#Debug & test
{
merge_node_list_pred=merge_node_list[complete.cases(merge_node_list[,c(3,6)]),]
pred_dif = merge_node_list_pred[merge_node_list_pred$MF!=merge_node_list_pred$Predict_formula,]
pred_dif_cat1 = pred_dif[pred_dif$category==1,]

cat2 = merge_node_list[merge_node_list$category==2&(is.na(merge_node_list$Predict_formula)),]
cat3 = merge_node_list[merge_node_list$category==3,]
cat2 = merge_node_list[merge_node_list$category==2,]
illegal_formula = merge_node_list[grepl("Illegal_formula",merge_node_list$Predict_formula),]
}

sf_unlist = bind_rows(sf)

write.csv(merge_node_list, "HMDB_node_list2.csv", row.names = F)
write.csv(edge_list, "HMDB_edge_list.csv", row.names = F)
write.csv(sf_unlist,"HMDB_predicted_formula.csv",row.names = F)
# lapply(mylist, function(x) write.table( data.frame(x), 'test.csv'  , append= T, sep=',' ))


##Network analysis
time["network_analysis"]=Sys.time()
#Generate network graph using node_list and edge_list above
merge_edge_list = edge_list[edge_list$edge_massdif_score>0.99,]
g <- graph_from_data_frame(d = merge_edge_list, vertices = merge_node_list, directed = FALSE)
colors <- c("white", "red", "orange", "blue", "dodgerblue", "cyan")
V(g)$color = colors[merge_node_list$category+1]
E(g)$color = colors[merge_edge_list$linktype+1]
merge_node_list["degree"]=degree(g, mode = c("all"))

edge_size= edge_membership = list()
for (i in unique(merge_edge_list$linktype)){
temp_edge_list=merge_edge_list[merge_edge_list$linktype!=i,]
temp_g <- graph_from_data_frame(d = temp_edge_list, vertices = merge_node_list, directed = FALSE)
temp_clu = components(temp_g)
edge_size[[i]]=table(temp_clu$csize)
edge_membership[[i]]=temp_clu$membership
}
#test indicates which cluster the nodes belongs to
test=bind_cols(edge_membership)
test["ID"]=1:nrow(test)
edge_size[[5]]

fun_group_1["drop_out_size"]=unlist(lapply(edge_size, function(x) max(as.numeric(names(x)))))
fun_group_1["drop_out_dif"]=fun_group_1["drop_out_size"]-fun_group_1[1,"drop_out_size"]

##Drop out CO

test2 = test[test$V1!=test$V5&test$V1==1,]

test3=merge_library[merge_library$ID %in% test2$ID,]

#subnetwork criteria 
clu = components(g)
main_network = igraph::groups(clu)[table(clu$membership)>1000]
subnetwork = igraph::groups(clu)[table(clu$membership)>30&table(clu$membership)<1000]
g_subnetwork_list = lapply(subnetwork, make_ego_graph, graph=g, order=diameter(g), mode="all")
sub_nodelist = merge_library[merge_library$ID %in% subnetwork[[2]],]
#write.csv(sub_nodelist, "maybe wrong entry in HMDB.csv")

g_sub = graph_from_data_frame(d = edge_list_sub, vertices = merge_node_list[merge_node_list$ID %in% c(edge_list_sub[,1], edge_list_sub[,2]),], directed = FALSE)
colors <- c("white", "red", "orange", "blue", "dodgerblue", "cyan")
V(g_sub)$color = colors[vertex.attributes(g_sub)$category+1]
E(g_sub)$color = colors[edge_list_sub$category+1]

#Basic graph characteristics, distance, degree, betweeness
{
#distance
farthest_vertices(g) 
#Degree
g.degree <- degree(g, mode = c("all"))
table(g.degree)
hist(g.degree)
# which.max(g.degree)
node_list_0degree = node_list[g.degree==0,]
# node_list$Exact_Mass[g.degree==0]
# node_list$Exact_Mass[71]
# node_list$MF[71]
#Betweenness
g.b <- betweenness(g, directed = T)
colors <- c("white", "red", "orange", "blue", "dodgerblue", "cyan")
V(g)$color <- colors[merge_node_list$category+2]
plot(g, 
     vertex.label = NA,
     edge.color = 'black',
     vertex.size = sqrt(degree(g, mode = c("all")))+1,
     edge.arrow.size = 0.05,
     layout = layout_nicely(g))
}

#Analyze the network/subgraph of specific node
{
interested_node = "1941"
g_intrest <- make_ego_graph(g,2, nodes = interested_node, mode = c("all"))[[1]]
dists = distances(g_intrest, interested_node)
colors <- c("black", "red", "orange", "blue", "dodgerblue", "cyan")

#V(g_intrest)$color <- colors[dists+1]
plot(g_intrest,
     #vertex.color = 'white',
     #vertex.label = vertex.attributes(g_intrest)$Predict_formula,
     #vertex.label = vertex.attributes(g_intrest)$MF,
     vertex.label = vertex.attributes(g_intrest)$MF,
     #vertex.label = vertex.attributes(g_intrest)$ID,
     vertex.label.color = "red",
     vertex.label.cex = 1,
     edge.color = 'black',
     edge.label = fun_group_1$fun_group[edge.attributes(g_intrest)$linktype],
     vertex.size = 10,
     edge.arrow.size = .05,
     main = paste("Subnetwork of node", interested_node)
)
}

merge_node_list[which(merge_node_list$MF=="C6H13O4PS2"),]
merge_node_list[which(merge_node_list$Predict_formula=="C15H4O5"),]

#Analysis of Subnetwork 
{
clu=components(g)
#table(clu$csize)
#subnetwork criteria 
main_network = igraph::groups(clu)[table(clu$membership)>1000]
subnetwork = igraph::groups(clu)[table(clu$membership)>30&table(clu$membership)<1000]
g_subnetwork_list = lapply(subnetwork, make_ego_graph, graph=g, order=diameter(g), mode="all")
sub_nodelist = merge_library[merge_library$ID %in% subnetwork[[1]],]
#write.csv(sub_nodelist, "maybe wrong entry in HMDB.csv")


for (i in 1:length(subnetwork)){
  plot(g_subnetwork_list[[i]][[1]],
       #vertex.color = 'white',
       vertex.label = vertex.attributes(g_subnetwork_list[[i]][[1]])$MF,
       #vertex.label = vertex.attributes(g_subnetwork_list[[i]][[1]])$Exact_Mass,
       vertex.label.color = "red",
       vertex.label.cex = 1,
       edge.color = 'black',
       edge.label = fun_group_1$fun_group[edge.attributes(g_subnetwork_list[[i]][[1]])$linktype],
       vertex.size = 10,
       edge.arrow.size = .05,
       main = paste("Subnetwork",names(subnetwork)[[i]])
  )
}
merge_node_list[subnetwork[["52"]],]
}



Match_HMDB = function(Exact_mass, data_known, mode){
H_mass = 1.00782503224
e_mass = 0.00054857990943
temp_mass = Exact_mass*abs(mode) - (H_mass-e_mass)*mode
temp_table = data_known[abs(data_known$Exact_Mass-temp_mass)<0.002,]
return(temp_table)
}


# #Interactive network w/ threejs
# install.packages("htmlwidgets ")
# library(threejs)
# library(htmlwidgets)
# 
# graphjs(g, main="Network!", bg="gray10", showLabels=F, stroke=F, 
#         curvature=0.1, attraction=0.9, repulsion=0.8, opacity=0.9)


time["end"]=Sys.time()
time_used=as.numeric()
for (i in 2:ncol(time)){
time_used[i-1] = round(time[1,i]-time[1,i-1],2)
}
names(time_used)=colnames(time)[1:(ncol(time)-1)]
time_used



# #Calculate mz difference
# {
# merge_node_list["MF_measured_mz_dif"]=as.numeric()
# merge_node_list["Predict_measured_mz_dif"]=as.numeric()
# for(i in 1: nrow(node_list)){
#   if(!is.na(merge_node_list$`MF`[i])){
#     merge_node_list$`MF_measured_mz_dif`[i]=get.formula(merge_node_list$MF[i])@mass-merge_node_list$Exact_Mass[i]
#   }
#   if((!is.na(merge_node_list$`Predict_formula`[i]))&(!grepl("Illegal_formula", merge_node_list$Predict_formula[i]))){
#     merge_node_list$`Predict_measured_mz_dif`[i]=get.formula(merge_node_list$Predict_formula[i])@mass-merge_node_list$Exact_Mass[i]
#   }
# }
# 
# data_analysis = merge_node_list[merge_node_list$category!=0&!(is.na(merge_node_list$MF)),]
# data_analysis = data_analysis[with(data_analysis, order(`MF_measured_mz_dif`)),]
# test_logic=data_analysis$MF!=data_analysis$Predict_formula
# test_logic[is.na(test_logic)]=T
# data_analysis_dif = data_analysis[test_logic,]
# }
#Plot machine measurement difference
# {
# cat1 = merge_node_list[merge_node_list$category==1,]
# pred_same_cat1 = cat1[(cat1$MF==cat1$Predict_formula)&(!is.na(cat1$Predict_formula)),]
# hist(abs(pred_same_cat1$MF_measured_mz_dif)*10^3)
# hist(abs(pred_same_cat1$MF_measured_mz_dif/pred_same_cat1$Exact_Mass)*10^6)
# }
# 
