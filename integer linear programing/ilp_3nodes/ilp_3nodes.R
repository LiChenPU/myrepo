Sys.setlocale(category = "LC_ALL", locale = "Chinese")

library(rcdk)
library(Matrix)
library(dplyr)
#install.packages("slam")
library(slam)
library(Rglpk)
## Simple mixed integer linear program.
## maximize: 3 x_1 + 1 x_2 + 3 x_3
## subject to: -1 x_1 + 2 x_2 + x_3 <= 4
## 4 x_2 - 3 x_3 <= 2
## x_1 - 3 x_2 + 2 x_3 <= 3
## x_1, x_3 are non-negative integers
## x_2 is a non-negative real number

{
obj <- c(rep(0,4), .1,0.2)

num_lib_nodes = 2
num_unknown_nodes = 1
num_unknown_formula = list()
num_unknown_formula[[1]] = 2
num_edges = 2

node_ls = node_list[1:3,]
edge_ls = edge_list[1:2,]

edge_ls$node1=c(1,2)
edge_ls$node2=c(3,3)
edge_ls=


triplet_library_node = data.frame(i=1:num_lib_nodes, j=1:num_lib_nodes, v=1)
triplet_unknown_nodes_ls = list()
temp_j=num_lib_nodes+1
for(n in 1:num_unknown_nodes){
  temp_i = num_lib_nodes+n
  triplet_unknown_nodes_ls[[n]]=data.frame(i=rep(temp_i, num_unknown_formula[[n]]), j=temp_j:(temp_j+num_unknown_formula[[n]]-1), v=1)
  temp_j = temp_j+num_unknown_formula[[n]]
}
triplet_unknown_node = bind_rows(triplet_unknown_nodes_ls)
triplet_edge_ls_edge = data.frame(i=(1:num_edges)+num_unknown_nodes+num_lib_nodes,
                                  j=(1:num_edges)+temp_j-1,
                                  v=2)
triplet_edge_ls_node = data.frame(i=rep((1:num_edges)+num_unknown_nodes+num_lib_nodes,2),
                                  #j=c(edge_ls$node1,edge_ls$node2),
                                  j=c(1,2,3,4),
                                  v=-1)
triplet_df = rbind(triplet_library_node,
                   triplet_unknown_node,
                   triplet_edge_ls_edge,
                   triplet_edge_ls_node
                   )
mat = simple_triplet_matrix(i=triplet_df$i,
                            j=triplet_df$j,
                            v=triplet_df$v)


dir <- c(rep("==",num_lib_nodes+num_unknown_nodes), rep("<=", num_edges))

rhs = c(rep(1,num_lib_nodes+num_unknown_nodes),rep(0,num_edges))

types <- c(rep("B",6))
max <- TRUE
Rglpk_solve_LP(obj, mat, dir, rhs, types = types, max = max)

}

{
  env <- openEnvCPLEX()
  prob <- initProbCPLEX(env)
  chgProbNameCPLEX(env, prob, "sample")
  nc <- max(mat$j)
  nr <- max(mat$i)
  obj <- c(rep(0, nrow(unknown_formula)), edge_info_sum$edge_score)
  rhs = c(rep(1,num_unknown_nodes),rep(0,nrow(edge_info_sum)))
  sense <- c(rep("E",num_unknown_nodes), rep("L", nrow(edge_info_sum)))
  lb <- rep(0, nc)
  ub <- rep(2, nc)
  # cn <- c("x1", "x2", "x3")
  # rn <- c("q1", "q2", "q3")
  
  
  triplet_df=triplet_df[with(triplet_df,order(i)),]
  cnt=as.vector(table(triplet_df$i))
  beg=vector()
  beg[1]=0
  for(i in 2:length(cnt)){beg[i]=beg[i-1]+cnt[i-1]}
  ind=triplet_df$i-1
  val = mat$v
  
  
  
  
  copyLpwNamesCPLEX(env, prob, nc, nr, CPX_MAX, obj, rhs, sense,
                    beg, cnt, ind, val, lb, ub, NULL, NULL, NULL)
  lpoptCPLEX(env, prob)
  test=solutionCPLEX(env, prob)
  writeProbCPLEX(env, prob, "prob.lp")
  
  lp <- initProbCPLEX(env)
  readCopyProbCPLEX(env, lp, "prob.lp")
  
  delProbCPLEX(env, prob)
  closeEnvCPLEX(env)
}



