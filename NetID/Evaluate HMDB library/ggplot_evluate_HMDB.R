library(ggplot2)
library(ggrepel)
library(ggpubr)
library(igraph)
library(RColorBrewer)

# Color setup ##
# display.brewer.pal(4, "Set3")
my_palette = c(brewer.pal(4, "Set3"), rep("#666666", 50))

colfunc = colorRampPalette(c(my_palette[3], "white"))
colfunc(5)[1]

## Figure 1B ####
# {
#   # Use all nodes and all correct connects
#   g = g_basic
#   
#   clu = components(g)
#   edge_size=table(clu$csize)
#   
#   
#   mainnetwork = igraph::groups(clu)[table(clu$membership)>1000]
#   subnetwork = igraph::groups(clu)[table(clu$membership)>10&table(clu$membership)<300]
#   
#   sub_nodelist = HMDB_clean2[subnetwork[[1]],]
#   g_subnetwork = make_ego_graph(subnetwork[[1]][1], graph = g, order = diameter(g), mode = "all")[[1]]
#   
#   
#   
#   temp_vertex_color = igraph::as_data_frame(g, "vertices") %>%
#     mutate(membership = clu$membership) %>%
#     mutate(color = case_when(
#       clu$csize[membership] == max(clu$csize) ~ my_palette[1],
#       clu$csize[membership] != 1 ~ my_palette[2],
#       clu$csize[membership] == 1 ~ my_palette[5]
#     )) %>%
#     pull(color)
#   
#   # pdf("full_all_connect_graph2.pdf")
#   # plot.igraph(g,
#   #      vertex.color = temp_vertex_color,
#   #      vertex.label = "",
#   #      vertex.frame.color = "#666666",
#   #      # vertex.label.color = "black",
#   #      # vertex.label.cex = 1,
#   #      edge.color = 'grey',
#   #      # edge.label = edge.attributes(g_subnetwork)$linktype,
#   #      vertex.size = sqrt(degree(g, mode = c("all")))+0.5,
#   #      # vertex.size = 2,
#   #      # edge.arrow.size = .05,
#   #      layout = layout_with_graphopt(g)
#   #      # main = paste("Subnetwork of node", interested_node)
#   # )
#   # dev.off()
# }

## Figure 1C ####
{
  ## Get results of
  # {
  #   # full rules
  #   g_edge_size
  #   # basic 15 rules
  #   g_basic_edge_size
  #   ## DQ 1910 metabolites
  #   ### biotransform_all 77 57 1776
  #   ### biotransform_basic15 137  75 1698
  #   ## Full 10304 metabolites 
  #   ### biotransform_all 318  324 9662
  #   ### biotransform_basic15 615  934 8755
  #   temp_size = g_edge_size
  #   unconnected = as.numeric(temp_size["1"])
  #   main_network = as.numeric(names(temp_size)[length(temp_size)])
  #   rest_network = nrow(merge_node_list) -unconnected-main_network
  #   print(c(unconnected, rest_network, main_network))
  # }
  
  
  ## ggplots
  connection_summary = data.frame(`Detected` = c(1698, 75, 137), 
                                  `All` = c(8755, 934, 615), 
                                  category = c("Main network", "Subnetworks", "Unconnected")) %>%
    gather(key = "metabolites", value = "number", -category) %>%
    mutate(metabolites = gsub("\\.", " ", metabolites))
  
  
  # pdf("bar_HMDB_connectivity.pdf", 
  #     width = 4,
  #     height = 3)
  # ggplot(connection_summary, aes(y = number, x = reorder(metabolites, -number), fill = forcats::fct_rev(category), label = number)) + 
  #   geom_bar(stat = "identity", 
  #            position = "fill" # make percentage graph
  #            ) +
  #   geom_text(size = 3,position = position_fill(vjust = 0.5)) +
  #   labs(# title = "Connectivity of HMDB metabolite formulas",
  #        x = NULL,
  #        y = "Fraction") + 
  #   guides(fill = guide_legend(
  #     title = "Metabolite status",
  #     reverse = F
  #     )) + 
  #   scale_y_continuous(expand = c(0,0),
  #                      labels = scales::percent,
  #                      breaks = scales::pretty_breaks(n = 8)
  #                      ) +
  #   scale_x_discrete(limits = c("Detected", "All")) +
  #   scale_fill_manual(values = c("Main network" = my_palette[1],
  #                     "Subnetworks" = my_palette[2],
  #                     "Unconnected" = my_palette[5])) +
  #   theme_classic(base_size = 12 # edit font size for all non-data text
  #                 ) +
  #   theme(plot.title = element_text(hjust = 0.5))
  # dev.off()
}



## Figure 1D-E - Plot step network ####
{
  # Plot step network
  {
    g = g_main
    
    node_list = main_node_list
    # node_list = unknown_formula
    edge_list = main_edge_list
    
    step_info = as.data.frame(table(node_list$steps), stringsAsFactors = F) %>%
      transmute(Step = Var1, No.metabolites = Freq) %>%
      mutate(Step = as.numeric(Step)) %>%
      mutate(Total.metabolites = cumsum(No.metabolites))  
    step_info2 = step_info %>%
      filter(!Step>9) %>%
      add_row(Step = ">9", No.metabolites=sum(step_info$No.metabolites[step_info$Step>9]), Total.metabolites = max(step_info$Total.metabolites))
    
    library(RColorBrewer)
    display.brewer.all(n=30, exact.n=FALSE)
    # color_set1 = c(brewer.pal(9, "Set1"),rep("#999999", 50))
    # color_set3_3color = rep(brewer.pal(3, "Set3"), c(3,3,50))
    # color_set1_3color = rep(c(brewer.pal(3, "Set1"),"#999999"), c(1,1,1,50))
    # color_set2_3color = rep(c(brewer.pal(3, "Set2"),"#999999"), c(1,2,3,50))
    
    color_spectral = c(brewer.pal(9, "Spectral"), rep("#666666", 50))
    
    vertex_color = color_spectral[main_node_list$steps]
    
    pdf("step_graph_equalsize.pdf")
    # for(i in 1:1){
      plot.igraph(g_main,
                  vertex.color = vertex_color,
                  vertex.label = "",
                  vertex.frame.color = "black",
                  # vertex.label.color = "black",
                  # vertex.label.cex = 1,
                  edge.color = 'grey',
                  # edge.label = edge.attributes(g_subnetwork)$linktype,
                  # vertex.size = sqrt(degree(g_main, mode = c("all")))+0.5,
                  vertex.size = 2,
                  # edge.arrow.size = .05,
                  edge.size = 0.001,
                  layout = layout_with_graphopt(g_main)
                  # main = paste("Subnetwork of node", interested_node)
      )
    # }
    dev.off()
    pdf("step_bar_graph.pdf", 
        width = 4,
        height = 3)
    step_info2$Step = as.factor(step_info2$Step)
    step_info2$Step  = factor(step_info2$Step, levels=c(levels(step_info2$Step)[-1],levels(step_info2$Step)[1]))
    ggplot(step_info2, aes(x = Step)) +
      geom_bar(aes(y = No.metabolites), stat = 'identity', fill = color_spectral[1:nrow(step_info2)]) + 
      # scale_x_discrete(limits = temp_position) +
      geom_point(aes(y = Total.metabolites/4)) + 
      geom_line(aes(y = Total.metabolites/4), group=1) +
      geom_text(aes(y = Total.metabolites/4, label = Total.metabolites), vjust = -.5) +
      scale_y_continuous(# sec.axis = sec_axis(~.*4, name = "Total.metabolites"), 
        limits = c(0,500), expand = c(0,0)) +
      # theme(legend.title = element_blank()) + 
      theme(legend.position = "right") +
      theme_classic()
    dev.off()
  }
}

## Figure 2B-C - HMDB connection and formulas bar graph summary  ####
{
  
  formula_optimized_correct = result_summary[[2]]$formula[1]
  formula_optimized_wrong = result_summary[[2]]$formula[2] - formula_optimized_correct
  formula_propagated_correct = result_summary[[2]]$formula[3]
  formula_propagated_wrong = result_summary[[2]]$formula[4] - formula_propagated_all
  
  connection_optimized_correct = result_summary[[2]]$connection[1]
  connection_optimized_wrong = result_summary[[2]]$connection[2] - connection_optimized_correct
  connection_propagated_corrects = result_summary[[2]]$connection[3]
  connection_propagated_wrong = result_summary[[2]]$connection[4] - connection_propagated_all
    
  
  formula_summary = data.frame(`Optimized` = c(formula_optimized_correct, formula_optimized_wrong),
                               `Pre-optimized` =  c(formula_propagated_correct, formula_propagated_wrong),
                               category = c("Correct", "Incorrect")) %>%
    gather(key = "cohorts", value = "number", -category) %>%
    mutate(cohorts = gsub("\\.", "-", cohorts))
  connection_summary = data.frame(`Optimized` = c(connection_optimized_correct, connection_optimized_wrong),
                               `Pre-optimized` =  c(connection_propagated_corrects, connection_propagated_wrong),
                               category = c("Correct", "Incorrect")) %>%
    gather(key = "cohorts", value = "number", -category) %>%
    mutate(cohorts = gsub("\\.", "-", cohorts))

  # pdf("bar_HMDB_formulas.pdf",
  #     width = 4,
  #     height = 4)
  # dev.new(width = 4, height = 3, unit = "in")
  figure_2B = ggplot(connection_summary, aes(y = number, x = cohorts, fill = forcats::fct_rev(category), label = number)) +
    geom_bar(stat = "identity",
             position = "fill" # make percentage graph
             ) +
    geom_text_repel(size = 3.5,position = position_fill(vjust = 0.5), max.iter=1,arrow=T) +
    labs(x = NULL,
         title = "Connection assignment",
         y = "Percentage") +
    guides(fill = guide_legend(
      title = NULL,
      reverse = F
      )) +
    scale_y_continuous(expand = c(0,0),
                       labels = scales::percent,
                       breaks = scales::pretty_breaks(n = 8)
                       ) +
    expand_limits(y = 1.05) +
    scale_fill_manual(values = c("Correct" = my_palette[1],
                                 "Incorrect" = my_palette[5])) +
    scale_x_discrete(limits = c("Pre-optimized", "Optimized")) +
    theme_classic(base_size = 14 # edit font size for all non-data text
                  ) +
    theme(plot.title = element_text(size = 14, hjust = 0.5),
          plot.margin = margin(0.5,0.5,0.5,0.5,"cm"),
          axis.text.x = element_text(angle = 30, hjust = .5, vjust = .5))
  
  figure_2C = ggplot(formula_summary, aes(y = number, x = cohorts, fill = forcats::fct_rev(category), label = number)) +
    geom_bar(stat = "identity",
             position = "fill" # make percentage graph
    ) +
    geom_text_repel(size = 3.5,position = position_fill(vjust = 0.5), max.iter=1,arrow=T) +
    labs(x = NULL,
         title = "Formula assignment",
         y = "Percentage") +
    guides(fill = guide_legend(
      title = NULL,
      reverse = F
    )) +
    scale_y_continuous(expand = c(0,0),
                       labels = scales::percent,
                       breaks = scales::pretty_breaks(n = 8)
    ) +
    expand_limits(y = 1.05) +
    scale_fill_manual(values = c("Correct" = my_palette[1],
                                 "Incorrect" = my_palette[5])) +
    scale_x_discrete(limits = c("Pre-optimized", "Optimized")) +
    theme_classic(base_size = 14 # edit font size for all non-data text
    ) +
    theme(plot.title = element_text(size = 14, hjust = 0.5),
          plot.margin = margin(0.5,0.5,0.5,0.5,"cm"),
          axis.text.x = element_text(angle = 30, hjust = .5, vjust = .5))
  # dev.off()
}

## Figure 2D - Heuristic bar graph in 1 ppm error ####
{
  ppm1_top1 = result_summary[[2]]$brute_force[1]
  ppm1_top3 = result_summary[[2]]$brute_force[2]
  total_assignment =result_summary[[2]]$formula[2]
  
  heursitic_bar_summary = data.frame(cohorts = c("Top1", "Top3"),
                                     number = c(ppm1_top1, ppm1_top3))
  
  # pdf("bar_heuristic_top1+3.pdf",
  #     width = 3,
  #     height = 4)
  figure_2D = ggplot(heursitic_bar_summary, aes(y = number/total_assignment, x = cohorts, fill = cohorts, label = number)) +
    geom_bar(stat = "identity",
             position = "identity" 
    ) +
    geom_text_repel(size = 3.5,position = position_stack(vjust = 0.5), max.iter=1,arrow=T) +
    labs(x = NULL,
         title = "Heuristic formula assignment",
         y = "Percentage") +
    guides(fill = guide_legend(
      title = NULL,
      reverse = F
    )) +
    scale_y_continuous(expand = c(0,0),
                       labels = scales::percent,
                       breaks = scales::pretty_breaks(n = 8)
    ) +
    expand_limits(y = 1.05) +
    scale_fill_manual(values = c(my_palette[5], my_palette[3])) +
    # scale_x_discrete(limits = c("Pre-optimized", "Optimized")) +
    theme_classic(base_size = 14 # edit font size for all non-data text
    ) +
    theme(plot.title = element_text(size = 14, hjust = 0.5),
          plot.margin = margin(0.5,0.5,0.5,0.5,"cm"),
          legend.position = "none"
          # axis.text.x = element_text(angle = 30, hjust = .5, vjust = .5)
          )
  # dev.off()
  
  
  
}


## Figure 2E - ppm_error formula assignment summary ####
{
  ppm_errors = sapply(result_summary, function(x) return(x$sigma[1]) )
  optimized_correct = sapply(result_summary, function(x) return(x$formula[1]) )
  optimized_wrong = sapply(result_summary, function(x) return(x$formula[2])) - optimized_correct
  
  ppm_error_assignment_summary = data.frame(Correct = optimized_correct,
                               Incorrect =  optimized_wrong,
                               ppm_errors = ppm_errors) %>%
    gather(key = "cohorts", value = "number", -ppm_errors) %>%
    mutate(cohorts = gsub("\\.", "-", cohorts))

  
  # pdf("bar_ppm_error_assignment.pdf",
  #     width = 5,
  #     height = 4)
  # dev.new(width = 4, height = 3, unit = "in")
  figure_2E = ggplot(ppm_error_assignment_summary, aes(y = number, x = factor(ppm_errors), fill = forcats::fct_rev(cohorts), label = number)) +
    geom_bar(stat = "identity",
             position = "stack" # make percentage graph
    ) +
    geom_text_repel(size = 3.5,position = position_stack(vjust = 0.5), max.iter=1,arrow=T) +
    labs(x = "Gaussian noise level (ppm)",
         title = "Formula assignment accuracy",
         y = "# of formula assignment") +
    guides(fill = guide_legend(
      title = NULL,
      reverse = F
    )) +
    scale_y_continuous(expand = c(0,0),
                       # labels = scales::percent,
                       breaks = scales::pretty_breaks(n = 8)
    ) +
    expand_limits(y = 1800) +
    scale_fill_manual(values = c("Correct" = my_palette[1],
                                 "Incorrect" = my_palette[5])) +
    theme_classic(base_size = 14 # edit font size for all non-data text
    ) +
    theme(plot.title = element_text(hjust = 0.5),
          plot.margin = margin(0.5,0.5,0.5,0.5,"cm")
          )
  # print(figure_2e)
  # dev.off()
  
}


## Figure 2F - ppm_error summary accuracy ####
{
  ppm_errors = sapply(result_summary, function(x) return(x$sigma[1]) )
  optimized_correct = sapply(result_summary, function(x) return(x$formula[1]) )
  optimized_all = sapply(result_summary, function(x) return(x$formula[2]))
  top1 = sapply(result_summary, function(x) return(x$brute_force[1]) )
  top3 = sapply(result_summary, function(x) return(x$brute_force[2]) )
  
  HMDB_ppm_summary = data.frame(ppm_errors = ppm_errors,
                                NetID = optimized_correct/optimized_all,
                                `Heuristic top1` = top1 / optimized_all,
                                `Heuristic top3` = top3 / optimized_all) %>%
    gather(key = "cohorts", value = "number", -ppm_errors) %>%
    mutate(cohorts = gsub("\\.", " ", cohorts)) %>%
    mutate(color = case_when(
      cohorts == "NetID" ~ "red",
      cohorts == "Heuristic top1" ~ "black",
      cohorts == "Heuristic top3" ~ "grey"
    ))
  
  # num_x = 4
  
  # pdf("bar_error_summary.pdf",
  #     width = 4,
  #     height = 4)
  # dev.new(width = 1, height = 1, unit = "in")
  figure_2F = ggplot(HMDB_ppm_summary, aes(y = number, x = factor(ppm_errors), group = forcats::fct_rev(cohorts), color = cohorts)) +
    geom_line(stat = "identity",
              size = 1.5
              # linetype = rep(c("solid", rep("dashed", 2)),2),
              # color = rep(c("red", "black", "grey"), num_x)
    ) + 
    geom_point(shape = 16,
               size = 4) + 
    labs(x = "Gaussian noise level (ppm)",
         title = "Formula assignment accuracy",
         y = "Accuracy percentage") +
    guides(color = guide_legend(
      title = NULL,
      reverse = T
    )) +
    scale_y_continuous(expand = c(0,0),
                       labels = scales::percent,
                       limits = c(0,1.05),
                       breaks = scales::pretty_breaks(n = 8)
    ) +
    # scale_x_discrete(limits = c("0.5", "1")) + 
    # scale_colour_manual(values=c("red", "black", "grey")) +
    scale_colour_manual(values=c("NetID" = my_palette[1],"Heuristic top1" = my_palette[5], "Heuristic top3" = my_palette[3])) +
    theme_classic(base_size = 14 # edit font size for all non-data text
    ) +
    theme(plot.title = element_text(hjust = 0.5),
          plot.margin = margin(0.5,0.5,0.5,0.5,"cm"))
  
  # dev.off()
  
}









## Figure 3C - Yeast dataset assignment summary ####
{
  yeast_formula_summary = data.frame(`Metabolites` = c(227, 1),
                               `Artifacts` =  c(293, 19),
                               category = c("Correct", "Incorrect")) %>%
    gather(key = "cohorts", value = "number", -category) %>%
    mutate(cohorts = gsub("\\.", "-", cohorts))
  
  pdf("yeast_formulas.pdf",
      width = 5,
      height = 1.5)
  # dev.new(width = 4, height = 3, unit = "in")
  figure_3C = ggplot(yeast_formula_summary, aes(y = number, x = cohorts, fill = forcats::fct_rev(category), label = number)) +
    geom_bar(stat = "identity",
             position = "fill" # make percentage graph
    ) +
    coord_flip() +
    geom_text_repel(size = 3.5,position = position_fill(vjust = 0.5), max.iter=1,arrow=T) +
    labs(x = NULL,
         title = "Formula assignment",
         y = "Percentage") +
    guides(fill = guide_legend(
      title = NULL,
      reverse = F
    )) +
    scale_y_continuous(expand = c(0,0),
                       labels = scales::percent,
                       breaks = scales::pretty_breaks(n = 5)
    ) +
    expand_limits(y = 1.05) +
    scale_fill_manual(values = c("Correct" = my_palette[1],
                                 "Incorrect" = my_palette[5])) +
    # scale_x_discrete(limits = c("Metabolites", "Artifacts")) +
    theme_classic(base_size = 14 # edit font size for all non-data text
    ) +
    theme(plot.title = element_text(size = 14, hjust = 0.5),
          plot.margin = margin(0.5,0.5,0.5,0.5,"cm")
          # axis.text.x = element_text(angle = 0, hjust = .5, vjust = .5)
          )
  print(figure_3C)
  dev.off()
  
}
## Merge graphs ####
{

  ggarrange(
    figure_2B, figure_2C, figure_2D,
    common.legend = T, legend = "right",
    # labels = c("B", "C","D"),
    nrow = 1, ncol = 3,
    align = "hv"
  ) %>%
    ggexport(filename = "figure_2b-d.pdf", width = 10, height = 4)
  
  ggarrange(
    figure_2E,
    figure_2F,
    align = "hv",
    nrow = 1, ncol = 2,
    widths = c(2,2)
  ) %>%
    ggexport(filename = "figure_2e-f.pdf", width = 10, height = 4)
  
}

save(result_summary, file = "result_summary.RData")
