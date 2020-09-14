# Library ####
{
  library(ggplot2)
  library(ggrepel)
  library(ggpubr)
  library(igraph)
  library(RColorBrewer)
  library(dplyr)
  library(tidyr)
  library(scales)
  library(janitor)
  library(readxl)
  setwd(dirname(rstudioapi::getSourceEditorContext()$path))
  rel_path = "C:/Users/Li Chen/Desktop/Network paper/network_figure_R"
  source('R_shiny_functions.R')
}

# Color setup ####
{
  library(colorscience)
  ## Functions ####
  hslToRgb = function (h, s, l){
    if(s == 0){
      r = g = b = l # achromatic
    }else{
      hue2rgb = function (p, q, t){
        if(t < 0) t = t+1;
        if(t > 1) t = t-1;
        if(t < 1/6) return(p + (q - p) * 6 * t)
        if(t < 1/2) return(q)
        if(t < 2/3) return(p + (q - p) * (2/3 - t) * 6)
        return (p)
      }
      
      q = ifelse(l < 0.5, l * (1 + s), l + s - l * s)
      p = 2 * l - q;
      r = hue2rgb(p, q, h + 1/3);
      g = hue2rgb(p, q, h);
      b = hue2rgb(p, q, h - 1/3);
    }
    
    return(c(round(r * 255), round(g * 255), round(b * 255)))
  }
  rgb2dex = function(rgb){
    paste0("#", paste0(as.hexmode(rgb),collapse = ""))
  }
  change_HSL = function(dex_color, H_change=0, S_change=0, L_change=0){
    HSL = RGB2HSL(t(col2rgb(dex_color)))
    HSL_new = HSL
    if(H_change >= 0){
      HSL_new[1] = HSL[1] + (1-HSL[1]) * H_change
    } else {
      HSL_new[1] = HSL[1] + HSL[1] * H_change
    }
    if(S_change >= 0){
      HSL_new[2] = HSL[2] + (1-HSL[2]) * S_change
    } else {
      HSL_new[2] = HSL[2] + HSL[2] * S_change
    }
    if(L_change >= 0){
      HSL_new[3] = HSL[3] + (1-HSL[3]) * L_change
    } else {
      HSL_new[3] = HSL[3] + HSL[3] * L_change
    }
    
    RGB_new = hslToRgb(HSL_new[1],HSL_new[2],HSL_new[3])
    
    return(rgb2dex(RGB_new))
  }
  
  display.brewer.all()
  # display.brewer.pal(8, "Set3")
  # display.brewer.pal(8, "Set1")
  my_palette = c(brewer.pal(4, "Set3"), rep("#666666", 50))
  my_palette_full = brewer.pal(12, "Set3")
  colfunc = colorRampPalette(c(my_palette[3], "white"))
  colfunc(5)[1]
  my_palette2 = brewer.pal(8, "Accent")[c(8,1,3,5,7,2,4,6)]
  my_palette3 = c(brewer.pal(8, "RdYlBu")[c(1)], my_palette2[1])
  
  test = c(my_palette[1], change_HSL(my_palette[1], S_change = 1) )
  # test = my_palette3
  image(1:length(test), 1, as.matrix(1:length(test)), 
        col=test,
        xlab="", ylab = "", xaxt = "n", yaxt = "n", bty = "n")
  
}


# General functions ####
## scientific plot ####
lb=1000 # min abs value to be included in the plot
dblog_trans <- function(){
  trans_new(name='dblog', transform = function(x) (log10(abs(x)+lb)-log10(lb))*sign(x),
            inverse = function(x) sign(x)*(10^(abs(x)+log10(lb))-lb))
}
fancy_scientific <- function(l) {
  # turn in to character string in scientific notation
  l <- format(l, scientific = TRUE)
  # quote the part before the exponent to keep all the digits
  l <- gsub("^(.*)e", "e", l)
  # turn the 'e+' into plotmath format
  l <- gsub("e\\+", "10^", l)
  # return this as an expression
  parse(text=l)
}
## plot_1 - thiamine abundance across tissues ####
{
  setwd(rel_path)
  setwd("thiamine abundance") # Set folder name
  
  header = read_xlsx("plot_data.xlsx", sheet = "Sheet1", n_max = 1, col_names = F)
  cohort = read_xlsx("plot_data.xlsx", sheet = "Sheet1", n_max = 1) 
  raw_data = read_xlsx("plot_data.xlsx", sheet = "Sheet1", skip = 2, col_names = F)
  colnames(raw_data) = header
  
  dt_cohort = cohort %>%
    select(-1) %>%
    gather(key = "header", value = "cohort")
  
  dt_number = raw_data %>%
    gather(key = "header", value = "number", -1)
  
  dt_plot = merge(dt_cohort, dt_number) %>%
    dplyr::rename(category = formula) %>%
    filter(complete.cases(.)) %>%
    filter(category != "Thiamine+[O]")
  
  dt_summary_plot = dt_plot %>%
    group_by(category, cohort) %>%
    summarise(
      sd = sd(number),
      number = mean(number),
      n = n(),
      se = sd/sqrt(n)
    ) %>%
    mutate(upper = round(number + se, 0), 
           lower = round(number - se, 0))
  
  plot_1A = ggplot(dt_summary_plot, aes(y = number, x = cohort, fill = category)) +
  # plot_1A = ggplot(dt_summary_plot, aes(y = number, x =category , fill =cohort )) +
    geom_bar(stat = "identity",
             width = .8,
             position = position_dodge(), 
             colour = "#333333"
    ) +
    geom_errorbar(aes(ymin=lower, ymax=upper),
                  width = 0.2,
                  # size = 2,
                  # width = 0.5, 
                  position=position_dodge(.8)) +
    geom_point(data = dt_plot,
               position = position_jitterdodge(jitter.width = 0.1,
                                               dodge.width = 0.8),
               # alpha = 0.5,
               shape = 1,
               show.legend=FALSE
    ) + 
    labs(x = NULL,
         # title = "U13C-glucose",
         y = "TIC") +
    guides(fill = guide_legend(
      title = NULL,
      reverse = F
    )
    ) +
    scale_y_continuous(expand = c(0,0),
                       trans = 'dblog',
                       limit=c(0,2e7), 
                       breaks = c(0,10^4,10^5,10^6,10^7,10^8,10^9), 
                       labels = fancy_scientific) + 
    scale_fill_manual(values = c(my_palette2)) +
    # scale_color_manual(values = "grey") + 
    # scale_x_discrete(labels = c("Brain", "Kidney", "Liver", "Pancreas", "Urine")) +
    # facet_wrap(~cohort) +
    theme_classic(base_size = 12 # edit font size for all non-data text
    ) +
    theme(plot.title = element_text(size = 12, hjust = 0.5),
          plot.margin = margin(0.1,0.1,0.1,0.1,"inch"),
          axis.text.x = element_text(angle = 0, hjust = .5, vjust = .5)
          # axis.ticks.x = element_blank()
    )
  print(plot_1A)
  
  pdf("tissue_thiamine_summary.pdf",
      width = 5,
      height = 1.45)
  print(plot_1A)
  dev.off()
  
}

## plot_2 - isotope labeling for thiamine ####
{
  setwd(rel_path)
  datapath = "fig4/"
  
  ## Data ####
  {
    filename = "U13C_glucose + no thiamine_corrected_2.xlsx"
    thiamine_labeling1 = read_xlsx(paste0(datapath, filename), sheet = "summary") %>%
      select_if(~sum(!is.na(.)) > 0)
    colnames(thiamine_labeling1)[1:2] = c("category", "cohort")
    thiamine_labeling1 = thiamine_labeling1 %>%
      gather(key = "key", value = "number", -category, -cohort) %>%
      mutate(cohort = ifelse(cohort == "Fully_labeled", "Fully\nlabeled", cohort))
    thiamine_labeling_summary1 = thiamine_labeling1 %>%
      group_by(category, cohort) %>%
      summarise(
        sd = sd(number),
        number = mean(number)
      )
    
    datapath = "fig4/"
    filename = "13C + unlabeled thiamine_corrected_2.xlsx"
    thiamine_labeling2 = read_xlsx(paste0(datapath, filename), sheet = "summary") %>%
      select_if(~sum(!is.na(.)) > 0)
    colnames(thiamine_labeling2)[1:2] = c("category", "cohort")
    thiamine_labeling2 = thiamine_labeling2 %>%
      gather(key = "key", value = "number", -category, -cohort) %>%
      mutate(cohort = ifelse(cohort == "Fully_labeled", "Fully\nlabeled", cohort))
    thiamine_labeling_summary2 = thiamine_labeling2 %>%
      group_by(category, cohort) %>%
      summarise(
        sd = sd(number),
        number = mean(number)
      )
  }
  
  
  ## Plot ####
  {
    plot_2A = ggplot(thiamine_labeling_summary1, aes(y = number, x = cohort, fill = category)) +
      geom_bar(stat = "identity",
               width = .8,
               position = position_dodge(), 
               colour = "#333333"
      ) +
      geom_errorbar(aes(ymin=number-sd, ymax=number+sd),
                    width = 0.2,
                    # size = 2,
                    # width = 0.5, 
                    position=position_dodge(.8)) +
      geom_point(data = thiamine_labeling1,
                 position = position_jitterdodge(jitter.width = 0.1,
                                                 dodge.width = 0.8),
                 # alpha = 0.5,
                 shape = 1,
                 show.legend=FALSE
      ) + 
      labs(x = NULL,
           # title = "U13C-glucose",
           y = "Labeling fraction") +
      guides(fill = guide_legend(
        title = NULL,
        reverse = F
      )
      ) +
      scale_y_continuous(expand = c(0,0),
                         breaks = scales::pretty_breaks(n = 2)
      ) +
      expand_limits(y=c(0, 1.05)) +
      scale_x_discrete(limits = c("M0", "M2", "Fully\nlabeled")) +
      scale_fill_manual(values = c(my_palette2),
                        labels = c("Thiamine", 
                                   # "Thiamine+[O]", 
                                   "Thiamine+[C2H2O]", 
                                   "Thiamine+[C2H4O]")) + 
      scale_color_manual(values = "grey") + 
      # facet_wrap(~cohort) +
      theme_classic(base_size = 12 # edit font size for all non-data text
      ) +
      theme(plot.title = element_text(size = 12, hjust = 0.5),
            plot.margin = margin(0.5,0.5,0.5,0.5,"cm"),
            axis.text.x = element_text(angle = 0, hjust = .5, vjust = .7)
            # axis.ticks.x = element_blank()
      )
    print(plot_2A)
    
    plot_2B = ggplot(thiamine_labeling_summary2, aes(y = number, x = cohort, fill = category)) +
      geom_bar(stat = "identity",
               width = .8,
               position = position_dodge(), 
               colour = "#333333"
      ) +
      geom_errorbar(aes(ymin=number-sd, ymax=number+sd),
                    width = 0.2,
                    # size = 2,
                    # width = 0.5, 
                    position=position_dodge(.8)) +
      geom_point(data = thiamine_labeling2,
                 position = position_jitterdodge(jitter.width = 0.1,
                                                 dodge.width = 0.8),
                 shape = 1,
                 show.legend=FALSE) + 
      labs(x = NULL,
           # title = "U13C-glucose + 12C-thiamine",
           y = "Labeling fraction") +
      guides(fill = guide_legend(
        title = NULL,
        reverse = F
      )
      ) +
      scale_y_continuous(expand = c(0,0),
                         breaks = scales::pretty_breaks(n = 2)
      ) +
      expand_limits(y=c(0, 1.05)) +
      scale_x_discrete(limits = c("M0", "M2", "Fully\nlabeled")) +
      scale_fill_manual(values = c(my_palette2),
                        labels = c("Thiamine", 
                                   # "Thiamine+[O]", 
                                   "Thiamine+[C2H2O]", 
                                   "Thiamine+[C2H4O]")) + 
      scale_color_manual(values = "grey") + 
      # facet_wrap(~cohort) +
      theme_classic(base_size = 12 # edit font size for all non-data text
      ) +
      theme(plot.title = element_text(size = 12, hjust = 0.5),
            plot.margin = margin(0.5,0.5,0.5,0.5,"cm"),
            axis.text.x = element_text(angle = 0, hjust = .5, vjust = .7)
            # axis.ticks.x = element_blank()
      )
    print(plot_2B)
  }
  
  
  ## Output ####
  {
    ggpubr::ggarrange(
      plot_2A,
      plot_2B,
      common.legend = T, legend = "none",
      align = "hv",
      nrow = 1, ncol = 2
    ) %>%
      ggexport(filename = "figure_3C_new.pdf", width = 5.5, height = 1.6)
    
  }
  
}


## plot_2_sup1 - more isotope labeling for thiamine ####
{
  setwd(rel_path)
  setwd("plot_2_sup1 - more isotope for thiamine/")
  ## Data ####
  {
    # Spike in
    filename = "12C_thiamine_spike_in.xlsx"
    raw_data = read_xlsx(filename, sheet = "summary") %>%
      select_if(~sum(!is.na(.)) > 0)
    colnames(raw_data)[1:2] = c("category", "cohort")
    thiamine_labeling = raw_data %>%
      gather(key = "key", value = "number", -category, -cohort) %>%
      mutate(cohort = ifelse(cohort == "Fully_labeled", "Fully\nlabeled", cohort))
    
    thiamine_labeling = thiamine_labeling %>%
      filter(!grepl("C10H9", category))
    # U13C
    filename = "WL_thiamine+C4H6O3.xlsx"
    raw_data = read_xlsx(filename, sheet = "summary") %>%
      select_if(~sum(!is.na(.)) > 0)
    colnames(raw_data)[1:2] = c("category", "cohort")
    thiamine_labeling2 = raw_data %>%
      gather(key = "key", value = "number", -category, -cohort) %>%
      mutate(cohort = ifelse(cohort == "Fully_labeled", "Fully\nlabeled", cohort))
    
    thiamine_labeling2 = thiamine_labeling2 %>%
      filter(!grepl("C10H9", category))
  }
  
  
  
  ## Plot ####
  thiamine_isotope_ggplot = function(thiamine_labeling_select){
    thiamine_labeling_normalize = thiamine_labeling_select %>%
      group_by(category, key) %>%
      summarise(sum = sum(number)) %>%
      merge(thiamine_labeling_select) %>%
      mutate(number = number / sum)
    
    thiamine_labeling_summary = thiamine_labeling_normalize %>%
      group_by(category, cohort) %>%
      summarise(
        sd = sd(number),
        number = mean(number)
      )
    
    plot_2A = ggplot(thiamine_labeling_summary, aes(y = number, x = cohort, fill = category)) +
      geom_bar(stat = "identity",
               width = .8,
               position = position_dodge(), 
               colour = "#333333"
      ) +
      geom_errorbar(aes(ymin=number-sd, ymax=number+sd),
                    width = 0.2,
                    # size = 2,
                    # width = 0.5, 
                    position=position_dodge(.8)) +
      geom_point(data = thiamine_labeling_normalize,
                 position = position_jitterdodge(jitter.width = 0.1,
                                                 dodge.width = 0.8),
                 # alpha = 0.5,
                 shape = 1,
                 show.legend=FALSE
      ) + 
      labs(x = NULL,
           # title = "U13C-glucose",
           y = "Labeling fraction") +
      guides(fill = guide_legend(
        title = NULL,
        reverse = F
      )
      ) +
      scale_y_continuous(expand = c(0,0),
                         breaks = scales::pretty_breaks(n = 2)
      ) +
      expand_limits(y=c(0, 1.05)) +
      scale_x_discrete(limits = c("M0", "M4", "Fully\nlabeled")) +
      scale_fill_manual(values = c(my_palette2)
                        # labels = c("Thiamine",
                        #            # "Thiamine+[O]",
                        #            "Thiamine+[C2H2O]",
                        #            "Thiamine+[C2H4O]")
      ) + 
      scale_color_manual(values = "grey") + 
      # facet_wrap(~cohort) +
      theme_classic(base_size = 12 # edit font size for all non-data text
      ) +
      theme(plot.title = element_text(size = 12, hjust = 0.5),
            plot.margin = margin(0.5,0.5,0.5,0.5,"cm"),
            axis.text.x = element_text(angle = 0, hjust = .5, vjust = .7)
            # axis.ticks.x = element_blank()
      )
    print(plot_2A)
    return(plot_2A)
  }
  
  {
    thiamine_labeling_select = thiamine_labeling %>%
      filter(grepl("12C", key)) %>%
      filter(grepl("C4H6O3", category))
    C4H6O3_12C = thiamine_isotope_ggplot(thiamine_labeling_select)
    
    thiamine_labeling_select = thiamine_labeling %>%
      filter(grepl("12C", key)) %>%
      filter(grepl("C4H8O", category))
    C4H8O_12C = thiamine_isotope_ggplot(thiamine_labeling_select)
    
    thiamine_labeling_select = thiamine_labeling %>%
      filter(grepl("13C", key)) %>%
      filter(grepl("C4H6O3", category))
    C4H6O3_13C_SPIKE = thiamine_isotope_ggplot(thiamine_labeling_select)
    
    thiamine_labeling_select = thiamine_labeling %>%
      filter(grepl("13C", key)) %>%
      filter(grepl("C4H8O", category))
    C4H8O_13C_SPIKE = thiamine_isotope_ggplot(thiamine_labeling_select)
    
    thiamine_labeling_select = thiamine_labeling2 %>%
      filter(grepl("13C", key))
    C4H6O3_13C = thiamine_isotope_ggplot(thiamine_labeling_select)
    
  }
  
  ## Output ####
  {
    
    ggpubr::ggarrange(
      C4H6O3_12C,
      C4H6O3_13C_SPIKE, 
      C4H6O3_13C,
      C4H8O_12C,
      C4H8O_13C_SPIKE,
      common.legend = T, legend = "none",
      align = "hv",
      nrow = 2, ncol = 3
    ) %>%
      ggexport(filename = "plot_2_sup1.pdf", width = 8, height = 3.2)
  }
  
  

}
## plot_3 - overview newtork plot for yeast neg data ####
{
  setwd(dirname(rstudioapi::getSourceEditorContext()$path))
  setwd("./Unknown in IO RT/Sc_neg")
  filename = "20200704085443_output.rds" # Sc_neg
  
  # Global parameter in the background ####
  {
    dt = readRDS(filename)
    
    ilp_nodes = dt$ilp_nodes %>%
      mutate(medMz = signif(medMz, 7),
             medRt = round(medRt, 2),
             log10_inten = round(log10_inten, 2),
             ppm_error = round(ppm_error, 2))
    
    ilp_nodes_ilp = ilp_nodes %>%
      filter(ilp_result == 1) %>%
      filter(class != "Unknown") %>%
      filter(T)
    
    ilp_edges = dt$ilp_edges
    
    g_met = initiate_g_met(ilp_nodes, ilp_edges)
    
    core_met = ilp_nodes %>%
      filter(steps == 0) %>%
      filter(class == "Metabolite")
    
    core_nonmet = ilp_nodes %>%
      filter(steps %% 1 == 0) %>%
      filter(class != "Unknown") 
    
    g_nonmet = initiate_g_nonmet(ilp_nodes, dt$ilp_edges, dt$heterodimer_ilp_edges)
    
    # ilp_edges_annotate_met = igraph::as_data_frame(g_met, "edges")
    # ilp_edges_annotate_nonmet = igraph::as_data_frame(g_nonmet, "edges")
    core_rank = core_annotate(ilp_nodes, dt$FormulaSet_df, dt$LibrarySet)
    
  }
  
  # Network palette ####
  {
    
    new_purple = colorRampPalette(c(my_palette[3], "#FFFFFF"))(10)
    new_green = colorRampPalette(c(my_palette[1], "#8DFFC7"))(10)
    new_yellow = colorRampPalette(c(my_palette[2], "#FFFF00"))(10)
    my_palette_network = c(change_HSL(new_green[6], S_change = 0.8), 
                           change_HSL(new_yellow[3], S_change = 0.8),
                           change_HSL(new_purple[2], S_change = 0.8),
                           my_palette[4])
    
    
    test = my_palette_network
    image(1:length(test), 1, as.matrix(1:length(test)),
          col=test,
          xlab="", ylab = "", xaxt = "n", yaxt = "n", bty = "n")
  }
  
  # Further filtering ####
  {
    g_met2_node = igraph::as_data_frame(g_met, "vertices") %>%
      filter(ilp_result > 0.01) %>%
      mutate(color = ifelse(class == "Metabolite", my_palette_network[1], my_palette_network[2]))
    
    g_met2_edges = igraph::as_data_frame(g_met, "edges") %>%
      filter(from %in% g_met2_node$name & to %in% g_met2_node$name) %>%
      mutate(color = my_palette[1])
    
    g_met2 = graph_from_data_frame(g_met2_edges,
                                   vertices = g_met2_node,
                                   directed = F)
    
    g_nonmet2_node = igraph::as_data_frame(g_nonmet, "vertices") %>%
      filter(ilp_result > 0.01) %>%
      mutate(color = case_when(
        class == "Metabolite" ~ my_palette_network[1],
        class == "Putative Metabolite" ~ my_palette_network[2],
        class == "Artifact" ~ my_palette_network[3],
        class == "Unknown" ~ "#666666"
      ))
    g_nonmet2_edges = igraph::as_data_frame(g_nonmet, "edges") %>%
      filter(from %in% g_nonmet2_node$name & to %in% g_nonmet2_node$name,
             ilp_result > 0.01) %>%
      mutate(color = new_purple[6])
    g_nonmet2 = graph_from_data_frame(g_nonmet2_edges,
                                      vertices = g_nonmet2_node,
                                      directed = F)
    
    g_all_nodes = bind_rows(g_met2_node, g_nonmet2_node) %>%
      distinct(name, .keep_all = T) %>%
      # filter(log10_inten >= 5) %>%
      filter(class != "Unknown")
    g_all_edges = bind_rows(g_met2_edges, g_nonmet2_edges) %>%
      distinct() %>%
      filter(from %in% g_all_nodes$name & to %in% g_all_nodes$name)
  }
  
  # Sub-graph ####
  {
    # Problem 1 - many single node in the output
    # Reason - Heterodimer edge score enables nodes without regular connections
    # Solutions - Remove nodes do not have regular connections
    
    g_all_nodes = g_all_nodes %>%
      filter(name %in% g_all_edges$from | name %in% g_all_edges$to)
    
    g_all = graph_from_data_frame(g_all_edges,
                                  directed = F,
                                  vertices = g_all_nodes)
    
    # Problem 2 - Putative metabolites in network without connecting to metabolites
    # Reason - Edge score by connection with isotope/adduct etc.
    # Solutions - Only retain subnetwork where step=0 annotation exist
    g_sub=g_all
    clu=components(g_sub)
    #subnetwork criteria 
    mainnetwork = igraph::groups(clu)[table(clu$membership) == max(table(clu$membership))]
    subnetwork = igraph::groups(clu)[table(clu$membership) < length(mainnetwork[[1]])]
    
    
    step0 = g_all_nodes %>%
      filter(steps == 0) %>%
      pull(name)
    subnetwork_valid = sapply(subnetwork, function(x){
      any(x %in% step0)
    })
    nodes_invalid = unlist(subnetwork[!subnetwork_valid])
    
    g_all_nodes_valid = g_all_nodes %>%
      filter(!name %in% nodes_invalid)
    g_all_edges_valid = g_all_edges %>%
      filter(from %in% g_all_nodes_valid$name & to %in% g_all_nodes_valid$name)
    g_all_valid = graph_from_data_frame(g_all_edges_valid,
                                  directed = F,
                                  vertices = g_all_nodes_valid)

    nodes_metabolite = g_all_nodes %>%
      filter(class == "Metabolite") %>%
      pull(name)
    
    subnetwork_valid_nonmet = sapply(subnetwork, function(x){
      (any(x %in% step0)) & (all(x %in% nodes_metabolite))
    })
    test = subnetwork[subnetwork_valid_nonmet]
    
    test2 = g_all_nodes_valid %>%
      filter(name %in% test$`22`)
    
    
    # g_main = make_ego_graph(graph=g_sub,
    #                         order=diameter(g_sub),
    #                         nodes = mainnetwork[[1]][1],
    #                         mode="all")[[1]]
    
  }

  
  # Plots ####
  {
    g_interest = g_all_valid %>% # g_all
      # set_vertex_attr("color", value = my_palette[1]) %>%
      set_vertex_attr("size", value = 2) %>%
      set_vertex_attr("label", value = NA) %>%
      set_vertex_attr("frame.color", value = "#AAAAAA") %>%
      # set_edge_attr("color", value = "#CCCCCC") %>%
      set_edge_attr("width", value = 1)
    
    setwd(dirname(rstudioapi::getSourceEditorContext()$path))
    
    pdf("network_plot_graphopt_test.pdf",
        width = 40,
        height = 40)
    
    # Plot various layout
    # {
    #   layouts <- grep("^layout_", ls("package:igraph"), value=TRUE)[-1] 
    #   # Remove layouts that do not apply to our graph.
    #   layouts <- layouts[!grepl("bipartite|merge|norm|sugiyama|tree", layouts)]
    #   
    #   for (layout in layouts) {
    #     print(layout)
    #     l <- do.call(layout, list(g_main)) 
    #     plot(g_main, edge.arrow.mode=0, layout=l, main=layout) }
    # }
    
    plot.igraph(g_interest,
                # layout = layout_with_fr # Faster but not beautiful
                # layout = layout_with_lgl
                layout = layout_with_graphopt
            
    )
    dev.off()
  }
}


## plot_4 - statistics of global networks ####
{
  setwd(dirname(rstudioapi::getSourceEditorContext()$path))
  datapath = "./Unknown in IO RT/Sc_neg"
  # datapath = "./Unknown in IO RT/Sc_pos"
  # datapath = "./WL_liver_neg/200711"
  # datapath = "./WL_liver_pos/200711"
  setwd(datapath)
  filename = "20200704085443_output.rds" # Sc_neg
  # filename = "20200704005808_output.rds" # Sc_pos
  # filename = "20200709234041_output.rds" # Liver_neg
  # filename = "20200711005313_output.rds" # Liver_pos
  
  filename_wl = "../Lu-Table-S4-final_Sc_neg.xlsx"
  # filename_wl = "../Lu-Table-S4-final_Sc_pos.xlsx"
  # filename_wl = "../pks_liver_neg_buff_2020-02-21.xlsx"
  # filename_wl = "../pks_liver_pos_buff_2020-02-21.xlsx"
  
  # Global parameter in the background ####
  {
    dt = readRDS(filename)
    
    ilp_nodes = dt$ilp_nodes %>%
      mutate(medMz = signif(medMz, 7),
             medRt = round(medRt, 2),
             log10_inten = round(log10_inten, 2),
             ppm_error = round(ppm_error, 2))
    
    ilp_edges = dt$ilp_edges
    
    g_met = initiate_g_met(ilp_nodes, dt$ilp_edges)
    
    core_met = ilp_nodes %>%
      filter(steps == 0) %>%
      filter(class == "Metabolite")
    
    core_nonmet = ilp_nodes %>%
      filter(steps %% 1 == 0) %>%
      filter(class != "Unknown")
    
    g_nonmet = initiate_g_nonmet(ilp_nodes, dt$ilp_edges, dt$heterodimer_ilp_edges)
    
    # ilp_edges_annotate_met = igraph::as_data_frame(g_met, "edges")
    # ilp_edges_annotate_nonmet = igraph::as_data_frame(g_nonmet, "edges")
    core_rank = core_annotate(ilp_nodes, dt$FormulaSet_df, dt$LibrarySet)
    
    data(isotopes)
    
  }
  
  
  # Node ####
  {
    
    
    WL = readxl::read_xlsx(filename_wl,
                           # sheet = "Sheet1",
                           guess_max = 1e6
    ) %>%
      dplyr::rename(medRt = rt,
                    medMz = mz) %>%
      dplyr::rename(Formula = Formula...32,
                    Feature = Feature...33,
                    Background = Background...25) %>%
      filter(T)
    
    
    test = ilp_nodes %>% 
      filter(ilp_result > 0.01)
    test2 = merge(WL, test, by.x = "Index", by.y = "Input_id", all.x = T, suffixes = c("",".y")) %>%
      dplyr::select(colnames(WL), node_id, class, path, formula, ppm_error, parent_id, parent_formula, transform, category, mass)
    
    test3 = test2 %>%
      mutate(Formula = ifelse(is.na(Formula), "", Formula)) %>%
      mutate(Formula = check_chemform(isotopes, Formula)$new_formula)
    
    test3_filter = test3 %>%
      filter(is.na(Background)) %>%
      filter(T) 
    
    tabyl(test3_filter, class, Background)
    # print(c(nrow(test3), nrow(test3)-nrow(test3_filter)))

    
    # Plot ####
    # {
    #   class_info = as.data.frame(table(test3_filter$class), stringsAsFactors = F) %>%
    #     transmute(Class = Var1, `Number of peaks` = Freq) 
    #   
    #   plot4_node = ggplot(class_info, aes(x = Class)) +
    #     geom_bar(aes(y = `Number of peaks`), stat = 'identity', fill = my_palette[1:4]) + 
    #     # scale_x_discrete(limits = temp_position) +
    #     # geom_point(aes(y = Total.metabolites/4)) + 
    #     # geom_line(aes(y = Total.metabolites/4), group=1) +
    #     # geom_text(aes(y = Total.metabolites/4, label = Total.metabolites), vjust = -.5) +
    #     scale_y_continuous(# sec.axis = sec_axis(~.*4, name = "Total.metabolites"), 
    #       limits = c(0,4000), expand = c(0,0)) +
    #     theme(legend.position = "right") + 
    #     theme_classic(base_size = 11 # edit font size for all non-data text
    #     ) +
    #     theme(plot.title = element_text(size = 11, hjust = 0.5),
    #           plot.margin = margin(0.5,0.5,0.5,0.5,"cm"),
    #           legend.position = "none"
    #           # axis.text.x = element_text(angle = 30, hjust = .5, vjust = .5)
    #     )
    #   
    # }
    
  }
  # Edge ####
  {
    ilp_edges_filter = ilp_edges %>%
      arrange(-ilp_result) %>%
      distinct(edge_id, .keep_all = T) %>%
      filter(ilp_nodes1 %in% test$ilp_node_id, 
             ilp_nodes2 %in% test$ilp_node_id) %>%
      filter(!(category != "Biotransform" & ilp_result == 0))
      
  
  }
  # Network ####
  
  # Output ####
  {
    print(c(nrow(test3_filter)))
    print(tabyl(test3_filter$class))
    print(nrow(ilp_edges_filter))
    print(nrow(ilp_edges_filter %>% filter(category == "Biotransform")))
    print(nrow(ilp_edges_filter %>% filter(category != "Biotransform")))
  }
  
}
## plot_5 - validate CN num
{
  setwd(dirname(rstudioapi::getSourceEditorContext()$path))
  datapath = "./Unknown in IO RT/Sc_neg"
  setwd(datapath)
  filename = "20200704085443_output.rds" # Sc_neg
  filename_wl = "../Lu-Table-S4-final_Sc_neg.xlsx"
  # Global parameter in the background ####
  {
    dt = readRDS(filename)
    
    ilp_nodes = dt$ilp_nodes %>%
      mutate(medMz = signif(medMz, 7),
             medRt = round(medRt, 2),
             log10_inten = round(log10_inten, 2),
             ppm_error = round(ppm_error, 2))
    
    ilp_edges = dt$ilp_edges
    
    g_met = initiate_g_met(ilp_nodes, dt$ilp_edges)
    
    core_met = ilp_nodes %>%
      filter(steps == 0) %>%
      filter(class == "Metabolite")
    
    core_nonmet = ilp_nodes %>%
      filter(steps %% 1 == 0) %>%
      filter(class != "Unknown")
    
    g_nonmet = initiate_g_nonmet(ilp_nodes, dt$ilp_edges, dt$heterodimer_ilp_edges)
    
    # ilp_edges_annotate_met = igraph::as_data_frame(g_met, "edges")
    # ilp_edges_annotate_nonmet = igraph::as_data_frame(g_nonmet, "edges")
    core_rank = core_annotate(ilp_nodes, dt$FormulaSet_df, dt$LibrarySet)
    
    data(isotopes)
    
  }
  
  
  # Node ####
  {
    WL = readxl::read_xlsx(filename_wl,
                           # sheet = "Sheet1",
                           guess_max = 1e6
    ) %>%
      dplyr::rename(medRt = rt,
                    medMz = mz) %>%
      dplyr::rename(Formula = Formula...32,
                    Feature = Feature...33,
                    Background = Background...25) %>%
      filter(T)
    
    
    test = ilp_nodes %>% 
      filter(ilp_result > 0.01)
    test2 = merge(WL, test, by.x = "Index", by.y = "Input_id", all.x = T, suffixes = c("",".y")) %>%
      dplyr::select(colnames(WL), node_id, class, path, formula, ppm_error, parent_id, parent_formula, transform, category, mass)
    
    test3 = test2 %>%
      mutate(Formula = ifelse(is.na(Formula), "", Formula)) %>%
      mutate(Formula = check_chemform(isotopes, Formula)$new_formula)
    
    test3_filter = test3 %>%
      filter(is.na(Background)) %>%
      filter(class != "Unknown") %>%
      filter(C != 0 | N != 0) %>%
      filter(T) 
  }
  
  ## Evaluatoin ####
  {
    test4 = test3_filter %>%
      filter(sig > log10(1e5)) %>%
      mutate(new_C = sapply(formula,elem_num_query, "C")) %>%
      mutate(new_N = sapply(formula,elem_num_query, "N")) %>%
      mutate(CN_match = new_C == C & new_N == N)
    tabyl(test4, class, CN_match)
    
    test4_filter = test4 %>%
      # filter(class == "Putative Metabolite") %>%
      filter(class == "Metabolite") %>%
      filter(!CN_match) %>%
      filter(T)
    
    
  }
}
## plot_5_1 - glucosyl-taurine chromatogram ####
{
  setwd(rel_path)
  setwd("plot_5 - glucosyltaurine")
  ## Data ####
  {
    # Spike in
    filename = "glucosyl-taurine-raw-data-0902.xlsx"
    raw_data = read_xlsx(filename, sheet = "chromatogram") %>%
      select_if(~sum(!is.na(.)) > 0)
  }
  
  ## Plot ####
  {
    liver_LC = ggplot(raw_data, aes(x = liver)) +
      # geom_point(aes(y = liver_TIC_normalized)) + 
      geom_line(aes(y = liver_TIC_normalized),
                color=my_palette3[1], size=1) +
      # geom_text(aes(y = Total.metabolites/4, label = Total.metabolites), vjust = -.5) +
      labs(x = "RT (min)",
           # title = "U13C-glucose",
           y = "Normalized ion counts") +
      scale_y_continuous(limits = c(0,100), expand = c(0,0),
                         breaks = scales::pretty_breaks(n = 2)
      ) +
      
      scale_x_continuous(limits = c(9.5,13.5)) + 
      # theme(legend.title = element_blank()) + 
      theme(legend.position = "right") +
      theme_classic(base_size = 12 # edit font size for all non-data text
      ) +
      theme(plot.title = element_text(size = 12, hjust = 0.5),
            plot.margin = margin(0.5,0.5,0.5,0.5,"cm"),
            axis.text.x = element_text(angle = 0, hjust = .5, vjust = .7)
            # axis.ticks.x = element_blank()
      )
    
    synthesis_LC = ggplot(raw_data, aes(x = synthesis)) +
      # geom_point(aes(y = liver_TIC_normalized)) + 
      geom_line(aes(y = synthesis_TIC_normalized),
                color=my_palette3[2], size=1) +
      # geom_text(aes(y = Total.metabolites/4, label = Total.metabolites), vjust = -.5) +
      labs(x = "RT (min)", 
           # title = "U13C-glucose",
           y = "Normalized ion counts") +
      scale_y_continuous(limits = c(0,100), expand = c(0,0),
                         breaks = scales::pretty_breaks(n = 2)
      ) +
      scale_x_continuous(limits = c(9.5,13.5)) + 
      # theme(legend.title = element_blank()) + 
      theme(legend.position = "right") +
      theme_classic(base_size = 12 # edit font size for all non-data text
      ) +
      theme(plot.title = element_text(size = 12, hjust = 0.5),
            plot.margin = margin(0.5,0.5,0.5,0.5,"cm"),
            axis.text.x = element_text(angle = 0, hjust = .5, vjust = .7)
            # axis.ticks.x = element_blank()
      )
  }
    
  ## Output ####
  {
    ggpubr::ggarrange(
      liver_LC,
      synthesis_LC,
      common.legend = T, legend = "none",
      align = "hv",
      nrow = 2, ncol = 1
    ) %>%
      ggexport(filename = "plot_5.pdf", width = 3.5, height = 3.2)
  }
  
  
  
  
}
## plot_5_2 - glucosyl-taurine MS2 ####
{
  setwd(rel_path)
  setwd("plot_5 - glucosyltaurine")
  
  ## Function ####
  {
    my_centroid = function(s, topn = 10, ms_dif_ppm=20, ms2report_cutoff=0.001){
      
      ms_dif_ppm=ms_dif_ppm/10^6
      ##Group MS groups
      s = s %>%
        arrange(mz)
      mzs = s$mz
      count = 1
      MZ_group = rep(1,(length(mzs)))
      for(i in 2:length(mzs)){
        if(mzs[i]-mzs[i-1]>mzs[i-1]*ms_dif_ppm){
          count = count+1
        }
        MZ_group[i]=count
      }
      s["MZ_group"]=MZ_group
      
      
      s2 = s %>%
        group_by(MZ_group) %>%
        summarise(Centroid_intensity = max(intensity),
                  Centroid_mz = sum(mz*intensity)/sum(intensity)) %>%
        filter(complete.cases(.))
      
      s3 = s2 %>% 
        filter(Centroid_intensity > ms2report_cutoff * max(Centroid_intensity)) %>%
        # select(c("Centroid_mz", "Centroid_intensity")) %>%
        distinct(.keep_all = T) %>%
        arrange(-Centroid_intensity) %>%
        slice(1:topn)
      
      return(s3)
    }
  }
  ## Data ####
  {
    filename = "glucosyl-taurine-raw-data-0902.xlsx"
    raw_data = read_xlsx(filename, sheet = "MS2-posi") %>%
      select_if(~sum(!is.na(.)) > 0)
    raw_data_neg = read_xlsx(filename, sheet = "MS2-nega") %>%
      select_if(~sum(!is.na(.)) > 0)
    
    ## Change the raw_data into raw_data_neg to get second graph
    liver_MS2 = raw_data %>%
      select(1,2) %>%
      filter(complete.cases(.))
    colnames(liver_MS2) = c("mz","intensity") 
    liver_MS2_normalized = my_centroid(liver_MS2, topn = 10) %>%
      mutate(normalized = Centroid_intensity / max(Centroid_intensity) * 100)
    
    synthesis_MS2 = raw_data %>%
      select(4,5) %>%
      filter(complete.cases(.))
    colnames(synthesis_MS2) = c("mz","intensity")
    synthesis_MS2_normalized = my_centroid(synthesis_MS2, topn = 10) %>%
      mutate(normalized = Centroid_intensity / max(Centroid_intensity) * 100)
    
    liver_MS2_neg = raw_data_neg %>%
      select(1,2) %>%
      filter(complete.cases(.))
    colnames(liver_MS2_neg) = c("mz","intensity") 
    liver_MS2_neg_normalized = my_centroid(liver_MS2_neg, topn = 10) %>%
      mutate(normalized = Centroid_intensity / max(Centroid_intensity) * 100)
    
    synthesis_MS2_neg = raw_data_neg %>%
      select(4,5) %>%
      filter(complete.cases(.))
    colnames(synthesis_MS2_neg) = c("mz","intensity")
    synthesis_MS2_neg_normalized = my_centroid(synthesis_MS2_neg, topn = 10) %>%
      mutate(normalized = Centroid_intensity / max(Centroid_intensity) * 100)
  }
  
  ## Plot ####
  {
    ms2_pos = ggplot(liver_MS2_normalized,aes(x=Centroid_mz, y= normalized)) +
      geom_linerange(aes(ymin = 0, ymax = normalized),
                     color=my_palette3[1]) +
      geom_linerange(data = synthesis_MS2_normalized, 
                     aes(ymin = -normalized, ymax = 0),
                     color=my_palette3[2]) +
      geom_hline(yintercept=0
                 # linetype="dashed", color = "red", size=2
                 ) + 
      labs(x = "m/z",
           # title = "U13C-glucose + 12C-thiamine",
           y = "Relative intensity") + 
      expand_limits(x = c(0,300),
                    y = c(0,100)) + 
      scale_y_continuous(expand = c(0,0),
                         breaks = scales::pretty_breaks(n = 4),
                         labels = c(100, 50, 0, 50, 100)
      ) +
      scale_x_continuous(expand = c(0,0)) +
      theme_classic(base_size = 12 # edit font size for all non-data text
      ) + 
      theme(plot.title = element_text(size = 12, hjust = 0.5),
            plot.margin = margin(0.2,0.2,0.2,0.2,"inch"),
            axis.text.x = element_text(angle = 0, hjust = .5, vjust = .7)
            # axis.ticks.x = element_blank()
      )
    
    ms2_neg = ggplot(liver_MS2_neg_normalized,aes(x=Centroid_mz, y= normalized)) +
      geom_linerange(aes(ymin = 0, ymax = normalized),
                     color=my_palette3[1]) +
      geom_linerange(data = synthesis_MS2_neg_normalized, 
                     aes(ymin = -normalized, ymax = 0),
                     color=my_palette3[2]) +
      geom_hline(yintercept=0
                 # linetype="dashed", color = "red", size=2
      ) + 
      labs(x = "m/z",
           # title = "U13C-glucose + 12C-thiamine",
           y = "Relative intensity") + 
      expand_limits(x = c(0,300),
                    y = c(0,100)) + 
      scale_y_continuous(expand = c(0,0),
                         breaks = scales::pretty_breaks(n = 4),
                         labels = c(100, 50, 0, 50, 100)
      ) +
      scale_x_continuous(expand = c(0,0)) +
      theme_classic(base_size = 12 # edit font size for all non-data text
      ) + 
      theme(plot.title = element_text(size = 12, hjust = 0.5),
            plot.margin = margin(0.2,0.2,0.2,0.2,"inch"),
            axis.text.x = element_text(angle = 0, hjust = .5, vjust = .7)
            # axis.ticks.x = element_blank()
      )
    
    
  }
  
  ## Output ####
  {
    ggpubr::ggarrange(
      ms2_pos,
      ms2_neg,
      common.legend = T, legend = "none",
      align = "hv",
      nrow = 2, ncol = 1
    ) %>%
      ggexport(filename = "plot_5_2.pdf", width = 3.5, height = 3.8)
  }
}
## plot_5_3 - isotope labeling for glucosyltaurine ####
{
  setwd(rel_path)
  setwd("plot_5 - glucosyltaurine")
  
  ## Data ####
  {
    filename = "glucosyl-taurine-raw-data-0902.xlsx"
    raw_data = read_xlsx(filename, sheet = "labeling") %>%
      select_if(~sum(!is.na(.)) > 0)
    colnames(raw_data)[1:2] = c("category", "cohort")
    
    
    # Filter data #
    labeling = raw_data %>%
      gather(key = "key", value = "number", -category, -cohort) %>%
      # Add an enter to name
      mutate(category = gsub(" ", "\\\n", category)) %>%
      # Filter out unwanted data
      # filter(!grepl("C0", cohort)) %>%
      group_by(category, key) %>%
      mutate(cum_number = cumsum(number),
             normalized = number/max(cum_number),
             cum_normalized = cumsum(normalized)) %>%
      filter(T)
    
    labeling_summary = labeling %>%
      group_by(category, cohort) %>%
      summarise(
        sd = sd(normalized),
        mean = mean(normalized)
      ) %>%
      ungroup() %>%
      group_by(category) %>%
      mutate(cum_mean = cumsum(mean)) %>%
      filter(T)
    
    labeling_merge = merge(labeling, labeling_summary) %>%
      mutate(datapoint = cum_mean - mean + normalized)
    
    
    
  }
  
  ## Plot ####
  {
    labeling_plot = ggplot(labeling_summary, aes(y = mean, x = category, fill = cohort)) +
      geom_bar(stat = "identity",
               width = .8,
               # position = position_identity(),
               position = position_stack(reverse=T),
               colour = "#333333"
      ) +
      geom_errorbar(aes(ymin=cum_mean-sd, ymax=cum_mean+sd),
                    width = 0.8,
                    # size = 2,
                    # width = 0.5,
                    position=position_dodge(0)) +
      # geom_point(data = labeling_merge,
      #            aes(y = datapoint),
      #            position = position_jitter(width = 0.2),
      #            # alpha = 0.5,
      #            shape = 1,
      #            show.legend=FALSE
      # ) +
      guides(fill = guide_legend(
        title = NULL,
        reverse = T
      )) +
      labs(x = NULL,
           # title = "U13C-glucose",
           y = "Labeling fraction") +
      scale_y_continuous(limits = c(0, 1.05),
                         expand = c(0,0),
                         breaks = c(0, 0.25, 0.5, 0.75, 1)
      ) +
      coord_cartesian(ylim = c(0.5, 1.05)) + # Zoom in without discarding data
      scale_x_discrete(limits = c("serum\nglucose", 
                                  "liver\nglucose",
                                  "liver\nG6P",
                                  "liver\nglucosyl-taurine"
                                  )) +
      scale_fill_manual(values = c(my_palette2)) + 
      scale_color_manual(values = "grey") + 
      # facet_wrap(~cohort) +
      theme_classic(base_size = 12 # edit font size for all non-data text
      ) +
      theme(plot.title = element_text(size = 12, hjust = 0.5),
            plot.margin = margin(0.5,0.5,0.5,0.5,"cm"),
            axis.text.x = element_text(angle = 0, hjust = .5, vjust = .7)
            # axis.ticks.x = element_blank()
      )
    print(labeling_plot)
  }
  
  
  ## Output ####
  {
    ggpubr::ggarrange(
      labeling_plot,
      
      common.legend = T, legend = "right",
      # align = "hv",
      nrow = 1, ncol = 1
    ) %>%
      ggexport(filename = "plot_5_3_zoom.pdf", width = 4.5, height = 2.5)
  }
  
}


## plot_5_4 - concentration in tissue ####
{
  setwd(rel_path)
  setwd("plot_5 - glucosyltaurine")
  
  ## Data ####
  {
    filename = "glucosyl-taurine-raw-data-0902.xlsx"
    raw_data = read_xlsx(filename, sheet = "tissue concentration") %>%
      select_if(~sum(!is.na(.)) > 0)
    colnames(raw_data)[1] = c("category")
    
    
    # Filter data #
    df = raw_data %>%
      gather(key = "key", value = "number", -category) %>%
      # Add an enter to name
      mutate(category = gsub(" ", "\\\n", category)) %>%
      # Filter out unwanted data
      # filter(!grepl("C0", cohort)) %>%
      # filter(grepl("liver|kidney|brain|pancreas|serum", category)) %>%
      filter(T)
    
    df_summary = df %>%
      group_by(category) %>%
      summarise(
        sd = sd(number),
        mean = mean(number),
        n = n(),
        se = sd/sqrt(n)
      ) %>%
      filter(T)
      
    df_merge = merge(df, df_summary)
    
    
    
  }
  
  ## Plot ####
  {
    df_plot = ggplot(df_summary, aes(y = mean, x = category)) +
      geom_bar(stat = "identity",
               width = .8,
               # position = position_identity(),
               position = position_stack(reverse=T),
               # fill = my_palette2[1],
               fill = "#888888",
               colour = "#333333"
      ) +
      geom_errorbar(aes(ymin=mean-se, ymax=mean+se),
                    width = 0.2,
                    # size = 2,
                    # width = 0.5,
                    position=position_dodge(0)) +
      geom_point(data = df_merge,
                 aes(y = number),
                 position = position_jitter(width = 0.2),
                 # alpha = 0.5,
                 shape = 1,
                 show.legend=FALSE
      ) +
      guides(fill = guide_legend(
        title = NULL,
        reverse = T
      )) +
      labs(x = NULL,
           # title = "U13C-glucose",
           y = "μM") +
      scale_y_continuous(expand = c(0,0),
                         breaks = scales::pretty_breaks(n = 4)
                         # breaks = c(0, 0.25, 0.5, 0.75, 1)
      ) +
      coord_cartesian(ylim = c(0, 205)) + # Zoom in without discarding data
      # scale_x_discrete(limits = c("serum\nglucose", "liver\nglucosyl-taurine")) +
      scale_fill_manual(values = c(my_palette2)) + 
      scale_color_manual(values = "grey") + 
      # facet_wrap(~cohort) +
      theme_classic(base_size = 12 # edit font size for all non-data text
      ) +
      theme(plot.title = element_text(size = 12, hjust = 0.5),
            plot.margin = margin(0.5,0.5,0.5,0.5,"cm"),
            axis.text.x = element_text(angle = 45, hjust = .5, vjust = .7)
            # axis.ticks.x = element_blank()
      )
    print(df_plot)
  }
  
  
  ## Output ####
  {
    ggpubr::ggarrange(
      df_plot,
      
      common.legend = T, legend = "right",
      # align = "hv",
      nrow = 1, ncol = 1
    ) %>%
      ggexport(filename = "plot_5_4_lighter.pdf", width = 5.5, height = 2.5)
  }
  
}


## plot_5_sup1 - kinetic of glucosyltaurine ####
{
  setwd(rel_path)
  setwd("plot_5 - glucosyltaurine")
  
  ## Data ####
  {
    filename = "glucosyl-taurine-raw-data-0902.xlsx"
    std_mix = read_xlsx(filename, sheet = "supplemental_mix of std") %>%
      select_if(~sum(!is.na(.)) > 0) %>%
      mutate(category = "std_mix")
    liver_spike_U13Cglc = read_xlsx(filename, sheet = "supplemental_13C6 glucose") %>%
      select_if(~sum(!is.na(.)) > 0) %>%
      mutate(category = "liver_spike_U13Cglc")
    liver = read_xlsx(filename, sheet = "supplemental_liver extract") %>%
      select_if(~sum(!is.na(.)) > 0) %>%
      mutate(category = "liver")
    
    condition1 = bind_rows(std_mix[,c("time", "404020-pH7", "category")], 
                           liver_spike_U13Cglc[,c("time", "404020-pH7", "category")], 
                           liver[,c("time", "404020-pH7", "category")])
    condition2 = bind_rows(std_mix[,c("time", "50:50 MeOH:H2O", "category")], 
                           liver_spike_U13Cglc[,c("time", "50:50 MeOH:H2O", "category")], 
                           liver[,c("time", "50:50 MeOH:H2O", "category")])
  }
  ## Plot ####
  {
    my_palette4 = c(change_HSL(brewer.pal(8, "RdYlBu")[c(2)], S_change = -0.5),
                    change_HSL(brewer.pal(8, "RdYlBu")[c(2)], S_change = 1),
                    brewer.pal(11, "RdYlBu")[c(10)],
                    change_HSL(brewer.pal(11, "RdYlBu")[c(10)], S_change = 1)
                    )
    image(1:length(my_palette4), 1, as.matrix(1:length(my_palette4)), 
          col=my_palette4,
          xlab="", ylab = "", xaxt = "n", yaxt = "n", bty = "n")
    
    condition1_plot = ggplot(condition1, aes(x = time, y = `404020-pH7`, group = category)) +
      geom_point(color = my_palette4[2], aes(shape = category)) + 
      scale_shape_manual(values=c(16, 2, 4))+ 
      geom_line(color = my_palette4[1]) + 
      labs(x = "hours",
           # title = "U13C-glucose",
           y = "μM") +
      guides(shape = guide_legend(
        title = NULL,
        reverse = F
      )) +
      scale_y_continuous(expand = c(0,0),
                         breaks = scales::pretty_breaks(n = 4)
                         # breaks = c(0, 0.25, 0.5, 0.75, 1)
      ) +
      coord_cartesian(ylim = c(0, 4)) + # Zoom in without discarding data
      # facet_wrap(~cohort) +
      theme_classic(base_size = 12 # edit font size for all non-data text
      ) +
      theme(plot.title = element_text(size = 12, hjust = 0.5),
            plot.margin = margin(0.5,0.5,0.5,0.5,"cm"),
            axis.text.x = element_text(angle = 0, hjust = .5, vjust = .7)
            # axis.ticks.x = element_blank()
      )
    print(condition1_plot)
    
    condition2_plot = ggplot(condition2, aes(x = time, y = `50:50 MeOH:H2O`, group = category)) +
      geom_point(color = my_palette4[4], aes(shape = category)) + 
      scale_shape_manual(values=c(16, 2, 4))+ 
      geom_line(color = my_palette4[3]) + 
      labs(x = "hours",
           # title = "U13C-glucose",
           y = "μM") +
      guides(shape = guide_legend(
        title = NULL,
        reverse = F
      )) +
      scale_y_continuous(expand = c(0,0),
                         breaks = scales::pretty_breaks(n = 4)
                         # breaks = c(0, 0.25, 0.5, 0.75, 1)
      ) +
      coord_cartesian(ylim = c(0, 4)) + # Zoom in without discarding data
      # facet_wrap(~cohort) +
      theme_classic(base_size = 12 # edit font size for all non-data text
      ) +
      theme(plot.title = element_text(size = 12, hjust = 0.5),
            plot.margin = margin(0.5,0.5,0.5,0.5,"cm"),
            axis.text.x = element_text(angle = 0, hjust = .5, vjust = .7)
            # axis.ticks.x = element_blank()
      )
      
    print(condition2_plot)
      
    condition1_plot = condition1_plot + 
    scale_y_continuous(expand = c(0,0),
                       breaks = scales::pretty_breaks(n = 2)
                       # breaks = c(0, 0.25, 0.5, 0.75, 1)
    )
    condition2_plot = condition2_plot + 
      scale_y_continuous(expand = c(0,0),
                         breaks = scales::pretty_breaks(n = 2)
                         # breaks = c(0, 0.25, 0.5, 0.75, 1)
      )
    
    condition1_zoom = condition1_plot +
      coord_cartesian(ylim = c(0, 0.06)) +
      scale_y_continuous(breaks = scales::pretty_breaks(n = 3))
    
    condition2_zoom = condition2_plot +
      coord_cartesian(ylim = c(0, 0.06)) +
      scale_y_continuous(breaks = scales::pretty_breaks(n = 3)) +
      coord_flip()
      
  }
  
  ## Output ####
  {
    ggpubr::ggarrange(
        condition1_plot,
        condition2_plot,
        condition1_zoom,
        condition2_zoom,
        common.legend = T, legend = "right",
        align = "hv",
        nrow = 2, ncol = 2
      ) %>%
        ggexport(filename = "plot_5_sup1_4.pdf", width = 7, height = 3)
    # ggpubr::ggarrange(
    #   condition1_plot,
    #   condition2_plot,
    #   common.legend = T, legend = "right",
    #   align = "hv",
    #   nrow = 1, ncol = 2
    # ) %>%
    #   ggexport(filename = "plot_5_sup1.pdf", width = 6, height = 2)
    # 
    # ggpubr::ggarrange(
    #   condition1_zoom,
    #   condition2_zoom,
    #   common.legend = T, legend = "right",
    #   align = "hv",
    #   nrow = 1, ncol = 2
    # ) %>%
    #   ggexport(filename = "plot_5_sup1_2.pdf", width = 5, height = 1.5)
  }
  
}
## plot_6 - network property of yeats neg data ####
{
  setwd(rel_path)
  setwd("plot_6_network property")
  filename = "20200704085443_output.rds" # Sc_neg
  
  # Global parameter in the background ####
  {
    dt = readRDS(filename)
    
    ilp_nodes = dt$ilp_nodes %>%
      mutate(medMz = signif(medMz, 7),
             medRt = round(medRt, 2),
             log10_inten = round(log10_inten, 2),
             ppm_error = round(ppm_error, 2))
    
    ilp_nodes_ilp = ilp_nodes %>%
      filter(ilp_result == 1) %>%
      filter(class != "Unknown") %>%
      filter(T)
    
    ilp_edges = dt$ilp_edges
    
    g_met = initiate_g_met(ilp_nodes, ilp_edges)
    
    core_met = ilp_nodes %>%
      filter(steps == 0) %>%
      filter(class == "Metabolite")
    
    core_nonmet = ilp_nodes %>%
      filter(steps %% 1 == 0) %>%
      filter(class != "Unknown") 
    
    g_nonmet = initiate_g_nonmet(ilp_nodes, dt$ilp_edges, dt$heterodimer_ilp_edges)
    
    # ilp_edges_annotate_met = igraph::as_data_frame(g_met, "edges")
    # ilp_edges_annotate_nonmet = igraph::as_data_frame(g_nonmet, "edges")
    core_rank = core_annotate(ilp_nodes, dt$FormulaSet_df, dt$LibrarySet)
    
  }
  
  # Further filtering ####
  {
    g_met2_node = igraph::as_data_frame(g_met, "vertices") %>%
      filter(ilp_result > 0.01)
    
    g_met2_edges = igraph::as_data_frame(g_met, "edges") %>%
      filter(from %in% g_met2_node$name & to %in% g_met2_node$name)
    
    g_met2 = graph_from_data_frame(g_met2_edges,
                                   vertices = g_met2_node,
                                   directed = F)
    
    g_nonmet2_node = igraph::as_data_frame(g_nonmet, "vertices") %>%
      filter(ilp_result > 0.01)
    g_nonmet2_edges = igraph::as_data_frame(g_nonmet, "edges") %>%
      filter(from %in% g_nonmet2_node$name & to %in% g_nonmet2_node$name,
             ilp_result > 0.01) 
    g_nonmet2 = graph_from_data_frame(g_nonmet2_edges,
                                      vertices = g_nonmet2_node,
                                      directed = F)
    
    g_all_nodes = bind_rows(g_met2_node, g_nonmet2_node) %>%
      distinct(name, .keep_all = T) %>%
      # filter(log10_inten >= 5) %>%
      filter(class != "Unknown")
    g_all_edges = bind_rows(g_met2_edges, g_nonmet2_edges) %>%
      distinct() %>%
      filter(from %in% g_all_nodes$name & to %in% g_all_nodes$name)
  }
  
  # Sub-graph ####
  {
    # Problem 1 - many single node in the output
    # Reason - Heterodimer edge score enables nodes without regular connections
    # Solutions - Remove nodes do not have regular connections
    
    g_all_nodes = g_all_nodes %>%
      filter(name %in% g_all_edges$from | name %in% g_all_edges$to)
    
    g_all = graph_from_data_frame(g_all_edges,
                                  directed = F,
                                  vertices = g_all_nodes)
    
    # Problem 2 - Putative metabolites in network without connecting to metabolites
    # Reason - Edge score by connection with isotope/adduct etc.
    # Solutions - Only retain subnetwork where step=0 annotation exist
    g_sub=g_all
    clu=components(g_sub)
    #subnetwork criteria 
    mainnetwork = igraph::groups(clu)[table(clu$membership) == max(table(clu$membership))]
    subnetwork = igraph::groups(clu)[table(clu$membership) < length(mainnetwork[[1]])]
    
    
    step0 = g_all_nodes %>%
      filter(steps == 0) %>%
      pull(name)
    subnetwork_valid = sapply(subnetwork, function(x){
      any(x %in% step0)
    })
    nodes_invalid = unlist(subnetwork[!subnetwork_valid])
    
    g_all_nodes_valid = g_all_nodes %>%
      filter(!name %in% nodes_invalid)
    g_all_edges_valid = g_all_edges %>%
      filter(from %in% g_all_nodes_valid$name & to %in% g_all_nodes_valid$name)
    g_all_valid = graph_from_data_frame(g_all_edges_valid,
                                        directed = F,
                                        vertices = g_all_nodes_valid)
    
    nodes_metabolite = g_all_nodes %>%
      filter(class == "Metabolite") %>%
      pull(name)
    
    subnetwork_valid_nonmet = sapply(subnetwork, function(x){
      (any(x %in% step0)) & (all(x %in% nodes_metabolite))
    })
    test = subnetwork[subnetwork_valid_nonmet]
    
    test2 = g_all_nodes_valid %>%
      filter(name %in% test$`22`)
    
    
    # g_main = make_ego_graph(graph=g_sub,
    #                         order=diameter(g_sub),
    #                         nodes = mainnetwork[[1]][1],
    #                         mode="all")[[1]]
    
  }
  
  
  
}
## Merge graphs ####
{
  
  
}

## Good ggplot graph ####
{
  ggplot(diamonds, aes(x=cut, y=price, group=cut))+
    geom_boxplot(aes(fill=cut))+scale_fill_brewer(palette="OrRd")
}
