# install.packages("visNetwork")
library(visNetwork)

nodes <- data.frame(id = 1:3)
edges <- data.frame(from = c(1,2), to = c(1,3))
visNetwork(nodes, edges, width = "100%")

nodes <- data.frame(id = 1:10, 
                    label = paste("Node", 1:10),                                 # add labels on nodes
                    group = c("GrA", "GrB"),                                     # add groups on nodes 
                    value = 1:10,                                                # size adding value
                    shape = c("square", "triangle", "box", "circle", "dot", "star",
                              "ellipse", "database", "text", "diamond"),                   # control shape of nodes
                    title = paste0(1:10,"<br>","node!"),         # tooltip (html or character)
                    color = c("darkred", "grey", "orange", "darkblue", "purple"),# color
                    shadow = c(FALSE, TRUE, FALSE, TRUE, TRUE))                  # shadow

head(nodes)
edges <- data.frame(from = sample(1:10, 8), to = sample(1:10, 8),
                    label = paste("Edge", 1:8),                                 # add labels on edges
                    length = c(100,500),                                        # length
                    arrows = c("to", "from", "middle", "middle;to"),            # arrows
                    dashes = c(TRUE, FALSE),                                    # dashes
                    title = paste("Edge", 1:8),                                 # tooltip (html or character)
                    smooth = c(FALSE, TRUE),                                    # smooth
                    shadow = c(FALSE, TRUE, FALSE, TRUE))                       # shadow
head(edges)

visNetwork(nodes, edges, width = "100%")


nodes <- data.frame(id = 1:5, group = c(rep("A", 2), rep("B", 3)))
edges <- data.frame(from = c(2,5,3,3), to = c(1,2,4,2))

visNetwork(nodes, edges, width = "100%") %>% 
  visNodes(shape = "square") %>%                        # square for all nodes
  visEdges(arrows ="to") %>%                            # arrow "to" for all edges
  visGroups(groupname = "A", color = "darkblue") %>%    # darkblue for group "A"
  visGroups(groupname = "B", color = "red")             # red for group "B"


nb <- 10
nodes <- data.frame(id = 1:nb, label = paste("Label", 1:nb),
                    group = sample(LETTERS[1:3], nb, replace = TRUE), value = 1:nb,
                    title = paste0("<p>", 1:nb,"<br>Tooltip !</p>"), stringsAsFactors = FALSE)

edges <- data.frame(from = trunc(runif(nb)*(nb-1))+1,
                    to = trunc(runif(nb)*(nb-1))+1,
                    value = rnorm(nb, 10), label = paste("Edge", 1:nb),
                    title = paste0("<p>", 1:nb,"<br>Edge Tooltip !</p>"))


visNetwork(nodes, edges, width = "100%") %>% visLegend()
visNetwork(nodes, edges, width = "100%", height="100%") %>% 
  visLegend(useGroups = FALSE, addNodes = data.frame(label = "Nodes", shape = "circle"), 
            addEdges = data.frame(label = "link", color = "black"))


visNetwork(nodes, edges, width = "100%") %>% 
  visOptions(highlightNearest = TRUE)

visNetwork(nodes, edges, width = "100%") %>% 
  visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE)

nodes$sel <- sample(c("sel1", "sel2"), nrow(nodes), replace = TRUE)
visNetwork(nodes, edges, width = "100%") %>%
  visOptions(selectedBy = "sel")



visNetwork(nodes, edges, width = "100%") %>% 
  visInteraction(navigationButtons = TRUE)

visNetwork(nodes, edges, width = "100%") %>% 
  visOptions(manipulation = TRUE)
