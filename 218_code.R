library(remotes)
#remotes::install_github("ropensci/rmangal")
library(rmangal)
help(rmangal)
library(rlang)
library(tidygraph)
library(ggraph)
library(igraph)
library(intergraph)
library(network)
library(ergm)
library(sna)
library(degreenet)
library(boot)

mgs <- search_datasets("Plant and animals interaction in Long Islands estuary")
mgn <- get_collection(mgs)
long_islands_ig <- as.igraph(mgn)
#long_islands_net <- as.network(long_islands_ig)

ilogit <- function(x){
   exp(x)/(1+exp(x))
}

###########################################################
# DATA PREPROCESSING PART 1: Cleaning up Vertex Attributes#
###########################################################

# delete the irrelevant vertex attributes
long_islands_ig <- delete_vertex_attr(long_islands_ig, "node_level")
long_islands_ig <- delete_vertex_attr(long_islands_ig, "network_id")
long_islands_ig <- delete_vertex_attr(long_islands_ig, "created_at")
long_islands_ig <- delete_vertex_attr(long_islands_ig, "updated_at")
long_islands_ig <- delete_vertex_attr(long_islands_ig, "taxonomy.ncbi")
long_islands_ig <- delete_vertex_attr(long_islands_ig, "taxonomy.tsn")
long_islands_ig <- delete_vertex_attr(long_islands_ig, "taxonomy.eol")
long_islands_ig <- delete_vertex_attr(long_islands_ig, "taxonomy.bold")
long_islands_ig <- delete_vertex_attr(long_islands_ig, "taxonomy.gbif")
long_islands_ig <- delete_vertex_attr(long_islands_ig, "taxonomy.col")
long_islands_ig <- delete_vertex_attr(long_islands_ig, "taxonomy.created_at")
long_islands_ig <- delete_vertex_attr(long_islands_ig, "taxonomy.updated_at")
long_islands_ig <- delete_vertex_attr(long_islands_ig, "taxonomy")
long_islands_ig <- delete_vertex_attr(long_islands_ig, "taxonomy.id")

vertex_attr(long_islands_ig)

# remove the N/A values for taxonomy_id and replace with 0 
z <- vertex_attr(long_islands_ig)["taxonomy_id"]
z <- rapply(z, f=function(x) ifelse(is.na(x),0,x), how="replace" )
vertex_attr(long_islands_ig)["taxonomy_id"] <- z


# remove the N/A values for taxonomy.name and replace with "None" 
z <- vertex_attr(long_islands_ig)["taxonomy.name"]
z <- rapply(z, f=function(x) ifelse(is.na(x), "None" ,x), how="replace" )
vertex_attr(long_islands_ig)["taxonomy.name"] <- z


# remove the N/A values for taxonomy.rank and replace with "None" 
z <- vertex_attr(long_islands_ig)["taxonomy.rank"]
z <- rapply(z, f=function(x) ifelse(is.na(x), "None" ,x), how="replace" )
vertex_attr(long_islands_ig)["taxonomy.rank"] <- z


###########################################################
# DATA PREPROCESSING PART 2: Cleaning up Edge Attributes###
###########################################################

# delete the irrelevant edge attributes
long_islands_ig <- delete_edge_attr(long_islands_ig, "date")
long_islands_ig <- delete_edge_attr(long_islands_ig, "direction")
long_islands_ig <- delete_edge_attr(long_islands_ig, "method")
long_islands_ig <- delete_edge_attr(long_islands_ig, "attr_id")
long_islands_ig <- delete_edge_attr(long_islands_ig, "value")
long_islands_ig <- delete_edge_attr(long_islands_ig, "public")
long_islands_ig <- delete_edge_attr(long_islands_ig, "network_id")
long_islands_ig <- delete_edge_attr(long_islands_ig, "created_at")
long_islands_ig <- delete_edge_attr(long_islands_ig, "updated_at")
long_islands_ig <- delete_edge_attr(long_islands_ig, "attribute.id")
long_islands_ig <- delete_edge_attr(long_islands_ig, "attribute.name")
long_islands_ig <- delete_edge_attr(long_islands_ig, "attribute.description")
long_islands_ig <- delete_edge_attr(long_islands_ig, "attribute.unit")
long_islands_ig <- delete_edge_attr(long_islands_ig, "attribute.created_at")
long_islands_ig <- delete_edge_attr(long_islands_ig, "attribute.updated_at")

vertex_attr(long_islands_ig)["name"] <- vertex_attr(long_islands_ig)["original_name"]

# using incident edges to determine if the organism acts as a herbivore
# based on this data, create a vertex attribute "herbivore" such that
# 0 = node does not function as a herbivore, 1 = node functions as a herbivore
# as a double vector. Likewise for detritivore.

incident_edges(long_islands_ig, 'plankton', mode="out") # sink node 
incident_edges(long_islands_ig, 'water plant', mode="out") # sink node
incident_edges(long_islands_ig, 'bay shrimp', mode="out") 
incident_edges(long_islands_ig, 'silversides', mode="out") 
incident_edges(long_islands_ig, 'mud snail', mode="out") 
incident_edges(long_islands_ig, 'clam', mode="out") 
incident_edges(long_islands_ig, 'billfish', mode="out")
incident_edges(long_islands_ig, 'eel', mode="out") 
incident_edges(long_islands_ig, 'blowfish', mode="out") 
incident_edges(long_islands_ig, 'minnow 1', mode="out") 
incident_edges(long_islands_ig, 'minnow 2', mode="out") 
incident_edges(long_islands_ig, 'fluke', mode="out") 
incident_edges(long_islands_ig, 'cricket', mode="out") 
incident_edges(long_islands_ig, 'mosquito', mode="out")
incident_edges(long_islands_ig, 'tern', mode="out")
incident_edges(long_islands_ig, 'Osprey', mode="out")
incident_edges(long_islands_ig, 'Green Heron', mode="out")
incident_edges(long_islands_ig, 'merganser', mode="out")
incident_edges(long_islands_ig, 'cormorant', mode="out")
incident_edges(long_islands_ig, 'gull', mode="out")
incident_edges(long_islands_ig, 'kingfisher', mode="out")
incident_edges(long_islands_ig, 'Red-winged Blackbird', mode="out")
incident_edges(long_islands_ig, 'organic debris', mode="out") # sink node
incident_edges(long_islands_ig, 'marsh plants', mode="out") # sink node

herbivore <- c(0, 0, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 
               0, 0, 0, 0, 0, 0, 0)
detritivore <- c(0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)

herbivore <- unlist(herbivore)
detritivore <- unlist(detritivore)
long_islands_ig <- set_vertex_attr(graph=long_islands_ig, name="herbivore", value=herbivore)
long_islands_ig <- set_vertex_attr(graph=long_islands_ig, name="detritivore", value=detritivore)

# we will then proceed to remove all of the sink nodes as these nodes do not 
# act as one of the designated consumer types: herbivore, predator, detritivore
# the sink nodes are: plankton (id=1), water plant (id=2), 
# organic debris (id=23), marsh plants (id=24) which we will now remove

long_islands_ig <- delete_vertices(long_islands_ig, "plankton")
long_islands_ig <- delete_vertices(long_islands_ig, "water plant")
long_islands_ig <- delete_vertices(long_islands_ig, "organic debris")
long_islands_ig <- delete_vertices(long_islands_ig, "marsh plants")


vlab <- vertex_attr(long_islands_ig, name="name")
elab <- edge_attr(long_islands_ig, name="type")

# we now have our predation network, with added vertex (nodal) attributes of herbivore/detritivore
plot(long_islands_ig, edge.arrow.size = 0.2, vertex.label = vlab,
     vertex.label.cex = 1,
     edge.label = elab, 
     edge.label.cex = 0.7,
     edge.label.color="darkred")

long_islands_net <- intergraph::asNetwork(long_islands_ig, matrix.type="edgelist")
plot(long_islands_net,vertex.col="herbivore")
plot(long_islands_net,vertex.col="detritivore")

vertex_attr(long_islands_ig)
edge_attr(long_islands_ig)

# out-degree network sequence
odeg <- degree(long_islands_net, cmode = "outdegree")
barplot(odeg, xlab = 'Organisms', ylab = 'Out Degree')
summary(odeg)

# in-degree network sequence
ideg <- sna::degree(long_islands_net, cmode = "indegree")
barplot(ideg, xlab = 'Organisms', ylab = 'In Degree')
summary(ideg)

# high out-degree means low in-degree (no one eats the big predators)
cor(ideg, odeg)

# Fit various degree distribution models:
# Pareto/Zipf law, Yule, Waring, Poisson, and Conway-Maxwell-Poisson models.
help(ayulemle)

m1 <- ayulemle(deg) # Yule model
m2 <- adpmle(odeg) # Pareto/Zipf model
m3 <- awarmle(odeg) # Waring model
m4 <- apoimle(odeg) # Poisson model
m5 <- acmpmle(odeg) # Conway-Maxwell-Poisson model

# Comparing fits of the degree distribution models

# Estimate the Yule Discrete Probability Distribution via maximum likelihood
llyuleall(m1$theta, odeg)
# Estimate the Pareto/Zipf Discrete Probability Distribution via maximum likelihood
lldpall(m2$theta, odeg)
# Estimate the Waring Discrete Probability Distribution via maximum likelihood
llwarall(m3$theta, odeg)
# Estimate the Poisson Discrete Probability Distribution via maximum likelihood
llpoiall(m4$theta, odeg)
# Estimate the Conway Maxwell Poisson Discrete Probability Distribution via maximum likelihood
llcmpall(m5$theta, odeg)

# Degree correlation is also known as the assortativity coefficient.
# It measures whether high-degree nodes tend to connect to other
# high-degree nodes, or in other words measures the homophily of nodes.
assortativity_degree(asIgraph(long_islands_net), directed = T)

# Network centrality measures
detach(package: igraph)
a <- sna::centralization(long_islands_net, degree, mode="digraph", diag=FALSE)
b <- sna::centralization(long_islands_net, closeness, mode="digraph", diag=FALSE)
c <- sna::centralization(long_islands_net, betweenness, mode="digraph", diag=FALSE)
d <- sna::centralization(long_islands_net, evcent, mode="digraph", diag=FALSE)
c(a, b, c, d)

# Interpretation model_a: We measure if an organism can get eaten by k=2 other predators using
# istar(). As such, we can determine if there are commonly consumed prey. We measure if a predator 
# can eat k=2 other organisms using ostar(). As such, we can determine if predators are clustered 
# in terms of eating habits/who they prey on. By Monte Carlo MLE Results, istar2 = 0.1331 (ilogit = 0.533226) so there is 
# evidence of commonly 
# consumed prey. ostar2 = 0.2413 (ilogit = 0.560034), so there is also evidence of clustered predators in terms of consumed
# prey. ostar2 is larger than istar2 indicating that clustered predators comsuming 2 different prey are 
# more likely than prey getting consumed by 2 different predators. 
# TO DO: Show the MCMC model plot and analyze the models and GOF figures

model_a <- ergm(long_islands_net~edges + istar(2) + ostar(2))
summary(model_a)
mcmc.diagnostics(model_a)
plot(gof(model_a))
# ilogit = inverse log likelihood
ilogit(0.1331) # istar2
ilogit(0.2413) # ostar2

##############################################################

# Interpretation model_b: We measure if an organism can get eaten by k=2 other predators using
# istar(). As such, we can determine if there are commonly consumed prey. We measure if a predator 
# can eat k=2 other organisms using ostar(). We measure differential homophily for within
# group ties where the groups are herbivore/not herbivore, and 
# detritivore/not detritivore. By Monte Carlo MLE Results, istar2 = 0.14117 (ilogit = 0.535234), ostar2 = 0.23938 (ilogit = 0.5595609),
# nodematch.herbivore.0 = -0.64545 (ilogit = 0.3440156), nodematch.herbivore.1 = -1.27698 (ilogit = 0.2180647), 
# nodematch.detritivore.0 = -0.03471 (ilogit = 0.4913234), 
# nodematch.detritivore.1 = -Inf (ilogit = 0). Interpretations for istar2 and ostar2 are the same as model_a. 
# We see that non-herbivores are somewhat unlikely to consume
# other non-herbivores, herbivores are unlikely to consume herbivores, non-detritivores are slightly unlikely
# to consume other non-detritivores, and detritivores are very unlikely to consume other detritivores. 
# Our data shows no evidence of detritivores consuming other detritivores, hence it is -Inf. 
# Herbivores by definition cannot consume other herbivores, however because organisms in our network can be multiple
# consumer types (ie. both herbivore and predator where consumption homophily exists), 
# nodematch.herbivore.1 is not -Inf. Non-herbivores can be both predators
# and/or detritivores and nodes can act as multiple consumer types. As such eventhough predators do consume other predators,
# detritivores do not consume other detritivores, and herbivores do not consume other herbivores. So, nodematch.herbivore.0 
# is slightly negative, but close to 0. Similarly, non-detritivores can be herbivores and predators where predators have 
# consumption homophily 
# with other predators, but herbivores do not have consumption homophily with other herbivores. So nodematch.detritivore.0 is
# slightly negative but close to 0. 
# TO DO: Show the MCMC model plot and analyze the models and GOF figures

model_b <- ergm(long_islands_net~edges + istar(2) + ostar(2) + 
                   nodematch("herbivore", diff=T) + nodematch("detritivore", diff=T))
summary(model_b)
mcmc.diagnostics(model_b)
plot(gof(model_b))
# ilogit = inverse log likelihood
ilogit(0.14117) # istar2
ilogit(0.23938) # ostar2
ilogit(-0.64545) # nodematch.herbivore.0
ilogit(-1.27698) # nodematch.herbivore.1
ilogit(-0.03471) # nodematch.detritivore.0
ilogit(-Inf) # nodematch.detritivore.1

##############################################################

# Interpretation model_c: We measure if an organism can get eaten by k=2 other predators using
# istar(). As such, we can determine if there are commonly consumed prey. We measure if a predator 
# can eat k=2 other organisms using ostar(). We measure differential homophily for within
# group ties where the groups are herbivore/not herbivore, and 
# detritivore/not detritivore. We measure node covariance of herbivores and detritivores in the
# network. As all ties in the network are of type predation, node covariance is used to
# determine if an organism being a herbivore or detritivore makes it more or less likely
# to also be a predator. By Monte Carlo MLE Results, istar2 = 0.133457 (ilogit = 0.5333148), ostar2 = 0.241320 (ilogit = 0.5600389),
# nodematch.herbivore.0 = -0.635441 (ilogit = 0.3462778), nodematch.herbivore.1 = -1.256960 (ilogit = 0.2214977), 
# nodematch.detritivore.0 = -0.038392 (ilogit = 0.4904032),
# nodematch.detritivore.1 = -Inf (ilogit = 0), nodecov.herbivore = -0.008025 (ilogit = 0.4979938), 
# nodecov.detritivore = -0.007792 (ilogit = 0.498052).
# Interpretations for istar2, ostar2, nodematch.herbivore, and nodematch.detritivore are the same 
# as model_a and model_b. Based on the results of nodecov.herbivore and nodecov.detritivore, 
# we can see that herbivores are less likely to also be predators (outgoing tie), but is not -Inf, 
# because there are organisms that are both herbivores and predators in the data. We can see that
# detritivores are slightly likely to also be predators, but close to 0 (random chance). 
# TO DO: Show the MCMC model plot and analyze the models and GOF figures


model_c <- ergm(long_islands_net~edges + istar(2) + ostar(2) + 
                   nodematch("herbivore", diff=T) + nodematch("detritivore", diff=T) +
                   nodecov("herbivore") + nodecov("detritivore"))
summary(model_c)
mcmc.diagnostics(model_c)
plot(gof(model_c))
# ilogit = inverse log likelihood
ilogit(0.133457) # istar2
ilogit(0.241320) # ostar2
ilogit(-0.635441) # nodematch.herbivore.0
ilogit(-1.256960) # nodematch.herbivore.1
ilogit(-0.038392) # nodematch.detritivore.0
ilogit(-Inf) # nodematch.detritivore.1
ilogit(-0.008025) # nodecov.herbivore
ilogit(-0.007792) # nodecov.detritivore

##############################################################

# Interpretation model_d: We measure if an organism can get eaten by k=2 other predators using
# istar(). As such, we can determine if there are commonly consumed prey. We measure if a predator 
# can eat k=2 other organisms using ostar(). We measure differential homophily for within
# group ties where the groups are herbivore/not herbivore, and 
# detritivore/not detritivore. We measure node covariance of herbivores and detritivores in the
# network. As all ties in the network are of type predation, node covariance is used to
# determine if an organism being a herbivore or detritivore makes it more or less likely
# to also be a predator. We measure if there exists cyclical triad eating patterns with
# predators in ctriple. A ttriple exists if there are a set of directed edges such that
# {(i -> j), (j -> k), (k -> i)}. We measure if there exists heirarchical triad eating patterns with
# organisms in ttriple. A ttriple exists if there are a set of edges such that
# {(i -> j), (j -> k), (i -> k)}. By Monte Carlo MLE Results: istar2 = 0.1810 (ilogit = 0.5451269), 
# ostar2 = 0.2990 (ilogit = 0.574198), nodematch.herbivore.0 = -0.5799 (ilogit = 0.3589556), 
# nodematch.herbivore.1 = -1.201 (ilogit = 0.2312974), 
# nodematch.detritivore.0 = -0.06312 (ilogit = 0.4842252), nodematch.detritivore.1 = -Inf (ilogit = 0), 
# nodecov.herbivore = -0.007876 (ilogit = 0.498031), 
# nodecov.detritivore = -0.04420 (ilogit = 0.4889518), ctriple = -Inf (ilogit = 0), and ttriple = -0.4874 (ilogit = 0.3805063). 
# The interpretations for
# istar2, ostar2, nodematch.herbivore, nodematch.detritivore, and nodecov.herbivore
# are the same as in model_a, model_b, and model_c. 
# There is no evidence of ctriples existing in the data, so it is -Inf. While the data does not contain
# ttriples (hence the negative ttriple value), it does show evidence of heirachical eating. So, we will
# instead use shortest path finding algorithms to determine types of heirachical eating consumption paths. 
# TO DO: Show the MCMC model plot and analyze the models and GOF figures

model_d <- ergm(long_islands_net~edges + istar(2) + ostar(2) + 
                   nodematch("herbivore", diff=T) + nodematch("detritivore", diff=T) +
                   nodecov("herbivore") + nodecov("detritivore") +
                   ctriple() + ttriple())
summary(model_d)
mcmc.diagnostics(model_d)
plot(gof(model_d))
# ilogit = inverse log likelihood
ilogit(0.1810) # istar2
ilogit(0.2990) # ostar2
ilogit(-0.5799) # nodematch.herbivore.0
ilogit(-1.201) # nodematch.herbivore.1
ilogit(-0.06312) # nodematch.detritivore.0
ilogit(-Inf) # nodematch.detritivore.1
ilogit(-0.007876) # nodecov.herbivore
ilogit(-0.04420) # nodecov.detritivore
ilogit(-Inf) # ctriple
ilogit(-0.4874) # ttriple

graph <- long_islands_ig
distances(graph, v = V(graph), to = V(graph), mode = c("out"), weights = NULL, algorithm = c("dijkstra"))
# tern -> billfish -> bay shrimp
# osprey -> billfish -> bay shrimp
# osprey -> billfish -> silversides
# osprey -> blowfish -> mud snail
# osprey -> blowfish -> clam
# green heron -> billfish -> bay shrimp
# green heron -> billfish -> silversides
# merganser -> blowfish -> mud snail
# merganser -> blowfish -> clam



