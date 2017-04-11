install.packages("igraph") 
install.packages("network")
install.packages("sna") 
install.packages("ndtv")

library("igraph")
library("network")
library("sna")
library("ndtv")

nodes<-read.csv("./networks/hawaii_vertices.csv")
stepedges<-read.csv("./networks/step_edges.csv",header=F)
colrs<-c("black","gray50")

#island model make all possible connections

edges<-t(combn(x=nodes$vertex,m=2))
# and the reverse
edges2<-rbind(edges,edges[,c(2,1)])

#make the n-island model
island<-graph.data.frame(edges2,nodes,directed=T)

circle<-layout_in_circle(island,nodes)

plot(island,layout=circle)


#make a stepping-stone model
stepstone<-graph.data.frame(stepedges,nodes,directed=T)
V(stepstone)$vertex.color<-colrs[V(stepstone)$sector]

plot(stepstone,layout=circle, edge.arrow.size=0.4, edge.curved=0.3, vertex.label=V(stepstone)$island, vertex.color=V(stepstone)$sector)

#make high-low

hilo<-contract(stepstone,mapping=c(1,1,1,1,2,2,2,2,2,2,2,2,2,2))
plot(hilo)
