# Need to run CompleteETV5Anlaysis.Rnw through calculation of 
# data_thresh (line 250)

I1=which(V(miceDESeqComplex$tumorg)$name %in% rownames(data_thresh))
HighBetween = induced_subgraph(miceDESeqComplex$tumorg,I1)
title1 = 'High Betweeness Subgraph'
#L = layout.auto(normalmax)
#L = layout.reingold(normalmax)
L1 = layout.fruchterman.reingold(HighBetween,niter=500,area=vcount(HighBetween)^2.3, repulserad=vcount(HighBetween)^2.8)
dev.new()
graphname = paste('Graphs/',title1,'2','.png') 
graphname=gsub(" ","",graphname,fixed=TRUE)   # remove spaces from graph name
png(graphname,width = 1600, height = 900)
plot(HighBetween,layout = L1,vertex.size=30,vertex.label.dist = 0,
     vertex.label.cex = .7, vertex.color='yellow',main = c(title1))
# make size of vertex proportional to betweeness
plot(HighBetween,layout = L,vertex.size=.00001*V(HighBetween)$between,
     vertex.label.dist = 0,vertex.label.cex = .8, 
     vertex.color='yellow',main = c(title1))
dev.off()

# For the next bit must have run  CompleteETV5Analysis.Rnw through line 405
Targs = tarGenesP$ETV5[c((tarGenesP$ETV5$padj < siglevel) & (!is.na(tarGenesP$ETV5$padj))),]
Tnames = Targs[,1]
I=which(V(miceDESeqComplex$tumorg)$name %in% Tnames)
E = which(V(miceDESeqComplex$tumorg)$name == "ETV5")
I = c(I,E)
colors = rep('yellow',length(I))
ETVtarg = induced_subgraph(miceDESeqComplex$tumorg,I)
colors[V(ETVtarg)$name == 'ETV5']='lavender'
L = layout.fruchterman.reingold(ETVtarg,niter=500,
                                area=vcount(ETVtarg)^2.3, 
                                repulserad=vcount(ETVtarg)^2.8)
title = "ETV and targets subgraph"
dev.new()
graphname = paste('Graphs/',title,'2','.png') 
graphname=gsub(" ","",graphname,fixed=TRUE)
png(graphname,width = 1600, height = 900)
plot(ETVtarg,layout = L,vertex.size=30,vertex.label.dist = 0,
     vertex.label.cex = .7, vertex.color=colors,main = c(title))
# make size of vertex proportional to betweeness
plot(ETVtarg,layout = L,vertex.size=.00001*V(ETVtarg)$between,
     vertex.label.dist = 0,vertex.label.cex = .8, 
     vertex.color='yellow',main = c(title))
dev.off()

dev.new()
title = "High Betweeness Nodes and ETV5 Targets"
graphname = paste('Graphs/',title,'2.pdf') 
pdf(graphname,width = 7, height = 7)
I2 = c(I,I1)
colors = rep('yellow',length(I2))
HiPlusETVtarg = induced_subgraph(miceDESeqComplex$tumorg,I2)
J = which(V(HiPlusETVtarg)$name %in% rownames(data_thresh))
colors[J]=rep('pink',length(rownames(data_thresh)))
colors[V(HiPlusETVtarg)$name == 'ETV5']='lavender'
 L2 = layout.auto(HiPlusETVtarg)
L2 = layout.fruchterman.reingold(HiPlusETVtarg, weights = 2*E(HiPlusETVtarg)$weight)
# L2 = layout_nicely(HiPlusETVtarg)
# L2 = layout.fruchterman.reingold(HiPlusETVtarg,niter=500,
#                                 area=vcount(HiPlusETVtarg)^2.3, 
#                                 repulserad=vcount(HiPlusETVtarg)^2.8)
plot(HiPlusETVtarg,layout = L2,vertex.size=22,vertex.label.dist = 0,
     vertex.label.cex = .7, vertex.color=colors,main = c(title),cex.main = 1.1)
dev.off()

# For a high-resolution tiff, used in PLosOne paper
#dev.new()
#tiff("Fig4.tiff", width = 6, height = 5, units = 'in', res = 500, compression = 'lzw')
L2 = layout.fruchterman.reingold(HiPlusETVtarg, weights = 3*E(HiPlusETVtarg)$weight)
tiff("Graphs/Fig4_3.tiff", width = 6, height = 5, units = 'in', res = 800)
plot(HiPlusETVtarg,layout = L2,vertex.size=16,vertex.label.dist = 0,
     vertex.label.cex = .35, vertex.color=colors)
dev.off()
