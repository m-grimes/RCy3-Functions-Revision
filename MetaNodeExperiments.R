# Metanode experiments

# Mark Grimes
source("New_RCy3_Functions.R")
#------------------------------------------------------------
# Use network from FigureVignette (or just skip this section and open MetaNodeTestNetwork.cys)
# First pick two nodes. 
node1 <-  "ALK"
node2 <-  "CTTN"
nodenames=c(node1, node2)
connected.path <- connectNodes.all(nodenames, essl.cfn.ig, essl.cfn)
path.nodes <- unique(c(as.character(connected.path$Gene.1), as.character(connected.path$Gene.2)))
# Get peptides from this network
extract.peptides <- function(nodename, edgefile=pepnet.edges) {
    peps.1 <- edgefile$Gene.1[grep(nodename, edgefile$Gene.1, fixed=TRUE)]
    peps.2 <- edgefile$Gene.2[grep(nodename, edgefile$Gene.2, fixed=TRUE)]
    return(as.character(unique(c(peps.1, peps.2))))
}
netpeps <- as.character(unlist(sapply(path.nodes, extract.peptides, essl.cccn.edges)))  
# extract.peptides defined above
pepnodes.df <- data.frame(netpeps, pep.nodes=(sapply(netpeps, function(x) unlist(strsplit(x, " "))[1])))
netpeps <- as.character(pepnodes.df[pepnodes.df$pep.nodes %in% path.nodes, 1])
ptm.cccn <-	filter.edges.0(netpeps, essl.cccnplus.edges) 
net.full <- mergeEdges(connected.path)
net.full$Alt.Weight <- net.full$Weight
net.gene.cf <- make.anynetcf(edge.df=net.full, data.file=ldgene.fc[, crizexpts.nc], geneatts=essl.netatts, ptmcccnatts=essl.ptmcccnatts, func.key=func.key, use=c("total", "mean", "median", "max"))
net.pep.cf <- make.anynetcf(edge.df=ptm.cccn, data.file=ld.fc[, crizexpts.nc], geneatts=essl.netatts, ptmcccnatts=essl.ptmcccnatts, func.key=func.key, use=c("total", "mean", "median"))
# combine node attribute cytoscape file
net.cf <- harmonize_cfs3(pepcf=net.pep.cf, genecf=net.gene.cf)
# make gene-peptide edges
net.gpe <- data.frame(Gene.1=net.pep.cf$Gene.Name, Gene.2=net.pep.cf$Peptide.Name, edgeType="peptide", Weight=100, Alt.Weight=1)
# net.gpe <- genepep.edges(ptm.cccn)
# combine edge files
net.full <- net.full[, c("Gene.1", "Gene.2", "Weight", "edgeType", "Alt.Weight")]
net.edges <- rbind(net.gpe, ptm.cccn, net.full)
# Note: to change "edgeType" to "interaction" according to new requirements for createNetworkFromDataFrames, use edgeType.to.interaction() or
# names(net.edges)[grep("edgeType", names(net.edges))] <- "interaction"
# Graph in RCy3
#set data frame colnames for translation/import into Cytoscape
colnames(net.edges)[1:3]<-c("source", "target", "interaction")
# Note: can use Alt.Weight for edge visualization because peptide edges are set to 100 for layouts that cluster them next to nodes. If so: 
# net.edges$Weight <- net.edges$Alt.Weight
colnames(net.cf)[1]<-'id'
net.suid <- createNetworkFromDataFrames(net.cf, net.edges, title=paste(paste(nodenames,  collapse=" to "), 1+length(getNetworkList())), collection = "Bromodomain Interactions")
#
# Try plotting with all the edges using the new methods
connected.path$Alt.Weight <- connected.path$Weight
net2.edges <- rbind(net.gpe, ptm.cccn, connected.path)
colnames(net2.edges)[1:3]<-c("source", "target", "interaction")
net.suid <- createNetworkFromDataFrames(net.cf, net2.edges, title=paste(paste(nodenames,  collapse=" to "), 1+length(getNetworkList())), collection = "Bromodomain Interactions")
#```
#This shows all the edges between interacting proteins when there are multiple interactions. 

# Apply style from New_RCy3_Functions

nodeDprops.RCy32(net.cf)
ratioProps.RCy32(net.cf, "Median")
edgeDprops.RCy32()
#Change the appearance of protein-PTM edges
allwedges <- getAllEdges()
prpepedges <- allwedges[grep("peptide", allwedges)]
setEdgeColorBypass(prpepedges, new.colors="#808080")
setEdgeLineWidthBypass(prpepedges, new.widths = 2)
setBackgroundColorDefault(col2hex('grey80'))
new.style <- "THE Test Style"
copyVisualStyle("default", new.style)
setVisualStyle(new.style)

#------------------------------------------------------------
# Note: Skip the above and just open MetaNodeTestNetwork.cys file
#
genes <- net.cf[grep("gene", net.cf$Node.ID), "id"]
selectNodes(net.cf[grep("gene", net.cf$Node.ID), "id"], by.col="id")
# Note: default is SUID
ptms <- net.cf[grep("peptide", net.cf$Node.ID), "id"]
selectNodes(ptms, by='id', pre=F)
nodedata <- getTableColumns("node")
edgedata <- getTableColumns("edge")
# Note that SUID is the first column. Selected nodes are T/F in 'selected'
deltacatnodes <- net.cf[grep("CTNND1", net.cf$Gene.Name), "id"]
selectNodes(deltacatnodes, by.col="id", preserve=F)
#------------------------------------------------------------
# Experiment with loop to collapse all groups into gene

for(i in 1:length(genes)) {
    print(genes[i])
    selectNodes(net.cf[grep(genes[i], net.cf$Gene.Name), "id"], by.col="id", preserve=F)
    createGroup(genes[i])
    collapseGroup(genes[i])
}
groups.1 <- listGroups()
# Wow!
# $groups
# [1] 284 418 312 354 451
metanodedata <- getTableColumns("node")
edgedata <- getTableColumns("edge")
# The quantities appear to be merged in some way, not added. The strings are conactenated.
# All the groups now appear as new networks in Cytoscape.
getGroupInfo("CTTN")
# Get all group's info
group.info <- list()
group.info <- lapply(listGroups()$groups, getGroupInfo)
# Okay, now we try expanding...
expandGroup("CTTN")
getGroupInfo("CTTN")
# Now returns error
RCy3::commandsPOST, HTTP Error Code: 500
url=http://localhost:1234/v1/commands/group/get
body={
    "node": "SUID:74",
    "network": "SUID:52" 
}
Error in commandsPOST(paste0("group get node=\"SUID:", group.suid, "\"",  : 
                                 java.lang.NullPointerException
# Try collapsing again
# Examine groups
group.info.2 <- lapply(listGroups()$groups, getGroupInfo)
# ERROR:
RCy3::commandsPOST, HTTP Error Code: 500
url=http://localhost:1234/v1/commands/group/get
body={
    "node": "SUID:",
    "network": "SUID:52" 
}
Error in commandsPOST(paste0("group get node=\"SUID:", group.suid, "\"",  : 
                                 java.lang.NullPointerException

listGroups()
# $groups
# [1] 284 418 312 354 451                           
# Same groups as above. 354 is CTTN
collapseGroup(354)
# This didn't work:
RCy3::commandsPOST, HTTP Error Code: 500
url=http://localhost:1234/v1/commands/group/collapse
body={
    "groupList": "354",
    "network": "SUID:52" 
}
Error in commandsPOST(paste0("group collapse groupList=\"", group.list,  : 
                                 java.lang.NullPointerException
#
collapseGroup("CTTN")
# Now Works!
# Try to do them all
expandGroup("ALK")
collapseGroup("ALK")
expandGroup(c("ALK", "NPM1"))
collapseGroup(c("ALK", "NPM1"))
expandGroup(genes)
collapseGroup(genes)
# Now working! 
#------------------------------------------------------------
# Note: three times during testing Cytoscape Hangs!  "Unsilencing event source: CTTN to ALK 1 default node"
#Terminal shows many errors in this form:
#   at org.eclipse.jetty.util.thread.QueuedThreadPool$3.run(QueuedThreadPool.java:543)
#   at java.lang.Thread.run(Thread.java:748)
# java.lang.NullPointerException 
#... 
