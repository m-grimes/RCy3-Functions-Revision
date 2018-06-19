# Metanode experiments

# Mark Grimes
source("New_RCy3_Functions.R")
# To source from either computer
here <- getwd()
computer <- paste("", unlist(strsplit(here, "/"))[2], unlist(strsplit(here, "/"))[3], collapse = NULL, sep = "/")
directory <- "/Dropbox/_Work/R_/_RCy3\ future\ discussion/RCy3-Functions-Revision/"
file <- "New_RCy3_Functions.R"
dir.file <- "/Dropbox/_Work/R_/_RCy3\ future\ discussion/RCy3-Functions-Revision/New_RCy3_Functions.R"
source(paste(computer, directory, file, sep=""))
source(paste(computer, dir.file, sep=""))
ls(pat="RCy32")
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
# _________________________________________________________________________________
# Make a table of these data for vignette
options(useFancyQuotes = FALSE)
cat(dQuote(net.cf$id), sep=", ")
cat(dQuote(net.cf$parent), sep=", ")
cat(dQuote(net.cf$Node.ID), sep=", ")
cat(dQuote(net.edges$source), sep=", ")
cat(dQuote(net.edges$target), sep=", ")
cat(net.edges$Weight, sep=", ")

# _________________________________________________________________________________

# Define node data 
net.nodes <- c("ALK", "ALK p Y1078", "ALK p Y1096", "ALK p Y1586", "CTNND1", "CTNND1 p Y193", "CTNND1 p Y217", "CTNND1 p Y228", "CTNND1 p Y241", "CTNND1 p Y248", "CTNND1 p Y302", "CTNND1 p Y904", "CTTN", "CTTN ack K107", "CTTN ack K124", "CTTN ack K147", "CTTN ack K161", "CTTN ack K235", "CTTN ack K390", "CTTN ack K87", "CTTN p S113", "CTTN p S224", "CTTN p Y104", "CTTN p Y154", "CTTN p Y162", "CTTN p Y228", "CTTN p Y334", "CTTN p Y421", "IRS1", "IRS1 p Y632", "IRS1 p Y941", "IRS1 p Y989", "NPM1", "NPM1 ack K154", "NPM1 ack K223", "NPM1 p S214", "NPM1 p S218")
net.genes <- sapply(net.nodes,  function (x) unlist(strsplit(x, " ",  fixed=TRUE))[1])
parent <- c("", "ALK", "ALK", "ALK", "", "CTNND1", "CTNND1", "CTNND1", "CTNND1", "CTNND1", "CTNND1", "CTNND1", "", "CTTN", "CTTN", "CTTN", "CTTN", "CTTN", "CTTN", "CTTN", "CTTN", "CTTN", "CTTN", "CTTN", "CTTN", "CTTN", "CTTN", "CTTN", "", "IRS1", "IRS1", "IRS1", "", "NPM1", "NPM1", "NPM1", "NPM1")
nodeType <- c("protein", "modification", "modification", "modification", "protein", "modification", "modification", "modification", "modification", "modification", "modification", "modification", "protein", "modification", "modification", "modification", "modification", "modification", "modification", "modification", "modification", "modification", "modification", "modification", "modification", "modification", "modification", "modification", "protein", "modification", "modification", "modification", "protein", "modification", "modification", "modification", "modification")
netnodes.df <- data.frame(id=net.nodes, Gene.Name=net.genes, parent, nodeType)
# Define edge data
source.nodes <- c("ALK", "ALK", "ALK", "CTNND1", "CTNND1", "CTNND1", "CTNND1", "CTNND1", "CTNND1", "CTNND1", "CTTN", "CTTN", "CTTN", "CTTN", "CTTN", "CTTN", "CTTN", "CTTN", "CTTN", "CTTN", "CTTN", "CTTN", "CTTN", "CTTN", "CTTN", "IRS1", "IRS1", "IRS1", "NPM1", "NPM1", "NPM1", "NPM1", "ALK p Y1096", "CTNND1 p Y193", "CTNND1 p Y193", "CTNND1 p Y228", "CTNND1 p Y904", "CTNND1 p Y217", "CTNND1 p Y241", "CTNND1 p Y248", "ALK p Y1078", "ALK p Y1096", "ALK p Y1586", "IRS1 p Y941", "CTTN ack K147", "CTTN ack K107", "CTTN ack K235", "CTTN ack K87", "CTTN ack K147", "CTTN ack K124", "CTTN ack K147", "CTTN ack K235", "CTTN ack K161", "CTTN ack K390", "NPM1 ack K223", "NPM1 ack K154", "NPM1 ack K223", "ALK", "CTNND1", "CTNND1", "CTTN", "IRS1")
target.nodes <- c("ALK p Y1078", "ALK p Y1096", "ALK p Y1586", "CTNND1 p Y193", "CTNND1 p Y217", "CTNND1 p Y228", "CTNND1 p Y241", "CTNND1 p Y248", "CTNND1 p Y302", "CTNND1 p Y904", "CTTN ack K107", "CTTN ack K124", "CTTN ack K147", "CTTN ack K161", "CTTN ack K235", "CTTN ack K390", "CTTN ack K87", "CTTN p S113", "CTTN p S224", "CTTN p Y104", "CTTN p Y154", "CTTN p Y162", "CTTN p Y228", "CTTN p Y334", "CTTN p Y421", "IRS1 p Y632", "IRS1 p Y941", "IRS1 p Y989", "NPM1 ack K154", "NPM1 ack K223", "NPM1 p S214", "NPM1 p S218", "ALK p Y1586", "CTNND1 p Y228", "CTNND1 p Y302", "CTNND1 p Y302", "CTTN p Y154", "CTTN p Y162", "CTTN p Y162", "CTTN p Y334", "IRS1 p Y632", "IRS1 p Y989", "IRS1 p Y989", "IRS1 p Y989", "CTTN p S113", "CTTN p S224", "CTTN p S224", "CTTN p S224", "CTTN p Y104", "CTTN p Y228", "CTTN p Y228", "CTTN p Y228", "CTTN p Y421", "CTTN p Y421", "NPM1 p S214", "NPM1 p S218", "NPM1 p S218", "IRS1", "CTTN", "IRS1", "NPM1", "NPM1")
Weight <- c(100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 0.8060606, 0.7575758, 0.7454545, 0.9393939, 0.8949096, 0.7329699, 0.7553845, 0.7866191, 0.775, 0.6969697, 0.7818182, 0.8424242, -0.7714286, -0.8385965, -0.5017544, -0.7473684, -0.5252838, -0.9428571, -0.8285714, -0.6713287, -0.5508772, -0.9428571, -0.8857143, -0.6310881, -0.8285714, 0.6123365, 2.115272, 0.002461723, 0.3354451, 0.5661711)
netedges.df <- data.frame(source=source.nodes, target=target.nodes, Weight)
net.suid <- createNetworkFromDataFrames(netnodes.df, netedges.df, title=paste(paste("Group Nodes Test"), 1+length(getNetworkList())), collection = "RCy3 Vignettes")
layoutNetwork('force-directed defaultSpringCoefficient=0.00001 defaultSpringLength=50 defaultNodeMass=5')


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
# Test graphNetworkPath.RCy32
test.edges <- connectNodes.all(c("EGFR", "SMARCA4"), ig.graph=essl.cfn.ig, edgefile=essl.cfn)
test.result <- graphNetworkPath.RCy32(nodenames=c("EGFR", "SMARCA4"), path.edges=test.edges, ptmedgefile=essl.cccnplus.edges, datacolumns=geldexpts.nc, geneatts=essl.netatts, ptmcccnatts=essl.ptmcccnatts, edgeMerge=FALSE)

nodedata <- getTableColumns("node")
edgedata <- getTableColumns("edge")

# Create function to collapse all or selective nodes

collapse.CCCN.nodes <- function(nodenames=NULL) {
    nodedata <- getTableColumns("node", columns = c("id", "Gene.Name", "parent", "Node.ID"))
    genes <- nodedata[grep("gene", nodedata$Node.ID), "id"]
    if (length(nodenames) > 0) {
        genes <- nodenames
    } 
    sapply(genes, function(x) createGroup(x, nodes=nodedata[grep(x, nodedata$Gene.Name), "id"], nodes.by.col = "id"))
    collapseGroup(genes)
}
 
collapse.CCCN.nodes( c("CDK1",    "CRK",     "DDX5"))
expandGroup(c("CDK1",    "CRK",     "DDX5"))
deleteGroup()
collapse.CCCN.nodes()
expandGroup("all")
