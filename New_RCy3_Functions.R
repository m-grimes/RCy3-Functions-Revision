# New versions of MGRCyFunctions.R linking RCy3 to Cytoscape
# author: Mark Grimes
# This is the best place to report issues:
# https://github.com/cytoscape/RCy3/issues
# Thanks!    - Alex
# Swagger: http://localhost:1234/v1/swaggerUI/swagger-ui/index.html?url=http://localhost:1234/v1/swagger.json#/
# ------------------------------------------------------------------------
# To UpDATE
# source("http://www.bioconductor.org/biocLite.R")
# biocLite("RCy3") 
library(devtools)
library(RCy3)
library(plyr)
library(dplyr)
library(RColorBrewer)
library(gplots)
library(igraph)

# There are several resources available:
browseVignettes('RCy3')
# https://github.com/cytoscape/RCy3/wiki/Upgrading-Existing-Scripts
# 

# Note: Start Cytoscape
# ------------------------------------------------------------------------
cytoscapePing ()
cytoscapeVersionInfo ()

# ------------------------------------------------------------------------
# Not ready yet:
# installApp('STRING')
# Cytoscape version 3.7 or greater is required. You are currently working with version 3.6.
sessionInfo()
#
# ------------------------------------------------------------------------
# First task: convert functions to display CFN and CCCN for https://cynetworkbrowser.umt.edu
# Note: For  For createNetworkFromDataFrames {RCy3}, in the node file, the first column must be named "id". 
# In edge files, the first two columns must be named "source" and "target"; and "edgeType" is no longer used, replaced with "interaction". For converting files we may wish to use this function:

edgeType.to.interaction <- function(net.edges) {
    names(net.edges)[1:2] <- c("source", "target")
    if (any(grepl("edgeType", names(net.edges)))) {
    names(net.edges)[grep("edgeType", names(net.edges))] <- "interaction"}
    return(net.edges) }
# 
# ✓
# Function to collapse multiple edges into directed and undirected merged edges. This is the old version that simply converts the column names at the end using the above function.
mergeEdges.dir <- function(edgefile) {
    # define directional and non-directional edges
    directed <- c("pp", "controls-phosphorylation-of", "controls-expression-of", "controls-transport-of",  "controls-state-change-of")
    undirected <- c("Physical interactions", "BioPlex", "in-complex-with",  'experiments',  'database',   "Pathway", "Predicted", "Genetic interactions", "correlation", "negative correlation", "positive correlation",  'combined_score', "merged" , "intersect", "peptide", 'homology', "Shared protein domains")
    # check for nodes in reversed orientation for undirected edges
    undir.edges <- edgefile[-(edgefile$edgeType %in% directed),]
    # NEW: simply sort/order the nodes
    # this works:
    # for(i in 1:dim(undir.edges)[1])	{
    #	undir.edges[i, 1:2] <- sort(as.vector(undir.edges[i,1:2]))
    #	}
    # Better: Working non-loop:
    undir.edges[, 1:2] <- t(apply(undir.edges[,1:2], 1, function (x) sort(x)))
    # merge undirected edges
    undir.merged <- ddply(undir.edges, .(Gene.1, Gene.2), numcolwise(sum), na.rm=TRUE)	
    undir.merged$edgeType <- dlply(undir.edges, .(Gene.1, Gene.2), function(x) paste0(x[, "edgeType"]))	
    undir.merged$edgeType <- sapply(undir.merged$edgeType, function(x) paste(c(x), collapse=", "))
    # undir.merged$Directed <- FALSE	
    # merge directed edges
    dir.edges <- edgefile[edgefile$edgeType %in% directed,]
    dir.merged <- ddply(dir.edges, .(Gene.1, Gene.2), numcolwise(sum), na.rm=TRUE)
    dir.merged$edgeType <- dlply(dir.edges, .(Gene.1, Gene.2), function(x) paste0(x[, "edgeType"]))	
    dir.merged$edgeType <- sapply(dir.merged$edgeType, function(x) paste(c(x), collapse=", "))
    # dir.merged$Directed <- TRUE
    edgefile.merged <- rbind(dir.merged, undir.merged)
    # Remove auto-phosphorylation loops
    edgefile.merged <- remove.autophos(edgefile.merged)
    # Convert to new RCy3.2 format
    edgefile.merged <- edgeType.to.interaction(edgefile.merged)
    return(edgefile.merged)
} 
# ✓
# New version to use the updated RCy3 2.x.x names
mergeEdges.RCy32 <- function(edgefile) {
    # define directional and non-directional edges
    directed <- c("pp", "controls-phosphorylation-of", "controls-expression-of", "controls-transport-of",  "controls-state-change-of")
    undirected <- c("Physical interactions", "BioPlex", "in-complex-with",  'experiments',  'database',   "Pathway", "Predicted", "Genetic interactions", "correlation", "negative correlation", "positive correlation",  'combined_score', "merged" , "intersect", "peptide", 'homology', "Shared protein domains")
    # check for nodes in reversed orientation for undirected edges
    undir.edges <- edgefile[-(edgefile$interaction %in% directed),]
    # NEW: simply sort/order the nodes
    # this works:
    # for(i in 1:dim(undir.edges)[1])	{
    #	undir.edges[i, 1:2] <- sort(as.vector(undir.edges[i,1:2]))
    #	}
    # Better: Working non-loop:
    undir.edges[, 1:2] <- t(apply(undir.edges[,1:2], 1, function (x) sort(x)))
    # merge undirected edges
    undir.merged <- ddply(undir.edges, .(source, target), numcolwise(sum), na.rm=TRUE)	
    undir.merged$interaction <- dlply(undir.edges, .(source, target), function(x) paste0(x[, "interaction"]))	
    undir.merged$interaction <- sapply(undir.merged$interaction, function(x) paste(c(x), collapse=", "))
    # undir.merged$Directed <- FALSE	
    # merge directed edges
    dir.edges <- edgefile[edgefile$interaction %in% directed,]
    dir.merged <- ddply(dir.edges, .(source, target), numcolwise(sum), na.rm=TRUE)
    dir.merged$interaction <- dlply(dir.edges, .(source, target), function(x) paste0(x[, "interaction"]))	
    dir.merged$interaction <- sapply(dir.merged$interaction, function(x) paste(c(x), collapse=", "))
    # dir.merged$Directed <- TRUE
    edgefile.merged <- rbind(dir.merged, undir.merged)
    # Remove auto-phosphorylation loops
    edgefile.merged <- remove.autophos(edgefile.merged)
    return(edgefile.merged)
} 
# ✓


# "Visual style setting methods have been updated to more consistently match Cytoscape's terminology."
# This function sets edge default settings. Use: edgeDprops.RCy32()
edgeDprops.RCy32 <- function() {
    setEdgeLineWidthDefault (3)
    setEdgeColorDefault ( "#FFFFFF")  # white
    setEdgeSelectionColorDefault ( "#FF69B4")  # hotpink
    edgecolors <- col2hex(c("red", "red", "magenta", "violet", "purple",  "green", "green2", "green3",  "aquamarine2", "cyan", "turquoise2", "cyan2", "lightseagreen", "gold",  "blue", "yellow", "slategrey", "darkslategrey", "grey", "black", "orange", "orange2"))
    edgecolorsplus <- col2hex(c("deeppink", "red", "red", "magenta", "violet", "purple",  "green", "green2", "green3",  "aquamarine2", "cyan", "turquoise2", "cyan2", "lightseagreen", "gold",  "blue", "yellow", "slategrey", "darkslategrey", "grey", "black", "orange", "orange2", "orangered2"))
    #  red; turquois; green; magenta; blue; violet; green;  bluegreen; black; gray; turquoiseblue; orange 
    edgeTypes <- c("pp", "controls-phosphorylation-of", "controls-expression-of", "controls-transport-of",  "controls-state-change-of", "Physical interactions", "BioPlex", "in-complex-with",  'experiments',  'database',   "Pathway", "Predicted", "Genetic interactions", "correlation", "negative correlation", "positive correlation",  'combined_score', "merged" , "intersect", "peptide", 'homology', "Shared protein domains") 
    # 22 edgeTypes            
    myarrows <- c ('Arrow', 'Arrow', 'Arrow', 'Arrow', "Arrow", 'None', 'None', 'None', 'None', 'None', 'None', 'None', 'None', 'None', 'None', 'None', 'None', 'None', 'None', 'None', 'None', 'None')
    setEdgeTargetArrowMapping( 'interaction', edgeTypes, myarrows, default.shape='None')  
    matchArrowColorToEdge('TRUE')
    setEdgeColorMapping( 'interaction', edgeTypes, edgecolors, 'd', default.color="#FFFFFF")
    # A grey background helps highlight some of the edges
    setBackgroundColorDefault("#949494") # grey 58
}
# ✓
# The following function is now working properly, though it is not clear why the Mapping fuction doesn't work. 
# log version looks good
setEdgeWidths.RCy32 <- function (edgefile, factor=1.2, log=TRUE)	{
    if (log==FALSE) {
        line.widths <- factor*abs(as.numeric(edgefile$Weight))
        }
    if (log==TRUE) {
        line.widths <- log(abs(edgefile$Weight)) + factor - min(log(abs(edgefile$Weight)))  
        }
    edgevalues <- getTableColumns('edge',c('name','Weight'))
    edgevalues['Weight'] <- line.widths
    setEdgeLineWidthBypass(edgevalues[['name']], edgevalues[['Weight']])
    # Not working:
    # setEdgeLineWidthMapping("Weight", table.column.values=edgefile$Weight, mapping.type = "c", widths=line.widths, default.width=1.2) 
    # for testing:
    # return(data.frame(edgefile$Weight, line.widths))
}
# ✓

# This function sets node default settings
nodeDprops.RCy32 <- function(nodefile) {
    #setBackgroundColorDefault("#949494") # grey 58
    setNodeShapeDefault( "ELLIPSE")
    setNodeColorDefault( '#F0FFFF') # azure1
    setNodeSizeDefault( 100) # for grey non-data nodes
    setNodeFontSizeDefault( 22)
    setNodeLabelColorDefault( '#000000')  # black
    setNodeBorderWidthDefault( 1.8)
    setNodeBorderColorDefault( '#888888')  # gray 
    molclasses <- c("unknown", "receptor tyrosine kinase",  "SH2 protein", "SH2-SH3 protein", "SH3 protein", "tyrosine kinase",  "SRC-family kinase",   "kinase", "phosphatase", "transcription factor", "RNA binding protein")
    #  NOTE getNodeShapes(cy) returns node shapes in random order!  Define manually 
    #	*12 for RCy2; 9 for RCy3
    # there are now 24 nodeType classes
    nodeshapes <- c("ELLIPSE","ROUND_RECTANGLE", "VEE", "VEE", "TRIANGLE", "HEXAGON", "DIAMOND", "OCTAGON", "OCTAGON", "PARALLELOGRAM", "RECTANGLE")
    setNodeSelectionColorDefault(  "#CC00FF") 
    setNodeShapeMapping ("nodeType", molclasses, nodeshapes, default.shape="ELLIPSE")
    setNodeBorderWidthMapping("nodeType", c("deacetylase","acetyltransferase","demethylase","methyltransferase","membrane protein", "receptor tyrosine kinase", "G protein-coupled receptor", "SRC-family kinase", "tyrosine kinase", "kinase", "phosphatase"), widths=c(4,12,4,12,8,16,16,12,12,12,14), 'd',default.width=4)
    cf<-nodefile
    if (length(cf[grep("SH2", cf$Domains), 1])>0 & !all(grep("SH2", cf$Domains) %in% which(cf$nodeType %in% molclasses))) {
        setNodeShapeBypass(cf[grep("SH2", cf$Domains) %w/o% which(cf$nodeType %in% molclasses), 1], nodeshapes[3])} 
    if (length(cf[grep("RNA", cf$nodeType), 1])>0) {
        setNodeShapeBypass(cf[grep("RNA", cf$nodeType), 1], nodeshapes[11])}
    if (length(cf[grep("transcription", cf$nodeType), 1])>0) {
        setNodeShapeBypass(cf[grep("transcription", cf$nodeType), 1], nodeshapes[10])}
    if (length(cf[grep("acetyl", cf$nodeType), 1])>0) {
        setNodeBorderColorBypass(cf[grep("acetyl", cf$nodeType), 1], "#FF8C00")} # darkorange
    if (length(cf[grep("methyl", cf$nodeType), 1])>0) {
        setNodeBorderColorBypass(cf[grep("methyl", cf$nodeType), 1], "#005CE6")} # blue
    if (length(cf[grep("membrane", cf$nodeType), 1])>0) {
        setNodeBorderColorBypass(cf[grep("membrane", cf$nodeType), 1], "#6600CC") # purple
        setNodeShapeBypass(cf[grep("membrane", cf$nodeType), 1], nodeshapes[2])} 
    if (length(cf[grep("kinase", cf$nodeType), 1])>0) {
        setNodeBorderColorBypass(cf[grep("kinase", cf$nodeType), 1], "#EE0000")} # red2
    if (length(cf[grep("phosphatase", cf$nodeType), 1])>0) {
        setNodeBorderColorBypass(cf[grep("phosphatase", cf$nodeType), 1], "#FFEC8B")} # lightgoldenrod1
    if (length(cf[grep("receptor", cf$nodeType), 1])>0) {
        setNodeBorderColorBypass(cf[grep("receptor", cf$nodeType), 1], "#BF3EFF") # darkorchid1
        setNodeShapeBypass(cf[grep("receptor", cf$nodeType), 1], nodeshapes[2])} 
    if (length(cf[grep("TM", cf$nodeType), 1])>0) {
        setNodeBorderColorBypass(cf[grep("TM", cf$Domains), 1], "#6600CC") # purple
        setNodeShapeBypass(cf[grep("TM", cf$Domains), 1], nodeshapes[2])}   
}
# ✓

# Ratio props
ratioProps.RCy32 <- function (nodefile, plotcol="Total") {
    if(!(plotcol %in% getTableColumnNames('node'))){
        print (getTableColumnNames('node'))
        cat("\n","\n","\t", "Which attribute will set node size and color?")
        plotcol <- as.character(readLines(con = stdin(), n = 1))
    }
    limits <- range(nodefile[, plotcol])
    node.sizes     = c (135, 130, 108, 75, 35, 75, 108, 130, 135)
    #	RATIO is plotted
    #	Blue is negative: Yellow positive, Green in middle
    #		
    size.control.points = c (-100.0, -15.0, -5.0, 0.0, 5.0, 15.0, 100.0)
    color.control.points = c (-100.0, -10.0, -5.0, -2.25, 0.0, 2.25, 5.0, 10.0, 100.0)
    if(limits[1] < min(size.control.points)) {
        size.control.points = c (limits[1], -15.0, -5.0, 0.0, 5.0, 15.0, 100.0)
        color.control.points = c (limits[1]-1, -10.0, -5.0, -2.25, 0.0, 2.25, 5.0, 10.0, 100.0)
    }
    if(limits[2] > max(size.control.points)) {
        size.control.points = c (limits[1], -15.0, -5.0, 0.0, 5.0, 15.0, limits[2])
        color.control.points = c (limits[1]-1, -10.0, -5.0, -2.25, 0.0, 2.25, 5.0, 10.0, limits[2]+1)
    }
    ratio.colors = c ('#0099FF', '#007FFF','#00BFFF', '#00CCFF', '#00FFFF', '#00EE00', '#FFFF7E', '#FFFF00', '#FFE600', '#FFD700', '#FFCC00')
    setNodeColorMapping (names(nodefile[plotcol]), color.control.points, ratio.colors, 'c')
    lockNodeDimensions('TRUE')
    setNodeSizeMapping (names(nodefile[plotcol]), size.control.points, node.sizes, 'c')
    setNodeSelectionColorDefault ( "#CC00FF") 
}
# ✓
# Intensity props
intensityprops.RCy32 <- function (nodefile, plotcol="Total.Phosphorylation") {
    setVisualStyle ("default")
    # print (getTableColumnNames ())
    node.sizes     = c (135, 130, 108, 75, 35, 75, 108, 130, 135)
    Intensity.Values <- nodefile[, plotcol]  # set to intensity or normalized intensity
    maxint <- max(Intensity.Values, na.rm=TRUE) 
    minint <- min(Intensity.Values, na.rm=TRUE)
    icolors <- c('#0099FF', '#007FFF','#00BFFF', '#00CCFF', '#00FFFF', '#FFFFFF', '#FFFF7E', '#FFFF00', '#FFE600', '#FFD700', '#FFCC00')
    #  Some Mappings to set the color and size depending on the values of intensity
    if (maxint>=abs(minint)) {
        color.control.points <- c(-(maxint+1), -(maxint/5), -(maxint/10), -(maxint*0.045), 0.0, (maxint*0.045), (maxint/10), (maxint/5), (maxint+1))
        setNodeColorMapping (names(nodefile[plotcol]),  color.control.points, icolors) 
        size.control.points = c (-(maxint+1), -(maxint*0.3), -(maxint/10), 0.0, (maxint/10), (maxint*0.3), (maxint+1))
    }
    if (maxint<abs(minint)) {
        color.control.points <- c((minint-1), (minint/5), (minint/10), (minint*0.045), 0.0, -(minint*0.045), -(minint/10), -(minint/5), -(minint-1))
        setNodeColorMapping (names(nodefile[plotcol]), color.control.points, icolors) 
        size.control.points = c ((minint-1), (minint*0.3), (minint/10), 0.0, -(minint/10), -(minint*0.3), abs(minint-1))
    }
    setNodeSizeMapping (names(nodefile[plotcol]), size.control.points, node.sizes)
}
# ✓







# Alex: Regarding your question about meta nodes: yes, there are grouping functions. Start by looking through the Swagger docs for 'group' functions: Help>Automation>CyREST Commands API (screenshot); e.g. see ?collapseGroup