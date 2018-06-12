# Node and edge visual properties test
     # Plus layout experiments
# Mark Grimes
# April 26, 2017
# Revised: June 11, 2018
# Tools:
# Cytoscape 3.6.1 (cytoscape.org) with cyREST app v 3.7.0, install using Apps->App Manager; and java 1.8 (including jdk and jre libraries; I'm currently running 1.8.0_171)
# RCy3
##### To install from bioconductor:
# source("http://www.bioconductor.org/biocLite.R")
# biocLite("RCy3") 
# Note that dependancies should also be installed
 # requires: library(RCy3) 
if(!"RCy3" %in% installed.packages()){
    source("https://bioconductor.org/biocLite.R")
    biocLite("RCy3")
}
library(RCy3)
#
library(gplots)
# This seems to help getting rid of unwanted class="factor"
options(stringsAsFactors=FALSE)
 # Turn on Cytoscape, preferrably using the terminal, then, in R:
 # ------------------------------------------------------------------------
 cytoscapePing ()
 cytoscapeVersionInfo ()
 # ------------------------------------------------------------------------
 # This establishes the connection. You should get: [1] "You are connected to Cytoscape!"
# Feature requests
# 1. Better layout contols:
     # Feature request: The ability to access all settings on layouts; 
     # e.g. to work only on selected nodes or choose the edge attribute
     # This worked previously in RCytoscape: setLayoutProperties(cy, "force-directed", list(selected_only=TRUE) )
     # Also: the ability to use Weight or Alt.Weight as layout spring attribute!
#
# 2. The ability use set_Rule functions, copyVisualStyle, and other functions without abolising all visual properties set by set_Direct functions. Examples of how this breaks are shown below. 
# 3. The ability to control metanodes with the "parent" attribute; e.g. expand and collapse.

# NODEKEY example
# cyPlot <- function (node.df, edge.df) # now part of RCy3

# Define node default property setting function

		
# Set up node key graph
	molclasses <- c("unknown", "receptor tyrosine kinase",  "SH2 protein", "SH3 protein", "tyrosine kinase",  "SRC-family kinase",   "kinase", "phosphatase", "transcription factor", "RNA binding protein")
	nodeTypes <- c(molclasses,"deacetylase","acetyltransferase","demethylase","methyltransferase","membrane protein", "receptor tyrosine kinase", "G protein-coupled receptor")
	nulledges <- data.frame(Gene.1="unknown", Gene.2="unknown", edgeType="unknown", Weight=0)
	nodetype.cf <- data.frame(Gene.Name=nodeTypes, nodeType=nodeTypes, size=130)
	nodetype.cf$Domains <- c("undefined", "TM", "SH2", "SH3", "tyrosine kinase", "tyrosine kinase", "kinase", "phosphatase", "zinc finger", "RNA binding", "undefined", "undefined","undefined","undefined", "TM", "TM", "TM")

# Note that sometimes the scale is compressed. The ability to control the tool panel would be useful!
     
# Edge key

interaction <- c("pp", "controls-phosphorylation-of", "controls-expression-of", "controls-transport-of",  "controls-state-change-of", "Physical interactions", "BioPlex", "in-complex-with",  'experiments',  'database',   "Pathway", "Predicted", "Genetic interactions", "correlation", "negative correlation", "positive correlation",  'combined_score', "merged" , "intersect", "peptide", 'homology', "Shared protein domains") 
	
edgecolors <- col2hex(c("red", "red", "magenta", "violet", "purple",  "green", "green2", "green3",  "aquamarine2", "cyan", "turquoise2", "cyan2", "lightseagreen", "gold",  "blue", "yellow", "slategrey", "darkslategrey", "grey", "black", "orange", "orange2"))

myarrows <- c ('Arrow', 'Arrow', 'Arrow', 'Arrow', "Arrow", 'No Arrow', 'No Arrow', 'No Arrow', 'No Arrow', 'No Arrow', 'No Arrow', 'No Arrow', 'No Arrow', 'No Arrow', 'No Arrow', 'No Arrow', 'No Arrow', 'No Arrow', 'No Arrow', 'No Arrow', 'No Arrow', 'No Arrow')
 
edge_key.tbl <- data.frame(interaction, colors=c("red", "red", "magenta", "violet", "purple",  "green", "green2", "green3",  "aquamarine2", "cyan", "turquoise2", "cyan2", "lightseagreen", "gold",  "blue", "yellow", "slategrey", "darkslategrey", "grey", "black", "orange", "orange2"), myarrows, Weight=c(10^-seq(0, 3, length.out = 4), seq(1, 11, 3), 10^-seq(0, 3, length.out = 4), 1, 0.75, -0.5, seq(0.5, 0.8, 0.1), 100, 1:2))
 
hist(edge_key.tbl$Weight, breaks=100, col="magenta", ylim=c(0, 50))
legend(5, 50, edgeTypes, cex=0.8, col= edgecolors,  lty=rep(1,24), lwd=4, bty="n");
# Bigger version (after expanding plot window vertically):
dev.new()
hist(edge_key.tbl$Weight, breaks=100, col="green", ylim=c(0, 50))
legend(5, 50, edgeTypes, cex=1.6, col= edgecolors,  lty=rep(1,24), lwd=8, bty="n");

# RCy3 and Cy test
# Define edge default property setting functions in New_RCy3_Functions.R


#edgefile <- data.frame(Gene.1=LETTERS[1:23], Gene.2=paste(LETTERS[1:23], 1, sep=""), edgeType=c("unknown", edgeTypes), Weight=c(exp(seq(-20, 5, 2)), 1:6, 10:13))
edgefile <- data.frame(source=LETTERS[1:22], target=paste(LETTERS[1:22], 1, sep=""))
edgefile <- cbind(edgefile, edge_key.tbl)
edgefile[20, 2] <- "T ack 100" # to represent PTM 
# The alternate weight attemps to drive negative correlations far away in certain layouts (force-directed; spring-embedded)
negs <- edgefile[edgefile $Weight<0,'Weight']
	frnegs <- exp (negs*20)
	# plot (negs,frnegs)
	edgefile $Alt.Weight <- edgefile $Weight
	edgefile[edgefile $Weight<0,'Alt.Weight'] <- frnegs	

nodefile <- data.frame(id=c(edgefile$source, edgefile$target), nodeType="unknown", Domains="unknown")
molclasses1 <- c("unknown", "receptor tyrosine kinase",  "SH2 protein", "SH2-SH3 protein", "SH3 protein", "tyrosine kinase",  "SRC-family kinase",   "kinase", "phosphatase", "transcription factor", "RNA binding protein")
molclasses2 <- c("deacetylase","acetyltransferase","demethylase","methyltransferase","membrane protein",  "G protein-coupled receptor", "RNA processing protein", "RNA splicing protein", "membrane protein", "acting-binding protein", "microtubule-associated protein")
nodefile[, "nodeType"] <- c(molclasses1, molclasses2, sample(c(molclasses1, molclasses2), 22))
nodefile$Domains <- sapply(nodefile$nodeType, function(x) unlist(strsplit(x, " "))[1])
nodefile[grep("membrane", nodefile$nodeType), "Domains"] <- "TM"
nodefile[grep("G protein", nodefile$nodeType), "Domains"] <- "TM"
nodefile[42, 2:3] <- "peptide"
ratios <- 10^(seq(-3, 4, length.out = 22))
nodefile$Total <- seq(440, 10, length.out=44)
nodefile$Ratio <- c(ratios, -ratios)
edgekey.suid <- createNetworkFromDataFrames(nodefile, edgefile, title=paste("Edge Key", 1+length(getNetworkList())), collection = "Keys")
# works
# Apply formatting function
edgeDprops.RCy32()
# Works: Edges should have different colors and directed edges, arrows.
# Apply edge width formatting fuction
setEdgeWidths.RCy32(edgefile=edgefile, factor=1.4, log=TRUE)
# Works: edge widths now proportional to Weight
# Now set node shapes and colors to indicate molecular classes:
nodeDprops.RCy32(nodefile)
# Now color and size nodes based on quantitative attributes
ratioProps.RCy32(nodefile, plotcol = "Ratio")
ratioProps.RCy32(nodefile, plotcol="Total")
intensityprops.RCy32(nodefile, plotcol="Total")
intensityprops.RCy32(nodefile, plotcol="Ratio")
# Note that the size and color gradients are slightly different for these two functions.
# Ratio differences are notcible at ratios > +/- 2.5; intensity differnces are more spread out depending on data.


	nodelayoutNetwork('force-directed defaultSpringCoefficient=0.00001 defaultSpringLength=50 defaultNodeMass=5')
	
	layoutNetwork("hierarchical")     
	
	# 
	 setDefaultBackgroundColor (edge.w, "#949494") # grey 58
	 # 
	#
	 setEdgeLabelRule (edge.w, 'edgeType')
	 # 
	 # 
	 # The following now works, but note what happens to the node key window - colors vanish.
	 	line.widths <- 4.5*abs(as.numeric(edgefile$Weight))
		setEdgeLineWidthRule(edge.w, "Weight", attribute.values=as.character(edgefile$Weight), line.widths, default.width=1.2)
	# Try a different function (after restarting Cytoscape and reloading graph)
	setEdgeLineWidthDirect(edge.w, getAllEdges(edge.w), line.widths)	
	# Works! But
	edgeDprops(edge.w)
	# Erases the edge widths, now all the same.
	# Note also
		redraw(edge.w)
		# Undoes the edge widths also
	data.frame(getAllEdges(edge.w), line.widths)
		 
# Note that any function that sets global properties removes the specifically set colors and shapes with setDirect functions in the last part of the nodeDprops. This includes setting e.g., background in Cytoscape itself. Or even changing the style name. Or 	redraw(windowobj)

# If the above is run in the same Cy session, the node colors are removed in the first window after the edge properties are set in the second window. 

# Function to get cy edge names from edgefile (to ensure the Weights match up)

getCyEdgeNames <- function(edgefile) {
	cyedges <- mapply(paste, edgefile $Gene.1, " (", edgefile $edgeType, ") ", edgefile $Gene.2, sep="")
	return(cyedges)
} 

# Functions to set edge widths using setEdgeLineWidthDirect

setEdgeWidths <- function (windowobj, factor=1)	{
	alledges <- getAllEdges(windowobj)
	edgeweights <- sapply(alledges, function(x) getEdgeAttribute(windowobj, x, "Weight"))
	line.widths <- factor*abs(as.numeric(edgeweights))
	setEdgeLineWidthDirect(windowobj, alledges, line.widths)	
}

setEdgeWidths.log <- function (windowobj, factor=1)	{
	alledges <- getAllEdges(windowobj)
	edgeweights <- sapply(alledges, function(x) getEdgeAttribute(windowobj, x, "Weight"))
	line.widths <- as.numeric(log(abs(edgeweights)) + factor - min(log(abs(edgeweights))))
	setEdgeLineWidthDirect(windowobj, alledges, unname(line.widths))	
}
	
	setEdgeWidths(edge.w, factor=3)
	setEdgeWidths.log(edge.w, factor=2)
# Note: both work directly. But visual properties are lost if style is saved under another name. 
	
# LAYOUT example 	
	# I'm going to create a real-world exmaple from scratch using the following to create a data frame
	# options(useFancyQuotes = FALSE); 
	# for string:cat(sQuote(edgefile[,1]), sep=", ");
	# for numerical data: cat(edgefile[,4], sep=", ")

# Funtion to graph node size and color
	ratioprops <- 	function (windowname, cytoscape.file, plotcol="Total") {
	     if(!(plotcol %in% getNodeAttributeNames(windowname))){
	          print (getNodeAttributeNames (windowname))
	          cat("\n","\n","\t", "Which attribute will set node size and color?")
	          plotcol <- as.character(readLines(con = stdin(), n = 1))}
	     setVisualStyle (windowname, "default")
	     limits <- range(cytoscape.file[, plotcol])
	     #
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
	     # displayGraph (windowname)
	     setNodeColorRule (windowname, names(cytoscape.file[plotcol]), color.control.points, ratio.colors, mode='interpolate')
	     setNodeSizeRule (windowname, names(cytoscape.file[plotcol]), size.control.points, node.sizes, mode='interpolate')
	     setDefaultNodeSelectionColor (windowname,  "#CC00FF") 
	     new.style.name = paste(plotcol, noquote("_ratio"),  sep="", collapse=NULL) 
	     # if (!(new.style.name %in% getVisualStyleNames(cy)))	{copyVisualStyle(windowname, 'default', new.style.name)} # broken
	     # setVisualStyle(windowname, new.style.name)
	     # return(windowname)
	} # end ratioprops
# Note that the commented out copyVisualStyle function is broken - another feature request!
	
	net.edges <- data.frame(Gene.1=c('ALK', 'ALK', 'ALK', 'CTNND1', 'CTNND1', 'CTNND1', 'CTNND1', 'CTNND1', 'CTNND1', 'CTNND1', 'CTTN', 'CTTN', 'CTTN', 'IRS1', 'IRS1', 'IRS1', 'ALK p Y1096', 'CTNND1 p Y193', 'CTNND1 p Y193', 'CTNND1 p Y228', 'CTNND1 p Y904', 'CTNND1 p Y217', 'CTNND1 p Y241', 'CTNND1 p Y248', 'ALK p Y1078', 'ALK p Y1096', 'ALK p Y1586', 'IRS1 p Y941', 'ALK', 'CTNND1', 'CTNND1', 'CTTN', 'IRS1'), Gene.2=c('ALK p Y1078', 'ALK p Y1096', 'ALK p Y1586', 'CTNND1 p Y193', 'CTNND1 p Y217', 'CTNND1 p Y228', 'CTNND1 p Y241', 'CTNND1 p Y248', 'CTNND1 p Y302', 'CTNND1 p Y904', 'CTTN p Y154', 'CTTN p Y162', 'CTTN p Y334', 'IRS1 p Y632', 'IRS1 p Y941', 'IRS1 p Y989', 'ALK p Y1586', 'CTNND1 p Y228', 'CTNND1 p Y302', 'CTNND1 p Y302', 'CTTN p Y154', 'CTTN p Y162', 'CTTN p Y162', 'CTTN p Y334', 'IRS1 p Y632', 'IRS1 p Y989', 'IRS1 p Y989', 'IRS1 p Y989', 'IRS1', 'CTTN', 'IRS1', 'NPM1', 'NPM1'), edgeType=c('peptide', 'peptide', 'peptide', 'peptide', 'peptide', 'peptide', 'peptide', 'peptide', 'peptide', 'peptide', 'peptide', 'peptide', 'peptide', 'peptide', 'peptide', 'peptide', 'positive correlation', 'positive correlation', 'positive correlation', 'positive correlation', 'positive correlation', 'positive correlation', 'positive correlation', 'positive correlation', 'positive correlation', 'positive correlation', 'positive correlation', 'positive correlation', 'Physical interactions, experiments', 'controls-state-change-of, controls-transport-of, Pathway, Physical interactions, experiments', 'Predicted', 'experiments, Physical interactions', 'experiments, Physical interactions'), Weight=c(100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 0.8060606, 0.7575758, 0.7454545, 0.9393939, 0.8949096, 0.7329699, 0.7553845, 0.7866191, 0.775, 0.6969697, 0.7818182, 0.8424242, 0.6123365, 2.115272, 0.002461723, 0.3354451, 0.5661711), Alt.Weight=c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.8060606, 0.7575758, 0.7454545, 0.9393939, 0.8949096, 0.7329699, 0.7553845, 0.7866191, 0.775, 0.6969697, 0.7818182, 0.8424242, 0.6123365, 2.115272, 0.002461723, 0.3354451, 0.5661711))
	
	net.cf <- data.frame(Node=c('ALK', 'ALK p Y1078', 'ALK p Y1096', 'ALK p Y1586', 'CTNND1', 'CTNND1 p Y193', 'CTNND1 p Y217', 'CTNND1 p Y228', 'CTNND1 p Y241', 'CTNND1 p Y248', 'CTNND1 p Y302', 'CTNND1 p Y904', 'CTTN', 'CTTN p Y154', 'CTTN p Y162', 'CTTN p Y334', 'IRS1', 'IRS1 p Y632', 'IRS1 p Y941', 'IRS1 p Y989', 'NPM1'), Gene.Name=c('ALK', 'ALK', 'ALK', 'ALK', 'CTNND1', 'CTNND1', 'CTNND1', 'CTNND1', 'CTNND1', 'CTNND1', 'CTNND1', 'CTNND1', 'CTTN', 'CTTN', 'CTTN', 'CTTN', 'IRS1', 'IRS1', 'IRS1', 'IRS1', 'NPM1'), nodeType=c('receptor tyrosine kinase', 'receptor tyrosine kinase', 'receptor tyrosine kinase', 'receptor tyrosine kinase', 'undefined', 'undefined', 'undefined', 'undefined', 'undefined', 'undefined', 'undefined', 'undefined', 'SH3 protein', 'SH3 protein', 'SH3 protein', 'SH3 protein', 'undefined', 'undefined', 'undefined', 'undefined', 'RNA processing protein'), Total=c(-209.4885, -11.15512, -7.712652, -9.917654, -783.7109, 0, -16.63529, 0, -95.21879, -35.03751, 0, -74.77153, -172.7889, -67.19445, -14.73001, -55.37009, -401.1118, -6.4122, -6.388159, -5.171309, 1.580289), parent=c('', 'ALK', 'ALK', 'ALK', '', 'CTNND1', 'CTNND1', 'CTNND1', 'CTNND1', 'CTNND1', 'CTNND1', 'CTNND1', '', 'CTTN', 'CTTN', 'CTTN', '', 'IRS1', 'IRS1', 'IRS1', ''), Node.ID=c('gene', 'peptide', 'peptide', 'peptide', 'gene', 'peptide', 'peptide', 'peptide', 'peptide', 'peptide', 'peptide', 'peptide', 'gene', 'peptide', 'peptide', 'peptide', 'gene', 'peptide', 'peptide', 'peptide', 'gene'))
	
	net.g <- cyPlot(net.cf, net.edges)
	net.w <- CytoscapeWindow(paste((c("ALK", "CTTN")),  collapse=" to "), net.g)
	displayGraph(net.w)
	layoutNetwork(net.w, "genemania-force-directed")
	#setEdgeLineWidthRule(net.w, edge.attribute.name="Weight", attribute.values=net.edges$Weight, line.widths=net.edges$Weight)
	edgeDprops(net.w)
	ratioprops(net.w, net.cf, "Total")  
	showGraphicsDetails(net.w, TRUE) 
	nodeDprops.new(net.w, net.cf)
	setEdgeWidths.log(net.w, factor=1.2)
	
     # Again the graph is compressed. 
	layoutNetwork(net.w, layout.name="kamada-kawai")
	# A bit better, not great.
	
# Layout property access
	# I would like the parent or gene nodes to be driving the layouts and the peptides to surround them. For this reason I set the edge weight from gene to peptide to be large. But my attempts to tweak existing layouts have not accomplished this goal.  A workaround could be 
	selectNodes(net.w, net.cf[which(net.cf$Node.ID=="gene"), "Node"]) 
	# The following use to work in RCytoscape. But NOT in RCy3
	setLayoutProperties(cy, "force-directed", list(selected_only=TRUE) )
	# selected_only is not a property in layout force-directed	
	setLayoutProperties(cy, "kamada-kawai", list(selected_only=TRUE) )
	# selected_only is not a property in layout kamada-kawai	
	layouts <- getLayoutNames(cy)
	getLayoutPropertyNames(cy, "genemania-force-directed")
	getLayoutPropertyNames(cy, "kamada-kawai")
	# Neither of these specify the quantitative edge attribute to drive the layout (very important!)
	# The following are my unsuccessful tweaks
	defalutgmvalues <- 	data.frame(name=getLayoutPropertyNames(cy, "genemania-force-directed"), value=getLayoutPropertyValue(cy, "genemania-force-directed", getLayoutPropertyNames(cy, "genemania-force-directed")))
	# Has to be manual without the ability to layout with selected_only=TRUE
	# Manually: genemania force directed with weight works okay, after manually resizing
	setLayoutProperties(cy, "genemania-force-directed", list(numIterations=100, defaultSpringCoefficient=0.1, defaultSpringLength=50, minNodeMass=0.001, maxNodeMass=2000, midpointEdges=250, curveSteepness=7.0e-03, isDeterministic=1, singlePartition=0, ignoreHiddenElements=1))
	layoutNetwork(net.w, layout.name= "genemania-force-directed")
	# another unsuccesful tweak
	defalutkkvalues <- 	data.frame(name=getLayoutPropertyNames(cy, "kamada-kawai"), value=getLayoutPropertyValue(cy, "kamada-kawai", getLayoutPropertyNames(cy, "kamada-kawai")))
	setLayoutProperties (cy, "kamada-kawai", list (m_averageIterationsPerNode=50, m_nodeDistanceStrengthConstant=8000, m_nodeDistanceRestLengthConstant=400, m_disconnectedNodeDistanceSpringStrength=0.0001, m_disconnectedNodeDistanceSpringRestLength=2000, m_anticollisionSpringStrength=0, m_layoutPass=50, singlePartition=10, unweighted=1, randomize=10)) 
	layoutNetwork(net.w, layout.name="kamada-kawai")
	# Not happy with this. 
# NOT RUN:	
> sessionInfo()
	R version 3.3.3 (2017-03-06)
	Platform: x86_64-apple-darwin13.4.0 (64-bit)
	Running under: macOS Sierra 10.12.4
	
	locale:
	     [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
	
	attached base packages:
	     [1] parallel  stats4    stats     graphics  grDevices utils     datasets  methods   base     
	
	other attached packages:
	     [1] RCy3_1.5.2            rgl_0.98.1            ndexr_0.99.6          tidyr_0.6.1           jsonlite_1.3         
	[6] data.table_1.10.4     glmnet_2.0-5          foreach_1.4.3         KEGGlincs_1.0.0       hgu133a.db_3.2.3     
	[11] KOdata_1.0.0          STRINGdb_1.14.0       preprocessCore_1.36.0 RefNet_1.10.1         shiny_1.0.1          
	[16] RCurl_1.95-4.8        bitops_1.0-6          AnnotationHub_2.6.5   PSICQUIC_1.12.0       httr_1.2.1           
	[21] biomaRt_2.30.0        Biostrings_2.42.1     XVector_0.14.1        GOstats_2.40.0        graph_1.52.0         
	[26] Category_2.40.0       GO.db_3.4.0           org.Hs.eg.db_3.4.0    AnnotationDbi_1.36.2  IRanges_2.8.2        
	[31] S4Vectors_0.12.2      Biobase_2.34.0        BiocGenerics_0.20.0   PMA_1.0.9             plyr_1.8.4           
	[36] impute_1.48.0         genlasso_1.3          igraph_1.0.1          Matrix_1.2-8          MASS_7.3-45          
	[41] devtools_1.12.0       gplots_3.0.1          Rtsne_0.11            tsne_0.1-3            vegan_2.4-2          
	[46] permute_0.9-4         RUnit_0.4.31          cluster_2.0.6         Hmisc_4.0-2           ggplot2_2.2.1        
	[51] Formula_1.2-1         survival_2.41-3       lattice_0.20-35       stringr_1.2.0        
	
	loaded via a namespace (and not attached):
	     [1] colorspace_1.3-2              htmlTable_1.9                 base64enc_0.1-3              
	[4] hash_2.2.6                    interactiveDisplayBase_1.12.0 sqldf_0.4-10                 
	[7] codetools_0.2-15              splines_3.3.3                 knitr_1.15.1                 
	[10] annotate_1.52.1               png_0.1-7                     backports_1.0.5              
	[13] lazyeval_0.2.0                acepack_1.4.1                 htmltools_0.3.5              
	[16] tools_3.3.3                   gtable_0.2.0                  Rcpp_0.12.10                 
	[19] RJSONIO_1.3-0                 gdata_2.17.0                  nlme_3.1-131                 
	[22] iterators_1.0.8               proto_1.0.0                   mime_0.5                     
	[25] gtools_3.5.0                  XML_3.98-1.6                  zlibbioc_1.20.0              
	[28] scales_0.4.1                  BiocInstaller_1.24.0          RBGL_1.50.0                  
	[31] KEGGgraph_1.32.0              RColorBrewer_1.1-2            curl_2.4                     
	[34] yaml_2.1.14                   memoise_1.0.0                 gridExtra_2.2.1              
	[37] rpart_4.1-10                  latticeExtra_0.6-28           stringi_1.1.3                
	[40] RSQLite_1.1-2                 genefilter_1.56.0             plotrix_3.6-4                
	[43] checkmate_1.8.2               caTools_1.17.1                chron_2.3-50                 
	[46] htmlwidgets_0.8               GSEABase_1.36.0               AnnotationForge_1.16.1       
	[49] magrittr_1.5                  R6_2.2.0                      DBI_0.6-1                    
	[52] gsubfn_0.6-6                  foreign_0.8-67                withr_1.0.2                  
	[55] mgcv_1.8-17                   KEGGREST_1.14.1               nnet_7.3-12                  
	[58] tibble_1.3.0                  KernSmooth_2.23-15            grid_3.3.3                   
	[61] digest_0.6.12                 xtable_1.8-2                  httpuv_1.3.3                 
	[64] munsell_0.4.3                 
