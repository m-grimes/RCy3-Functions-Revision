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
# Ratio differences are notcible at ratios > +/- 2.5; intensity differences are more spread out depending on data.
# End of new stuff.
#------------------------------------------------------------

	layoutNetwork('force-directed defaultSpringCoefficient=0.00001 defaultSpringLength=50 defaultNodeMass=5')
	

	
	net.edges <- data.frame(source=c('ALK', 'ALK', 'ALK', 'CTNND1', 'CTNND1', 'CTNND1', 'CTNND1', 'CTNND1', 'CTNND1', 'CTNND1', 'CTTN', 'CTTN', 'CTTN', 'IRS1', 'IRS1', 'IRS1', 'ALK p Y1096', 'CTNND1 p Y193', 'CTNND1 p Y193', 'CTNND1 p Y228', 'CTNND1 p Y904', 'CTNND1 p Y217', 'CTNND1 p Y241', 'CTNND1 p Y248', 'ALK p Y1078', 'ALK p Y1096', 'ALK p Y1586', 'IRS1 p Y941', 'ALK', 'CTNND1', 'CTNND1', 'CTTN', 'IRS1'), target=c('ALK p Y1078', 'ALK p Y1096', 'ALK p Y1586', 'CTNND1 p Y193', 'CTNND1 p Y217', 'CTNND1 p Y228', 'CTNND1 p Y241', 'CTNND1 p Y248', 'CTNND1 p Y302', 'CTNND1 p Y904', 'CTTN p Y154', 'CTTN p Y162', 'CTTN p Y334', 'IRS1 p Y632', 'IRS1 p Y941', 'IRS1 p Y989', 'ALK p Y1586', 'CTNND1 p Y228', 'CTNND1 p Y302', 'CTNND1 p Y302', 'CTTN p Y154', 'CTTN p Y162', 'CTTN p Y162', 'CTTN p Y334', 'IRS1 p Y632', 'IRS1 p Y989', 'IRS1 p Y989', 'IRS1 p Y989', 'IRS1', 'CTTN', 'IRS1', 'NPM1', 'NPM1'), interaction=c('peptide', 'peptide', 'peptide', 'peptide', 'peptide', 'peptide', 'peptide', 'peptide', 'peptide', 'peptide', 'peptide', 'peptide', 'peptide', 'peptide', 'peptide', 'peptide', 'positive correlation', 'positive correlation', 'positive correlation', 'positive correlation', 'positive correlation', 'positive correlation', 'positive correlation', 'positive correlation', 'positive correlation', 'positive correlation', 'positive correlation', 'positive correlation', 'Physical interactions, experiments', 'controls-state-change-of, controls-transport-of, Pathway, Physical interactions, experiments', 'Predicted', 'experiments, Physical interactions', 'experiments, Physical interactions'), Weight=c(100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 0.8060606, 0.7575758, 0.7454545, 0.9393939, 0.8949096, 0.7329699, 0.7553845, 0.7866191, 0.775, 0.6969697, 0.7818182, 0.8424242, 0.6123365, 2.115272, 0.002461723, 0.3354451, 0.5661711), Alt.Weight=c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.8060606, 0.7575758, 0.7454545, 0.9393939, 0.8949096, 0.7329699, 0.7553845, 0.7866191, 0.775, 0.6969697, 0.7818182, 0.8424242, 0.6123365, 2.115272, 0.002461723, 0.3354451, 0.5661711))
	
	net.cf <- data.frame(id=c('ALK', 'ALK p Y1078', 'ALK p Y1096', 'ALK p Y1586', 'CTNND1', 'CTNND1 p Y193', 'CTNND1 p Y217', 'CTNND1 p Y228', 'CTNND1 p Y241', 'CTNND1 p Y248', 'CTNND1 p Y302', 'CTNND1 p Y904', 'CTTN', 'CTTN p Y154', 'CTTN p Y162', 'CTTN p Y334', 'IRS1', 'IRS1 p Y632', 'IRS1 p Y941', 'IRS1 p Y989', 'NPM1'), Gene.Name=c('ALK', 'ALK', 'ALK', 'ALK', 'CTNND1', 'CTNND1', 'CTNND1', 'CTNND1', 'CTNND1', 'CTNND1', 'CTNND1', 'CTNND1', 'CTTN', 'CTTN', 'CTTN', 'CTTN', 'IRS1', 'IRS1', 'IRS1', 'IRS1', 'NPM1'), nodeType=c('receptor tyrosine kinase', 'receptor tyrosine kinase', 'receptor tyrosine kinase', 'receptor tyrosine kinase', 'undefined', 'undefined', 'undefined', 'undefined', 'undefined', 'undefined', 'undefined', 'undefined', 'SH3 protein', 'SH3 protein', 'SH3 protein', 'SH3 protein', 'undefined', 'undefined', 'undefined', 'undefined', 'RNA processing protein'), Total=c(-209.4885, -11.15512, -7.712652, -9.917654, -783.7109, 0, -16.63529, 0, -95.21879, -35.03751, 0, -74.77153, -172.7889, -67.19445, -14.73001, -55.37009, -401.1118, -6.4122, -6.388159, -5.171309, 1.580289), parent=c('', 'ALK', 'ALK', 'ALK', '', 'CTNND1', 'CTNND1', 'CTNND1', 'CTNND1', 'CTNND1', 'CTNND1', 'CTNND1', '', 'CTTN', 'CTTN', 'CTTN', '', 'IRS1', 'IRS1', 'IRS1', ''), Node.ID=c('gene', 'peptide', 'peptide', 'peptide', 'gene', 'peptide', 'peptide', 'peptide', 'peptide', 'peptide', 'peptide', 'peptide', 'gene', 'peptide', 'peptide', 'peptide', 'gene', 'peptide', 'peptide', 'peptide', 'gene'))
	
	net3.suid <- createNetworkFromDataFrames(net.cf, net.edges, title=paste("Example", 1+length(getNetworkList())), collection = "Keys")
	layoutNetwork(net.w, "genemania-force-directed")
	#setEdgeLineWidthRule(net.w, edge.attribute.name="Weight", attribute.values=net.edges$Weight, line.widths=net.edges$Weight)
	edgeDprops.RCy32()
	setEdgeWidths.RCy32(edgefile=net.edges, factor=1.4, log=TRUE)
	nodeDprops.RCy32(net.cf)
	ratioProps.RCy32(net.cf, "Total")  
	setEdgeWidths.log(net.w, factor=1.2)
#--------------------------------------------------------------------------------	
     # Alternative layout 
	layoutNetwork(layout.name="kamada-kawai")
	# A bit better, not great.
	
# Layout property access
	### I would like the parent or gene nodes to be driving the layouts and the peptides to surround them. For this reason I set the edge weight from gene to peptide to be large. But my attempts to tweak existing layouts have not accomplished this goal.  A workaround could be 
	selectNodes(net.cf[which(net.cf$Node.ID=="gene"), "id"], by="id") 
	# The following use to work in RCytoscape. But NOT in RCy3
	setLayoutProperties("force-directed", list(selected_only=TRUE) )
	# But: selected_only is not a property in layout force-directed	
	setLayoutProperties("kamada-kawai", list(selected_only=TRUE) )
	# selected_only is not a property in layout kamada-kawai	
	layouts <- getLayoutNames()
	getLayoutPropertyNames("genemania-force-directed")
	getLayoutPropertyNames("kamada-kawai")
### Neither of these specify the quantitative edge attribute to drive the layout. This is important because I would sometimes like to use Alt.Weight instead of Weight, for example to drive negative correlation edges far compared to postitive correlation edges. 
	# The following are my unsuccessful tweaks
	defalutgmvalues <- 	data.frame(name=getLayoutPropertyNames("genemania-force-directed"), value=getLayoutPropertyValue(cy, "genemania-force-directed", getLayoutPropertyNames(cy, "genemania-force-directed")))
	# Has to be manual without the ability to layout with selected_only=TRUE
	# Manually: genemania force directed with weight works okay, after manually resizing
	setLayoutProperties("genemania-force-directed", list(numIterations=100, defaultSpringCoefficient=0.1, defaultSpringLength=50, minNodeMass=0.001, maxNodeMass=2000, midpointEdges=250, curveSteepness=7.0e-03, isDeterministic=1, singlePartition=0, ignoreHiddenElements=1))
	layoutNetwork(layout.name= "genemania-force-directed")
	# Not bad, but not exacly heat I was trying to do.
	defalutkkvalues <- 	data.frame(name=getLayoutPropertyNames(cy, "kamada-kawai"), value=getLayoutPropertyValue(cy, "kamada-kawai", getLayoutPropertyNames(cy, "kamada-kawai")))
	setLayoutProperties ("kamada-kawai", list (m_averageIterationsPerNode=50, m_nodeDistanceStrengthConstant=8000, m_nodeDistanceRestLengthConstant=400, m_disconnectedNodeDistanceSpringStrength=0.0001, m_disconnectedNodeDistanceSpringRestLength=2000, m_anticollisionSpringStrength=0, m_layoutPass=50, singlePartition=10, unweighted=1, randomize=10)) 
	layoutNetwork(layout.name="kamada-kawai")
	# Not happy with this. 
# - - - - - - - - - - - - - - - - -
	# Idea: use groups!	
	genes <- net.cf[grep("gene", net.cf$Node.ID), "id"]
	sapply(genes, function(x) createGroup(x, nodes=net.cf[grep(x, net.cf$Gene.Name), "id"], nodes.by.col = "id"))
	collapseGroup(genes)
	# Move nodes around or
	layoutNetwork(layout.name="kamada-kawai")
	expandGroup(genes)
	# This could work for figures, still requires moving PTM nodes around after expanding, but it kind of serves the purpose. 
	# NOTE CRASH: if you move nodes around in the group networks, then try to expandGroup - Cytoscape hangs!
	
# NOT RUN:	
sessionInfo()
R version 3.5.0 (2018-04-23)
Platform: x86_64-apple-darwin15.6.0 (64-bit)
Running under: macOS High Sierra 10.13.5

Matrix products: default
BLAS: /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
LAPACK: /Library/Frameworks/R.framework/Versions/3.5/Resources/lib/libRlapack.dylib

locale:
    [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

attached base packages:
    [1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
    [1] BiocStyle_2.8.2      BiocInstaller_1.30.0 igraph_1.2.1         gplots_3.0.1         RColorBrewer_1.1-2  
[6] dplyr_0.7.5          plyr_1.8.4           RCy3_2.0.3           devtools_1.13.5      knitr_1.20          

loaded via a namespace (and not attached):
    [1] Rcpp_0.12.17        pillar_1.2.3        compiler_3.5.0      bindr_0.1.1         base64enc_0.1-3    
[6] bitops_1.0-6        tools_3.5.0         digest_0.6.15       jsonlite_1.5        evaluate_0.10.1    
[11] memoise_1.1.0       tibble_1.4.2        pkgconfig_2.0.1     rlang_0.2.1         graph_1.58.0       
[16] curl_3.2            yaml_2.1.19         parallel_3.5.0      bindrcpp_0.2.2      stringr_1.3.1      
[21] withr_2.1.2         httr_1.3.1          gtools_3.5.0        caTools_1.17.1      rprojroot_1.3-2    
[26] stats4_3.5.0        tidyselect_0.2.4    glue_1.2.0          R6_2.2.2            XML_3.98-1.11      
[31] rmarkdown_1.10      RJSONIO_1.3-0       gdata_2.18.0        purrr_0.2.5         magrittr_1.5       
[36] backports_1.1.2     htmltools_0.3.6     BiocGenerics_0.26.0 assertthat_0.2.0    KernSmooth_2.23-15 
[41] stringi_1.2.3 