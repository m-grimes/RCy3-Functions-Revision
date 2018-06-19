# RCy3-Functions-Revision
Conversion of RCy3 functions to the updated direct calls to Cytoscape using cyREST

Author: Mark Grimes
References: http://hs.umt.edu/dbs/labs/grimes/publications.php
The new fuctions are based on Alex Pico's cyREST work and suggestions.

### Note: For  For createNetworkFromDataFrames {RCy3}, in the node file, the first column must be named "id". 
### In edge files, the first two columns must be named "source" and "target"; and "edgeType" is no longer used, replaced with "interaction". For converting files we may wish to use the edgeTypeToInteraction() function


Functions entered so far: edgeType.to.interaction; mergeEdges.dir; mergeEdges.RCy32; edgeDprops.RCy32; nodeDprops.RCy32; setEdgeWidths.RCy32; intensityProps.RCy32; getCyEdgeNames.RCy32; extract.peptides.RCy32; graphNetworkPath.RCy32.

Also created a function to collapse all or selective nodes: collapse.CCCN.nodes.



