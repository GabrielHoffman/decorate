resDiffCorr = function(clusterValues, param, npermute=100){

#   clustValTab = table(clusterValues)
#   clustValTab = clustValTab[clustValTab >3]
#   clustIDLst = names(clustValTab)

#   sledRes = bplapply( clustIDLst, function( clstLabel ){

#     cat("\rCluster label", clstLabel)
#     peakIDs = names(clusterValues)[clusterValues == clstLabel]

#     if( length(peakIDs) > 3 ){

#       Y1 = t(vobj$E[peakIDs,metadata$Disease == 0])
#       Y2 = t(vobj$E[peakIDs,metadata$Disease == 1])

#       res = sLED(X=Y1, Y=Y2, npermute=npermute, verbose=FALSE)
#     }else{
#       res = NULL
#     }

#     res
#   }, BPPARAM=param)
#   names(sledRes) = clustIDLst

#   sledRes
# }