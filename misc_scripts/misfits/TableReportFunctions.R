#---------------------------------------------------------------------
# Script Name:   TableReportFunctions.R
# Description:   Functions to deal with TASSEL TableReport
# Author:        Brandon Monier & Ed buckler
# Created:       2018-11-26 at 11:14:36
# Last Modified: 2018-12-21 at 15:16:56
#--------------------------------------------------------------------

#--------------------------------------------------------------------
# Detailed Purpose:
#    The main purpose of this Rscript produce wrapper classes for
#    TASSEL classes
#--------------------------------------------------------------------




convertTableReportToDataFrame <- function(tableReport) {
  tableReportVectors <- J("net/maizegenetics/plugindef/GenerateRCode")$tableReportToVectors(tableReport)
  colNum <- length(tableReportVectors$columnNames)
  aDF <- data.frame(tableReportVectors$dataVector$get(0L))
  for(i in 2:(colNum)) {
    aDF[[i]] <- tableReportVectors$dataVector$get(i-1L)
  }
  colnames(aDF) <- tableReportVectors$columnNames
  aDF
}


#These should all be vectorized

tasselObjToRWrapper <- function(tasselObj) {
  #If already a TASSEL object just return it
  
}

rObjToTasselObj <- function(rObj) {
  #If already a TASSEL object just return it
  
}



