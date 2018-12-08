#--------------------------------------------------------------------
# Script Name:   AllGenerics.R
# Description:   Various tests with rJava
# Author:        Brandon Monier & Ed buckler
# Created:       2018-11-26 at 11:14:36
# Last Modified: 2018-12-03 at 17:58:46
#--------------------------------------------------------------------

#--------------------------------------------------------------------
# Detailed Purpose:
#    The main purpose of this Rscript produce wrapper classes for 
#    TASSEL S4 generics
#--------------------------------------------------------------------

setGeneric(
    name = "positions",
    def = function(object) {
        standardGeneric("positions")
    }
)

setGeneric(
  name = "taxa",
  def = function(object) {
    standardGeneric("taxa")
  }
)