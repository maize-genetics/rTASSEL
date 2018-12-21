# Functions for plotting data from a snpMatrix or data frame

## Function to calulate Fst for many ref ranges and return weighted Fst in dataframe
## Requires a list of variables of the snpMatrix class each of which corresponds to a 
## specific reference range, and a dataframe with rownames = taxa and a single column
## designating a group for each taxon.
refRangeFst <- function(snpMatrix_list, taxa_groups) {
  df_out <- matrix(ncol=1, nrow=length(snpMatrix_list))
  for (i in 1:length(snpMatrix_list)) {
    snpmx <- get(snpMatrix_list[i])
    snpFst <- Fst(snpmx, group = taxa_groups[, 1])
    meanFst <- weighted.mean(snpFst$Fst, snpFst$weight)
    df_out[i,] <- meanFst
  }
  return(df_out)
}


## Calculate distance matrix and create upgma tree for taxa in a snpMatrix
## Recommended to subset the snp matrix to a smaller number of taxa for easy viewing
## Function expects a variable of class snpMatrix
library(phangorn)
createUPGMATree <- function(snpmx_subset) {
  ibs <- ibsCount(snpmx_subset)
  distance <- ibsDist(ibs)
  tree <- upgma(distance)
}

## Plot per-line missingness from HapsInRefRange data frame
# char_count <- function(x, c) {nchar(x) - nchar(gsub(c, '', x)) }
# charcount = nchar(HapsInRefRange$Sequence)
# Ncount = char_count(HapsInRefRange$Sequence, 'N')
# tapply(X=Ncount/charcount, INDEX = HapsInRefRange$HapID, FUN = mean)
# missingness = tapply(X=Ncount/charcount, INDEX = HapsInRefRange$HapID, FUN = mean)




## Examples
# Simulate 100 snp matrices, each with 20 taxa and 1000 SNPs
for (i in 1:100) {
  randsnps <- random.snps(20, 1000)
  assign(paste("snpset", i, sep = "_"), randsnps)
}

# Plot Fst across simulated ref ranges
Fst_test = c(paste0("snpset_", 1:100))
group_test = as.data.frame(c(rep('group1', 10), rep('group2', 10)))
rownames(group_test) = rownames(snpset_1)
colnames(group_test) = c("Group")

testOut2 = refRangeFst(snpMatrix_list = Fst_test, taxa_groups = group_test)
plot(1:100, testOut2, xlab = "Reference Range", ylab = "Fst", type = 'b')

# UPGMA tree for single ref range
tree_test = createUPGMATree(snpset_1)
plot(tree_test)

