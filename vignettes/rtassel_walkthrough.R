## ----setup, include=FALSE-----------------------------------------------------
library(rTASSEL)

knitr::opts_chunk$set(
    fig.path='figure/graphics-',
    cache.path='cache/graphics-',
    fig.align='center',
    external=TRUE,
    echo=TRUE,
    warning=FALSE
    # fig.pos="H"
)

## ---- eval=FALSE, echo=TRUE---------------------------------------------------
#  options(java.parameters = c("-Xmx<memory>", "-Xms<memory>"))

## ---- eval=FALSE, echo=TRUE---------------------------------------------------
#  rTASSEL::startLogger(fullPath = NULL, fileName = NULL)

## ---- eval=TRUE, echo=TRUE----------------------------------------------------
# Load hapmap data
genoPathHMP <- system.file(
    "extdata",
    "mdp_genotype.hmp.txt",
    package = "rTASSEL"
)
genoPathHMP

# Load VCF data
genoPathVCF <- system.file(
    "extdata",
    "maize_chr9_10thin40000.recode.vcf",
    package = "rTASSEL"
)
genoPathVCF

## ---- eval=TRUE, echo=TRUE----------------------------------------------------
# Load in hapmap file
tasGenoHMP <- rTASSEL::readGenotypeTableFromPath(
    path = genoPathHMP
)

# Load in VCF file
tasGenoVCF <- rTASSEL::readGenotypeTableFromPath(
    path = genoPathVCF
)

## ---- eval=TRUE, echo=TRUE----------------------------------------------------
tasGenoHMP

## ---- eval=TRUE, echo=TRUE----------------------------------------------------
class(tasGenoHMP)
slotNames(tasGenoHMP)

## ---- eval=TRUE, echo=TRUE----------------------------------------------------
tasGenoHMP@jGenotypeTable

## ---- eval=TRUE, echo=TRUE----------------------------------------------------
# Read from phenotype path
phenoPath  <- system.file("extdata", "mdp_traits.txt", package = "rTASSEL")
phenoPath

# Load into rTASSEL `TasselGenotypePhenotype` object
tasPheno <- rTASSEL::readPhenotypeFromPath(
    path = phenoPath
)

# Inspect object
tasPheno

## ---- eval=TRUE, echo=TRUE----------------------------------------------------
# Create phenotype data frame
phenoDF <- read.table(phenoPath, header = TRUE)
colnames(phenoDF)[1] <- "Taxon"

# Inspect first few rows
head(phenoDF)

# Load into rTASSEL `TasselGenotypePhenotype` object
tasPhenoDF <- rTASSEL::readPhenotypeFromDataFrame(
    phenotypeDF = phenoDF,
    taxaID = "Taxon",
    attributeTypes = NULL
)

# Inspect new object
tasPhenoDF

## ---- eval=TRUE, echo=TRUE----------------------------------------------------
tasGenoPheno <- rTASSEL::readGenotypePhenotype(
    genoPathOrObj = tasGenoHMP,
    phenoPathDFOrObj = tasPheno
)
tasGenoPheno

## ---- eval=TRUE, echo=TRUE----------------------------------------------------
tasGenoPhenoDF <- rTASSEL::readGenotypePhenotype(
    genoPathOrObj = genoPathHMP,
    phenoPathDFOrObj = phenoDF,
    taxaID = "Taxon",
    attributeTypes = NULL
)
tasGenoPhenoDF

## ---- eval=TRUE, echo=TRUE----------------------------------------------------
## Get toy kinship data from package ----
kinshipPath <- system.file(
  "extdata", 
  "mdp_kinship.txt", 
  package = "rTASSEL"
)

## Read ----
rTASSEL::readTasselDistanceMatrix(kinshipPath)

## ---- eval=TRUE, echo=TRUE----------------------------------------------------
tasSumExp <- rTASSEL::getSumExpFromGenotypeTable(
    tasObj = tasGenoPheno
)
tasSumExp

## ---- eval=TRUE, echo=TRUE----------------------------------------------------
SummarizedExperiment::colData(tasSumExp)

## ---- eval=TRUE, echo=TRUE----------------------------------------------------
SummarizedExperiment::rowData(tasSumExp)

## ---- eval=TRUE, echo=TRUE----------------------------------------------------
SummarizedExperiment::rowRanges(tasSumExp)

## ---- eval=TRUE, echo=TRUE----------------------------------------------------
tasExportPhenoDF <- rTASSEL::getPhenotypeDF(
    tasObj = tasGenoPheno
)
tasExportPhenoDF

## ---- eval=TRUE, echo=TRUE----------------------------------------------------
tasGenoPhenoFilt <- rTASSEL::filterGenotypeTableSites(
    tasObj = tasGenoPheno,
    siteMinCount = 150,
    siteMinAlleleFreq = 0.05,
    siteMaxAlleleFreq = 1.0,
    siteRangeFilterType = "none"
)
tasGenoPhenoFilt

## ---- eval=TRUE, echo=TRUE----------------------------------------------------
tasGenoPheno

## ---- eval=TRUE, echo=TRUE----------------------------------------------------
tasKin <- rTASSEL::kinshipMatrix(tasObj = tasGenoPheno)

## ---- eval=TRUE, echo=TRUE----------------------------------------------------
tasDist <- rTASSEL::distanceMatrix(tasObj = tasGenoPheno)

## ---- eval=TRUE, echo=TRUE----------------------------------------------------
tasKin

## ---- eval=TRUE, echo=TRUE----------------------------------------------------
library(magrittr)

tasKin %>% colnames() %>% head()
tasKin %>% rownames() %>% head()

tasKin %>% dim()

tasKin %>% nrow()
tasKin %>% ncol()

## ---- eval=TRUE, echo=TRUE----------------------------------------------------
library(magrittr)

tasKinR <- tasKin %>% as.matrix()

## Inspect first 5 rows and columns ----
tasKinR[1:5, 1:5]

## ---- eval=TRUE, echo=TRUE----------------------------------------------------
library(magrittr)

## Create a dummy pairwise matrix object ----
set.seed(123)
m <- 10
s <- matrix(rnorm(100), m)
s[lower.tri(s)] <- t(s)[lower.tri(s)]
diag(s) <- 2

## Add sample IDs ----
colnames(s) <- rownames(s) <- paste0("s_", seq_len(m))

testTasselDist <- s %>% asTasselDistanceMatrix()
testTasselDist

## ---- eval=TRUE, echo=TRUE----------------------------------------------------
tasGenoHMP

pcaRes <- pca(tasGenoHMP)

## ---- eval=TRUE, echo=TRUE----------------------------------------------------
tasDist

mdsRes <- mds(tasDist)

## ---- eval=TRUE, echo=TRUE----------------------------------------------------
pcaRes

## ---- echo=TRUE, eval=TRUE----------------------------------------------------
# Read in phenotype data
phenoPathCov <- system.file("extdata", "mdp_phenotype.txt", package = "rTASSEL")
tasPhenoCov <- rTASSEL::readPhenotypeFromPath(phenoPathCov)

# Calculate BLUEs
tasBLUE <- rTASSEL::assocModelFitter(
    tasObj = tasPhenoCov,
    formula = . ~ .,                  # <- All data is used!
    fitMarkers = FALSE,
    kinship = NULL,
    fastAssociation = FALSE
)

# Return BLUE output
str(tasBLUE)

## ---- echo=TRUE, eval=TRUE----------------------------------------------------
# Calculate GLM
tasGLM <- rTASSEL::assocModelFitter(
    tasObj = tasGenoPheno,             # <- our prior TASSEL object
    formula = list(EarHT, dpoll) ~ .,  # <- only EarHT and dpoll are ran
    fitMarkers = TRUE,                 # <- set this to TRUE for GLM
    kinship = NULL,
    fastAssociation = FALSE
)

# Return GLM output
str(tasGLM)

## ---- echo=TRUE, eval=TRUE----------------------------------------------------
# Calculate MLM
tasMLM <- rTASSEL::assocModelFitter(
    tasObj = tasGenoPheno,             # <- our prior TASSEL object
    formula = EarHT ~ .,               # <- run only EarHT
    fitMarkers = TRUE,                 # <- set this to TRUE for GLM
    kinship = tasKin,                  # <- our prior kinship object
    fastAssociation = FALSE
)

# Return GLM output
str(tasMLM)

## ---- echo=TRUE, eval=TRUE----------------------------------------------------
# Read data - need only non missing data!
phenoPathFast <-system.file(
    "extdata",
    "mdp_traits_nomissing.txt",
    package = "rTASSEL"
)

# Creat rTASSEL object - use prior TASSEL genotype object
tasGenoPhenoFast <- rTASSEL::readGenotypePhenotype(
    genoPathOrObj = tasGenoHMP,
    phenoPathDFOrObj = phenoPathFast
)


# Calculate MLM
tasFAST <- rTASSEL::assocModelFitter(
    tasObj = tasGenoPhenoFast,         # <- our prior TASSEL object
    formula = . ~ .,                   # <- run all of the phenotype data
    fitMarkers = TRUE,                 # <- set this to TRUE for GLM
    kinship = NULL,
    fastAssociation = TRUE             # <- set this to TRUE for fast assoc.
)

# Return GLM output
str(tasFAST)

## ---- echo=TRUE, eval=TRUE----------------------------------------------------
phyloTree <- createTree(
    tasObj = tasGenoHMP,
    clustMethod = "Neighbor_Joining"
)

## ---- echo=TRUE, eval=TRUE----------------------------------------------------
phyloTree

## ---- echo=TRUE, eval=FALSE---------------------------------------------------
#  # Generate Manhattan plot for ear height trait
#  manhattanEH <- manhattanPlot(
#      assocStats = tasMLM$MLM_Stats,
#      trait      = "EarHT",
#      threshold  = 5
#  )
#  
#  manhattanEH

## ---- echo=TRUE, eval=FALSE---------------------------------------------------
#  # Filter genotype table by position
#  tasGenoPhenoFilt <- filterGenotypeTableSites(
#      tasObj              = tasGenoPheno,
#      siteRangeFilterType = "position",
#      startPos            = 228e6,
#      endPos              = 300e6,
#      startChr            = 2,
#      endChr              = 2
#  )
#  
#  # Generate and visualize LD
#  myLD <- ldPlot(
#      tasObj  = tasGenoPhenoFilt,
#      ldType  = "All",
#      plotVal = "r2",
#      verbose = FALSE
#  )
#  
#  myLD

## ---- echo=TRUE, eval=FALSE---------------------------------------------------
#  library(magrittr)
#  
#  tasGenoHMP %>% ldJavaApp(windowSize = 100)

## ---- echo=TRUE, eval=FALSE---------------------------------------------------
#  library(magrittr)
#  
#  tasGenoHMP %>%
#      filterGenotypeTableTaxa(
#        taxa = c("33-16", "38-11", "4226", "4722", "A188", "A214N")
#      ) %>%
#      treeJavaApp()

## ---- echo=TRUE, eval=FALSE---------------------------------------------------
#  tasCV <- genomicPrediction(
#      tasPhenoObj = tasGenoPheno,
#      kinship     = tasKin,
#      doCV        = TRUE,
#      kFolds      = 5,
#      nIter       = 1
#  )
#  head(tasCV)

