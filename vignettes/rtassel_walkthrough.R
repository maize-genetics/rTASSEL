## ----setup, include=FALSE--------------------------------------------------
knitr::opts_chunk$set(
    fig.path='figure/graphics-', 
    cache.path='cache/graphics-', 
    fig.align='center',
    external=TRUE,
    echo=TRUE,
    warning=FALSE
    # fig.pos="H"
)

## ---- eval=TRUE, echo=FALSE, message=FALSE, warning=FALSE------------------
rTASSEL::startLogger()

## ---- eval=FALSE, echo=TRUE------------------------------------------------
#  if (!require("devtools")) install.packages("devtools")
#  devtools::install_bitbucket(repo = "bucklerlab/rtassel", ref = "master")

## ---- eval=FALSE, echo=TRUE------------------------------------------------
#  library(rTASSEL)

## ---- eval=FALSE, echo=TRUE------------------------------------------------
#  rTASSEL::startLogger(fullPath = NULL, fileName = NULL)

## ---- eval=FALSE, echo=TRUE------------------------------------------------
#  options(java.parameters = c("-Xmx<memory>", "-Xms<memory>"))

## ---- eval=TRUE, echo=TRUE-------------------------------------------------
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

## ---- eval=TRUE, echo=TRUE-------------------------------------------------
# Load in hapmap file
tasGenoHMP <- rTASSEL::readGenotypeTableFromPath(
    path = genoPathHMP
)

# Load in VCF file
tasGenoVCF <- rTASSEL::readGenotypeTableFromPath(
    path = genoPathVCF
)

## ---- eval=TRUE, echo=TRUE-------------------------------------------------
tasGenoHMP

## ---- eval=TRUE, echo=TRUE-------------------------------------------------
class(tasGenoHMP)
slotNames(tasGenoHMP)

## ---- eval=TRUE, echo=TRUE-------------------------------------------------
tasGenoHMP@jGenotypeTable

## ---- eval=TRUE, echo=TRUE-------------------------------------------------
# Read from phenotype path
phenoPath  <- system.file("extdata", "mdp_traits.txt", package = "rTASSEL")
phenoPath

# Load into rTASSEL `TasselGenotypePhenotype` object
tasPheno <- rTASSEL::readPhenotypeFromPath(
    path = phenoPath
)

# Inspect object
tasPheno

## ---- eval=TRUE, echo=TRUE-------------------------------------------------
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

## ---- eval=TRUE, echo=TRUE-------------------------------------------------
tasGenoPheno <- rTASSEL::readGenotypePhenotype(
    genoPathOrObj = tasGenoHMP,
    phenoPathDFOrObj = tasPheno
)
tasGenoPheno

## ---- eval=TRUE, echo=TRUE-------------------------------------------------
tasGenoPhenoDF <- rTASSEL::readGenotypePhenotype(
    genoPathOrObj = genoPathHMP,
    phenoPathDFOrObj = phenoDF,
    taxaID = "Taxon",
    attributeTypes = NULL
)
tasGenoPhenoDF

## ---- eval=TRUE, echo=TRUE-------------------------------------------------
tasSumExp <- rTASSEL::getSumExpFromGenotypeTable(
    tasObj = tasGenoPheno
)
tasSumExp

## ---- eval=TRUE, echo=TRUE-------------------------------------------------
SummarizedExperiment::colData(tasSumExp)

## ---- eval=TRUE, echo=TRUE-------------------------------------------------
SummarizedExperiment::rowData(tasSumExp)

## ---- eval=TRUE, echo=TRUE-------------------------------------------------
SummarizedExperiment::rowRanges(tasSumExp)

## ---- eval=TRUE, echo=TRUE-------------------------------------------------
tasExportPhenoDF <- rTASSEL::getPhenotypeDF(
    tasObj = tasGenoPheno
)
tasExportPhenoDF

## ---- eval=TRUE, echo=TRUE-------------------------------------------------
tasGenoPhenoFilt <- rTASSEL::filterGenotypeTableSites(
    tasObj = tasGenoPheno,
    siteMinCount = 150,
    siteMinAlleleFreq = 0.05,
    siteMaxAlleleFreq = 1.0
)
tasGenoPhenoFilt

## ---- eval=TRUE, echo=TRUE-------------------------------------------------
tasGenoPheno

## ---- eval=TRUE, echo=TRUE-------------------------------------------------
tasKin <- rTASSEL::kinshipMatrix(tasObj = tasGenoPheno)

## ---- eval=TRUE, echo=TRUE-------------------------------------------------
# Get full R matrix
tasKinRMat <- rTASSEL::kinshipToRMatrix(tasKin)

# Inspect the first 5 rows and columns
tasKinRMat[1:5, 1:5]

## ---- eval=TRUE, echo=TRUE-------------------------------------------------
tasDist <- rTASSEL::distanceMatrix(tasObj = tasGenoPheno)

## ---- eval=TRUE, echo=TRUE-------------------------------------------------
# Get full R matrix
tasDistRMat <- rTASSEL::distanceToRMatrix(tasDist)

# Inspect the first 5 rows and columns
tasDistRMat[1:5, 1:5]

## ---- echo=TRUE, eval=TRUE-------------------------------------------------
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
tasBLUE

## ---- echo=TRUE, eval=TRUE-------------------------------------------------
# Calculate GLM
tasGLM <- rTASSEL::assocModelFitter(
    tasObj = tasGenoPheno,             # <- our prior TASSEL object
    formula = list(EarHT, dpoll) ~ .,  # <- only EarHT and dpoll are ran
    fitMarkers = TRUE,                 # <- set this to TRUE for GLM
    kinship = NULL,
    fastAssociation = FALSE
)

# Return GLM output
tasGLM

## ---- echo=TRUE, eval=TRUE-------------------------------------------------
# Calculate MLM
tasMLM <- rTASSEL::assocModelFitter(
    tasObj = tasGenoPheno,             # <- our prior TASSEL object
    formula = EarHT ~ .,               # <- run only EarHT
    fitMarkers = TRUE,                 # <- set this to TRUE for GLM
    kinship = tasKin,                  # <- our prior kinship object
    fastAssociation = FALSE
)

# Return GLM output
tasMLM

## ---- echo=TRUE, eval=TRUE-------------------------------------------------
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
tasFAST

