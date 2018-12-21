GWASpolyGenoFromSummarizedExperiment <- function(SummarizedExperimentObject){
  library(SummarizedExperiment)
  geno <- data.frame(markerName = paste("dummy", 1:length(ranges(SummarizedExperimentObject@rowRanges)), sep = "-"), # dummy name as current summarizeExperimentFromGenotypeTable doesn't keep
                     chr = seqnames(SummarizedExperimentObject@rowRanges),
                     pos = start(ranges(SummarizedExperimentObject@rowRanges)),
                     as.data.frame(SummarizedExperimentObject@assays$data@listData) # same as assay(SummarizedExperimentObject)
  )
  colnames(geno)[4:ncol(geno)] <- as.character(SummarizedExperimentObject$Sample)
  geno
}

convert.snp <- function(x) {#convert to {-1,0,1,NA}
  alleles <- na.omit(unique(x))
  y <- rep(NA,length(x))
  y[which(x==alleles[1])] <- -1
  y[which(x==alleles[2])] <- 1
  return(y)
}

# Coopting read.gwaspoly to create GWASpoly class object from tassel simplified experiment
#right now it makes dummy marker names as current summarizeExperimentFromGenotypeTable doesn't keep them
se_createGWASpolyObject <- function(ploidy, phenoDF, SummarizedExperimentObject, format, n.traits)
{
  if (format == "ACTG") {
    format <- "ACGT"
  }
  if (!is.element(format, c("AB", "numeric", "ACGT"))) {
    stop("Invalid genotype format.")
  }
  bases <- c("A", "C", "G", "T")
  get.ref <- function(x, format) {
    if (format == "numeric") {
      ref.alt <- c(0, 1)
    }
    if (format == "AB") {
      ref.alt <- c("A", "B")
    }
    if (format == "ACGT") {
      y <- paste(na.omit(x), collapse = "")
      ans <- apply(array(bases), 1, function(z, y) {
        length(grep(z, y, fixed = T))
      }, y)
      if (sum(ans) > 2) {
        stop("Error in genotype matrix: More than 2 alleles")
      }
      if (sum(ans) == 2) {
        ref.alt <- bases[which(ans == 1)]
      }
      if (sum(ans) == 1) {
        ref.alt <- c(bases[which(ans == 1)], NA)
      }
    }
    return(ref.alt)
  }

  ## unwrap tassel summarized experiment to expected GWASpoly genotype format
  geno <- data.frame(markerName = paste("dummy", 1:length(ranges(SummarizedExperimentObject@rowRanges)), sep = "-"), # dummy name as current summarizeExperimentFromGenotypeTable doesn't keep
                     chr = seqnames(SummarizedExperimentObject@rowRanges),
                     pos = start(ranges(SummarizedExperimentObject@rowRanges)),
                     as.data.frame(SummarizedExperimentObject@assays$data@listData), # same as assay(SummarizedExperimentObject)
                     stringsAsFactors = F
  )
  colnames(geno)[4:ncol(geno)] <- as.character(SummarizedExperimentObject$Sample)

  map <- data.frame(Marker = geno[, 1], Chrom = factor(geno[,
                                                            2], ordered = T), Position = geno[, 3], stringsAsFactors = F)
  markers <- as.matrix(geno[, -(1:3)])
  rownames(markers) <- geno[, 1]
  tmp <- apply(markers, 1, get.ref, format)
  map$Ref <- tmp[1, ]
  map$Alt <- tmp[2, ]
  if (is.element(format, c("AB", "ACGT"))) {
    M <- apply(cbind(map$Ref, markers), 1, function(x) {
      y <- gregexpr(pattern = x[1], text = x[-1], fixed = T)
      ans <- as.integer(lapply(y, function(z) {
        ifelse(z[1] < 0, ploidy, ploidy - length(z))
      }))
      return(ans)
    })
  }
  else {
    M <- t(markers)
  }
  gid.geno <- colnames(geno)[-(1:3)]
  rownames(M) <- gid.geno
  bad <- length(which(!is.element(na.omit(M), 0:ploidy)))
  if (bad > 0) {
    stop("Invalid marker calls.")
  }
  MAF <- apply(M, 2, function(x) {
    AF <- mean(x, na.rm = T)/ploidy
    MAF <- ifelse(AF > 0.5, 1 - AF, AF)
  })
  polymorphic <- which(MAF > 0)
  M <- M[, polymorphic]
  map <- map[polymorphic, ]
  map <- map[order(map$Chrom, map$Position), ]
  M <- M[, map$Marker]
  m <- nrow(map)
  cat(paste("Number of polymorphic markers:", m, "\n"))
  impute.mode <- function(x) {
    ix <- which(is.na(x))
    if (length(ix) > 0) {
      x[ix] <- as.integer(names(which.max(table(x))))
    }
    return(x)
  }
  missing <- which(is.na(M))
  if (length(missing) > 0) {
    cat("Missing marker data imputed with population mode \n")
    M <- apply(M, 2, impute.mode)
  }
  pheno <- phenoDF
  gid.pheno <- unique(pheno[, 1])
  gid <- intersect(gid.pheno, gid.geno)
  pheno <- pheno[is.element(pheno[, 1], gid), ]
  M <- M[gid, ]
  N <- length(gid)
  cat(paste("N =", N, "individuals with phenotypic and genotypic information \n"))
  n.fixed <- ncol(pheno) - n.traits - 1
  if (n.fixed > 0) {
    fixed <- data.frame(pheno[, (n.traits + 2):ncol(pheno)],
                        stringsAsFactors = F)
    fixed.names <- colnames(pheno)[(n.traits + 2):ncol(pheno)]
    colnames(fixed) <- fixed.names
    pheno <- data.frame(pheno[, 1:(1 + n.traits)], stringsAsFactors = F)
    cat(paste("Detected following fixed effects:\n", paste(fixed.names,
                                                           collapse = "\n"), "\n", sep = ""))
  }
  else {
    fixed <- data.frame(NULL)
  }
  traits <- colnames(pheno)[-1]
  cat(paste("Detected following traits:\n", paste(traits, collapse = "\n"),
            "\n", sep = ""))
  return(new("GWASpoly", map = map, pheno = pheno, fixed = fixed,
             geno = M, ploidy = ploidy))
}


#Function to unwrap GWASpoly result object to tidy DF
gwasPolyToDF <- function(gwasPolyRes=NA){

  traitGWASresults <- data.frame()
  for(Trait in names(gwasPolyRes@scores)){
    message(paste("getting results for trait:", Trait), appendLF = F)

    SigTreshold <- gwasPolyRes@threshold[rownames(gwasPolyRes@threshold)==Trait]

    for(Model in colnames(gwasPolyRes@scores[[Trait]])){
      message(paste(", for model", Model))
      oneTraitGWAS <- data.frame(Marker = rownames(gwasPolyRes@scores[[Trait]]),
                                 MarkerpVal = 10^-gwasPolyRes@scores[[Trait]][,Model],
                                 MarkerLogpVal = gwasPolyRes@scores[[Trait]][,Model],
                                 MarkerEffect = gwasPolyRes@effects[[Trait]][,Model],
                                 Trait, Model, SigTreshold)

      traitGWASresults <- rbind(traitGWASresults, oneTraitGWAS)
    }

  }
  traitGWASresults <- merge(traitGWASresults, gwasPolyRes@map, by = "Marker")
  traitGWASresults
}


tasselGWASPluginResToSingleDF <- function(gwasPluginResultDF, pluginName="fixedEffectLMPlugin"){
  gwasRes <- data.frame()
  if(pluginName == "fixedEffectLMPlugin"){
    gwasRes <- gwasPluginResultDF$"GLM_Stats_FromR"
    gwasRes <- merge(gwasRes, gwasPluginResultDF$"GLM_Genotypes_FromR", by = c("Trait", "Marker", "Pos", "Chr"))
  }
  NAs <- sum(is.na(gwasRes$p))
  warning(paste("Removed", NAs, "Marker-Trait entries due to having 0 df and NA p-value"))
  gwasRes <- gwasRes[-which(is.na(gwasRes$p)),]
  gwasRes
}

#function for generating ggplot of the qq kind
#can specify columns that include relevant required data
#includes option to plot fancy points with conditional format
qqPlot_multi_trait <- function(traitGWASresults, pValIDcol = "MarkerpVal",
                               sigTresholdIDcol = "SigTreshold", traitIDcol = "Trait",
                               markerEffectIDcol = "MarkerEffect", markerIDcol = "Marker", fancy = F){

  traitGWASres = data.frame(Trait =  traitGWASresults[,traitIDcol],
                            MarkerpVal = traitGWASresults[,pValIDcol],
                            MarkerEffect = traitGWASresults[,markerEffectIDcol],
                            Marker = traitGWASresults[,markerIDcol]
                            )

  if(class(sigTresholdIDcol) == "character"){
    traitGWASres$SigTreshold = traitGWASresults[,sigTresholdIDcol]
  }
  if(class(sigTresholdIDcol) == "numeric"){
    traitGWASres$SigTreshold = sigTresholdIDcol
  }

  if(sigTresholdIDcol == ""){
    message("Setting default significance treshold to Bonferroni correction of 0.05")
    traitGWASres$SigTreshold = 0.05/length(unique(traitGWASres$Marker))
  }

  #transform sigTreshold if on -log10 scale
  if(max(traitGWASres$SigTreshold)>1){
    traitGWASres$SigTreshold <- 10^-traitGWASres$SigTreshold
  }

  traitGWASres$Significant <- traitGWASres$MarkerpVal <= traitGWASres$SigTreshold
  traitGWASres$MarkerEffectSign <- as.factor(sign(traitGWASres$MarkerEffect))
  traitGWASres <- arrange(traitGWASres, Trait, MarkerpVal)
  traitGWASres$NormalpVal <- NaN
  for(Trait in unique(traitGWASres$Trait)){
    traitGWASres$NormalpVal[traitGWASres$Trait == Trait] <- ppoints(sum(traitGWASres$Trait == Trait))
  }

  g1 <- ggplot(traitGWASres) +
    geom_point(aes(x = -log10(NormalpVal),
                   y = -log10(MarkerpVal) )) +
    geom_line(aes(x=-log10(NormalpVal),
                  y=-log10(NormalpVal)), lty=2) +
    facet_wrap(~Trait, scales="free") +
    xlab(expression(paste("Expected -log"[10], "(p)", sep = ""))) +
    ylab(expression(paste("Observed -log"[10], "(p)", sep = ""))) +
    theme(legend.position = "bottom")

  if(fancy){
    g1 + geom_point(aes(x = -log10(NormalpVal),
                        y = -log10(MarkerpVal), size = abs(MarkerEffect), col = MarkerEffectSign, pch = Significant))
  } else g1
}

manhattan_trait_plot <-function(traitGWASresults, traitIDcol = "Trait",
                                positionIDcol = "Position", chromIDcol = "Chrom",
                                pValIDcol = "MarkerpVal", sigTresholdIDcol = "SigTreshold",
                                colorPal = "Dark2"){
  library(ggplot2)

  traitGWASres <- data.frame(Chrom = traitGWASresults[,chromIDcol],
                             Position = traitGWASresults[,positionIDcol],
                             pVal = traitGWASresults[,pValIDcol],
                             sigTreshold = traitGWASresults[,sigTresholdIDcol],
                             trait = traitGWASresults[,traitIDcol])

  traitGWASres$markerLogPVal <- abs(log(traitGWASres$pVal))

  ggplot(traitGWASres) +
    geom_point(aes(Position, markerLogPVal, col = trait), alpha=0.3) +
    geom_point(aes(Position, markerLogPVal, col = trait),
               data = traitGWASres[traitGWASres$markerLogPVal>traitGWASres$sigTreshold,]) +
    geom_hline(aes(yintercept = sigTreshold, col = trait), lty = "dashed") +
    facet_grid(~Chrom, scales = "free_x") + theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom") +
    scale_color_brewer(palette = colorPal) + xlab("Position by Chr") + ylab("-Log10Pval")
}

#function to read a gff file and extract gene annotations
#expects first item, as separated by semicolons, on last (9th) column to define the gene name with an id identifier
#e.g for content on 9th column with annotation names: "someGeneID=someGeneName-ver1.alt1;otherID=somethingElse;someMore=more things"
gffToAnnotationGR <- function(gffFile){
  library(stringr)
  library(GenomicRanges)

  gffDF <- read.table(file = gffFile, quote = "", sep = "\t")

  colnames(gffDF) <- c("Chr", "Source", "AnnotType", "Start", "End", "Other", "Strand", "Other2", "IDs")

  gffDF$Chr <- as.factor(gffDF$Chr)

  gffGenes <- gffDF[gffDF$AnnotType %in% c("gene", "Gene"),]

  gffGenes$Gene <- str_split(str_split(gffGenes$IDs, ";", simplify = T)[,1],":", simplify = T)[,2]

  gffGenesGR <- makeGRangesFromDataFrame(gffGenes, keep.extra.columns = T)

  gffGenesGR
}

#function to annotate a gwas results table with genomicRanges annotations.
#annotated by their closest annotation
annotate_gwasRes_byNearest <- function(gwasResDF, annotationGR, positionIDcol_gwas = "Position",  chromIDcol_gwas = "Chrom", outFmt = "data.frame"){
  library(GenomicRanges)

  if(outFmt %in% c("data.frame", "GenomicRanges")  == F){
    stop(paste("outFmt parameter must be 'data.frame' or 'GenomicRanges'"))
  }

  gwasResGR <- makeGRangesFromDataFrame(gwasResDF, seqnames.field = chromIDcol_gwas, start.field = positionIDcol_gwas, end.field = positionIDcol_gwas, keep.extra.columns = T)

  annotNearestToGWASranges <- nearest(gwasResGR, annotationGR, ignore.strand = T)

  distToNearest <- distanceToNearest(gwasResGR, annotationGR, ignore.strand = T)

  annotNearestToGWASranges <- annotationGR[annotNearestToGWASranges]

  annotNearestToGWASranges <- as.data.frame(annotNearestToGWASranges, row.names = NULL)

  annotNearestToGWASranges$distanceToNearestAnnot <- distToNearest@elementMetadata$distance

  gwasResDF_temp <- gwasResDF

  colnames(gwasResDF_temp) <- paste(colnames(gwasResDF), "gwas", sep = "_")

  annotNearestToGWASranges <- cbind(gwasResDF_temp, annotNearestToGWASranges)

  rm(gwasResDF_temp)


  if(outFmt == "data.frame"){
    message("Returning data.frame")
    return(annotNearestToGWASranges)
  }
  if(outFmt == "GenomicRanges"){
    message("Returning GenomicRanges")
    annotNearestToGWASranges <- makeGRangesFromDataFrame(annotNearestToGWASranges, keep.extra.columns = T)
    return(annotNearestToGWASranges)
  }

}


##function for creating manhattan plot with closest gene annotation on significant gwas hits.
##requires annotated gwas object. basically, each SNP assigned an annotation in the same row
manhattan_annot_plot <- function(annotatedGWASresults, traitIDcol = "Trait_gwas",
                                 genomicSeqIDcol =  "seqnames", posIDcol = "Position_gwas",
                                 pValIDcol = "MarkerpVal_gwas" , sigTresholdIDcol = "SigTreshold_gwas",
                                 markerNameIDcol = "Marker_gwas", annotNameIDcol = "Gene", annotTypeIDcol = "annotType",
                                 annotDistanceIDcol = "distanceToNearestAnnot", annotStartIDcol = "start",
                                 annotEndIDcol = "end", labelType = "composite", colorPal = "Dark2",
                                 zoomGR = "none", pointColorIDCol = "Trait_gwas",
                                 annotateEffect = F, effectSizeIDcol = "none", labelAnnots = F, annotationLabelAngle = 90
                                 )
  {

  library(ggplot2)
  library(ggrepel)
  library(plyr)

  if(labelType %in% c("annotationName", "annotationDistance", "composite") == F){
    stop("label parameter must be one of 'annotationName', or 'composite'")
  }

  #create temporary data.frame with remapped variable IDs
  annotatedGWASres <- data.frame(trait =  annotatedGWASresults[,traitIDcol],
                                     marker = annotatedGWASresults[,markerNameIDcol],
                                     seqnames =  annotatedGWASresults[,genomicSeqIDcol],
                                     position = annotatedGWASresults[,posIDcol],
                                     sig = annotatedGWASresults[,pValIDcol],
                                     annotationName = annotatedGWASresults[,annotNameIDcol],
                                     annotStart = annotatedGWASresults[,annotStartIDcol],
                                     annotEnd = annotatedGWASresults[,annotEndIDcol],
                                     annotationDistance = annotatedGWASresults[,annotDistanceIDcol],
                                     annotType = annotatedGWASresults[,annotTypeIDcol],
                                     pointColor =  annotatedGWASresults[,pointColorIDCol]
                                     )

  #add default annotation type in case not set
  if(annotTypeIDcol == ""){
    annotatedGWASres$annotType <- "default"
  }

  #add signifcance treshold
  if(class(sigTresholdIDcol) == "character"){
    annotatedGWASres$sigTreshold = annotatedGWASresults[,sigTresholdIDcol]
  }
  if(class(sigTresholdIDcol) == "numeric"){
    annotatedGWASres$sigTreshold = sigTresholdIDcol
  }
  if(sigTresholdIDcol == ""){
    message("Setting default significance treshold to Bonferroni correction of 0.05")
    annotatedGWASres$sigTreshold = 0.05/length(unique(annotatedGWASres$marker))
  }

  if(max(annotatedGWASres$sigTreshold)>1){
    message("Significant treshold larger than 1. Transforming to 10^-treshold.")
    annotatedGWASres$sigTreshold <- 10^-annotatedGWASres$sigTreshold
  }


  #set labels for plotting based on significance
  annotatedGWASres$annotationName_label <- as.character(annotatedGWASres$annotationName)

  annotatedGWASres$annotationName_label[annotatedGWASres$sig > annotatedGWASres$sigTreshold] <- "" #mark not significant
  annotatedGWASres$annotationName_label <- as.factor(annotatedGWASres$annotationName_label)

  annotatedGWASres$annotationDistance_label <- annotatedGWASres$annotationDistance
  annotatedGWASres$annotationDistance_label[annotatedGWASres$sig > annotatedGWASres$sigTreshold] <- "" #mark not significant


  if(sum(annotatedGWASres$sig <= annotatedGWASres$sigTreshold) > 10) warning("More than 10 significant SNPs")

  annotatedGWASres$sigLog10 = -log10(annotatedGWASres$sig)
  annotatedGWASres$sigTresholdLog10 = -log10(annotatedGWASres$sigTreshold)

  annotatedGWASres_sig <- annotatedGWASres[annotatedGWASres$sig <= annotatedGWASres$sigTreshold,]

  if(labelType == "annotationName"){
    annotatedGWASres$annotation_label <- annotatedGWASres$annotationName_label
  }

  if(labelType == "annotationDistance"){
    annotatedGWASres$annotation_label <- annotatedGWASres$annotationDistance_label
  }

  if(labelType == "composite"){
    annotatedGWASres$annotation_label <- paste(annotatedGWASres$annotationName_label,
                                               annotatedGWASres$annotationDistance_label, sep = ": ")
    annotatedGWASres$annotation_label <- gsub("^: $", "" ,annotatedGWASres$annotation_label)
  }


  if(annotateEffect !=F & effectSizeIDcol != "none"){
    annotatedGWASres$effectSize <- annotatedGWASresults[,effectSizeIDcol]
    annotatedGWASres_sig <- annotatedGWASres[annotatedGWASres$sig <= annotatedGWASres$sigTreshold,]
    annotatedGWASres_sig$effectArrowPos <- annotatedGWASres_sig$position
    annotatedGWASres_sig$effect <- "Negative"
    for(idx in 1:nrow(annotatedGWASres_sig)){
      if(annotatedGWASres_sig$effectSize[idx]>0) {
        annotatedGWASres_sig$effectArrowStart[idx] <- annotatedGWASres_sig$sigLog10[idx] * 0.975
        annotatedGWASres_sig$effectArrowEnd[idx] <- annotatedGWASres_sig$sigLog10[idx] * 1.025
        annotatedGWASres_sig[idx]$effect <- "Positive"
      } else{
        annotatedGWASres_sig$effectArrowStart[idx] <- annotatedGWASres_sig$sigLog10[idx] * 1.025
        annotatedGWASres_sig$effectArrowEnd[idx] <- annotatedGWASres_sig$sigLog10[idx] * 0.975
      }
    }
    arrowSeg <- geom_segment(aes(x=effectArrowPos, xend=effectArrowPos, y=effectArrowStart,
      yend=effectArrowEnd), arrow=arrow(length=unit(abs(0.02*annotatedGWASres_sig$effectSize), "npc")),
      data = annotatedGWASres_sig)
  }



  if (class(zoomGR) == "character") if (zoomGR == "none"){
    g1 <- ggplot(annotatedGWASres) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom") +
      geom_point(aes(position, sigLog10, col = pointColor), alpha = 0.3) +
      geom_point(aes(position, sigLog10, col = pointColor), data = annotatedGWASres_sig)  +
      geom_hline(aes(yintercept = sigTresholdLog10, col = trait), lty = "dashed") +
      scale_color_brewer(palette = colorPal) +
      xlab("Position by Chr") + ylab("-Log10Pval") +
      geom_label_repel(aes(x = position, y = sigLog10, label = annotation_label),
                       inherit.aes = T, box.padding = 2, hjust = 0.5, vjust = 0.5) +
      facet_wrap(~seqnames, scales = "free_x")

    if(annotateEffect !=F & effectSizeIDcol != "none"){
      g1 <- g1 + arrowSeg
    }

  }

  if(class(zoomGR) == "character") if (zoomGR == "auto"){
    zoomDF <- annotatedGWASres[annotatedGWASres$annotationName %in% unique(annotatedGWASres_sig$annotationName),]

    g1 <- ggplot(zoomDF) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom") +
      geom_rect(aes(xmin = annotStart, xmax = annotEnd, ymin=0, ymax = max(sig)/5,
                    fill = annotType), data = zoomDF[zoomDF$sig <= zoomDF$sigTreshold,], alpha = 0.5) +
      geom_point(aes(position, sigLog10, col = pointColor), alpha = 0.3) +
      geom_point(aes(position, sigLog10, col = pointColor), data = zoomDF[zoomDF$sig <= zoomDF$sigTreshold,])  +
      geom_hline(aes(yintercept = sigTresholdLog10, col = trait), lty = "dashed") +
      scale_color_brewer(palette = colorPal) +
      xlab("Position by Chr") + ylab("-Log10Pval") +
      geom_label_repel(aes(x = position, y = sigLog10, label = annotation_label),
                       inherit.aes = T, box.padding = 2, hjust = 0.5, vjust = 0.5) +
      facet_wrap(~paste(seqnames, annotationName, sep = ": "), scales = "free_x")
    if(annotateEffect !=F & effectSizeIDcol != "none"){
      g1 <- g1 + arrowSeg
    }

    if(labelAnnots){
      labelData <- ddply(zoomDF, .(annotationName), summarize,
                         xPos = mean(c(unique(annotStart),unique(annotEnd))),
                         sig=max(sig), seqnames=unique(seqnames))
      nAnnots <- nrow(labelData)
      if(nAnnots>20){
        message(paste(nAnnots, "annotations will be labeled, it might take a bit..."))
      }
      g1 <- g1 + geom_text_repel(aes(x = xPos, y = max(sigLog10, na.rm = T)/5,
                                     label = unique(annotationName)), angle = annotationLabelAngle,
                                  inherit.aes = T, box.padding = 2, data = labelData, hjust = 0.5, vjust = 0.5)
    }

  }


  if(class(zoomGR) == "GRanges"){
    annotatedGWASresGR <- makeGRangesFromDataFrame(annotatedGWASres, keep.extra.columns = T)
    zoomDF <- as.data.frame(subsetByOverlaps(annotatedGWASresGR, zoomGR, ignore.strand = T))

    if(annotateEffect !=F & effectSizeIDcol != "none"){
      zoomDF <- merge(zoomDF, annotatedGWASres_sig[,c("marker", "effectArrowPos", "effectArrowStart", "effectArrowEnd")],
                      by = "marker", all.x=T)
    }

    g1 <- ggplot(zoomDF) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom") +
      geom_rect(aes(xmin = start, xmax = end, ymin = 0, ymax = max(sigLog10, na.rm = T)/5,
                    fill = annotType), alpha = 0.5) +
      geom_point(aes(position, sigLog10, col = pointColor), alpha = 0.3) +
      geom_point(aes(position, sigLog10, col = pointColor), data = zoomDF[zoomDF$sig <= zoomDF$sigTreshold,])  +
      geom_hline(aes(yintercept = sigTresholdLog10, col = trait), lty = "dashed") +
      scale_color_brewer(palette = colorPal) +
      xlab("Position by Chr") + ylab("-Log10Pval") +
      geom_label_repel(aes(x = position, y = sigLog10, label = annotation_label, hjust = 0.5, vjust = 0.5),
                       inherit.aes = T, box.padding = 2) +
      facet_wrap(~paste(seqnames, sep = ": "), scales = "free_x")

    if(annotateEffect !=F & effectSizeIDcol != "none"){
      arrowSeg <- geom_segment(aes(x=effectArrowPos, xend=effectArrowPos, y=effectArrowStart,
                                   yend=effectArrowEnd), arrow=arrow(length=unit(0.02, "npc")),
                               data = zoomDF[zoomDF$sig <= zoomDF$sigTreshold,])
      g1 <- g1 + arrowSeg
    }

    if(labelAnnots){
      labelData <- ddply(zoomDF, .(annotationName), summarize,
                         xPos = mean(c(unique(start),unique(end))),
                         sig=max(sig), seqnames=unique(seqnames))
      nAnnots <- nrow(labelData)
      if(nAnnots>20){
        message(paste(nAnnots, "annotations will be labeled, it might take a bit..."))
      }
      g1 <- g1 + geom_text_repel(aes(x = xPos, y = max(sigLog10, na.rm = T)/5,
                                     label = unique(annotationName)), inherit.aes = T,
                                 box.padding = 2, data = labelData, angle = annotationLabelAngle,
                                 hjust = 0.5, vjust = 0.5)
    }

  }

  return(g1)

}
