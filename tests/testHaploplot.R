library(ggplot2)
library(GenomicRanges)
library(rJava)

# load rtassel
setwd(paste(fs::path_home(), "Code", "rtassel", sep = .Platform$file.sep))

# initialize java
.jinit(parameters="-Xmx6g")
path_jar_dir <- "inst/java"
.jaddClassPath(paste0(path_jar_dir, .Platform$file.sep, "sTASSEL.jar"))
.jaddClassPath(paste0(path_jar_dir, .Platform$file.sep, "phg.jar"))
.jaddClassPath(paste0(path_jar_dir, .Platform$file.sep, "kotlin-stdlib-1.3.10.jar"))
.jaddClassPath(paste0(path_jar_dir, .Platform$file.sep, "lib"))

haplotypeGraphBuilderPlugin <- function(configFile, myMethods) {
  plugin <- new(J("net.maizegenetics.pangenome.api.HaplotypeGraphBuilderPlugin"), .jnull(), FALSE)
  plugin$setParameter("configFile",toString(configFile))
  plugin$setParameter("methods",toString(myMethods))
  plugin$setParameter("includeSequences",toString(FALSE))
  plugin$build()
}

configFilePath <- "data/configSQLiteR.txt"
method <- "mummer4,*"
phg_hap_graph <- haplotypeGraphBuilderPlugin(configFilePath, method)
generateRforPHG <- .jnew("net.maizegenetics.pangenome/pipelineTests.GenerateRForPHG")
rr_vecs <- rJava::.jcall(generateRforPHG, "Lnet/maizegenetics/pangenome/pipelineTests/RefRangeVectors;", "graphToRefRangeVectors", phg_hap_graph, .jarray(integer()))
hap_vecs <- rJava::.jcall(generateRforPHG, "Lnet/maizegenetics/pangenome/pipelineTests/HaplotypesDataVectors;", "graphToHapsInRefRangeVectors", phg_hap_graph, .jarray(integer()), FALSE, FALSE)

rr_df <- data.frame(
  ref_range_id = rr_vecs$refRangeId,
  chrom = rr_vecs$chromosomes,
  range_start = rr_vecs$startPos,
  range_end = rr_vecs$endPos
)
hap_df <- data.frame(
  ref_range_id = hap_vecs$refRangeIds,
  line_name = hap_vecs$taxa
)

hap_df <- cbind(hap_df, rr_df[match(hap_df$ref_range_id, rr_df$ref_range_id), c("chrom", "range_start", "range_end")])
rm(rr_df, rr_vecs, hap_vecs, phg_hap_graph)

hap_df$chrom <- as.integer(hap_df$chrom)
hap_df$y <- match(hap_df$line_name, levels(hap_df$line_name)) # generate y values for haplotypes based on taxa
hap_df$anchor <- ifelse(as.integer(hap_df$ref_range_id) < 37804, "anchor", "interanchor")
hap_df <- hap_df[order(hap_df$range_start),] # order haps by ref range start

hap_g_ranges <- GRanges(
  seqnames = as.character(hap_df$chrom),
  ranges = IRanges(
    start = hap_df$range_start,
    end = hap_df$range_end,
    names = paste(paste0("rr", hap_df$ref_range_id), hap_df$line_name, sep = "_")
  ),
  strand = strand("+"),
  taxa = hap_df$line_name,
  y = hap_df$y,
  anchor = hap_df$anchor
)

window_start <- 1
window_end <- 100000
window_g_range <- GRanges(
  seqnames = as.character(1),
  ranges = IRanges(
    start = window_start,
    end = window_end,
    names = "haploplot_window"
  ),
  strand = strand("+")
)

y_breaks <- seq(levels(hap_df$line_name))
y_labels <- levels(hap_df$line_name)
rm(hap_df)

plot_hap_ranges <- subsetByOverlaps(hap_g_ranges, window_g_range)
plot_haps <- data.frame(
  start = start(plot_hap_ranges@ranges),
  end = end(plot_hap_ranges@ranges),
  y = plot_hap_ranges@elementMetadata$y,
  line = plot_hap_ranges@elementMetadata$taxa,
  anchor = plot_hap_ranges@elementMetadata$anchor
)

ggplot(plot_haps) +
  geom_rect(aes(xmin = pmax(start, window_start), ymin = y - 0.25, xmax = pmin(end, window_end), ymax = y + 0.25, fill = anchor)) +
  geom_vline(xintercept = plot_haps$end[plot_haps$end <= window_end] + 0.5, alpha = 0.5) +
  scale_fill_brewer(palette = 4, type = "qual") +
  scale_x_continuous(limits = c(window_start, window_end), labels = scales::comma) +
  scale_y_continuous(breaks = y_breaks, labels = y_labels, limits = c(0, length(y_labels) + 1)) +
  xlab("Position (bp)") +
  ylab("Taxa") +
  theme_minimal()
