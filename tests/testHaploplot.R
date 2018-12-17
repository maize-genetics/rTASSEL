library(ggplot2)
library(DBI)
library(GenomicRanges)

# Parameters
phg_path <- "/Users/tw493/Downloads/phgMaizeGzip_deltaG_depthAndNonMappedFix.db"
# ---

phg_db <- dbConnect(RSQLite::SQLite(), phg_path)
haps <- dbGetQuery(phg_db, 'SELECT gt.line_name, rr.chrom, rr.range_start, rr.range_end, rr.ref_range_id FROM haplotypes AS ht INNER JOIN genotypes AS gt ON ht.gamete_grp_id=gt.genoid INNER JOIN reference_ranges AS rr ON ht.ref_range_id=rr.ref_range_id')
dbDisconnect(phg_db)

haps$line_name <- as.factor(haps$line_name)
haps$chrom <- as.integer(haps$chrom)
haps$ref_range_id <- as.factor(haps$ref_range_id)
haps$y <- match(haps$line_name, levels(haps$line_name)) # generate y values for haplotypes based on taxa
haps$genic <- ifelse(as.integer(haps$ref_range_id) < 37804, "anchor", "interanchor")
haps <- haps[order(haps$range_start),] # order haps by ref range start

hap_g_ranges <- GRanges(
  seqnames = as.character(haps$chrom),
  ranges = IRanges(
    start = haps$range_start,
    end = haps$range_end,
    names = paste(paste0("rr", haps$ref_range_id), haps$line_name, sep = "_")
  ),
  strand = strand("+"),
  taxa = haps$line_name,
  y = haps$y
)
seqlengths(hap_g_ranges) <- aggregate(haps$range_end, by = list(haps$chrom), max)$x
genome(hap_g_ranges) <- "PHG_Maize"

window_g_range <- GRanges(
  seqnames = as.character(1),
  ranges = IRanges(
    start = 1,
    end = 10000,
    names = "haploplot_window"
  ),
  strand = strand("+")
)

y_breaks <- seq(levels(haps$line_name))
y_labels <- levels(haps$line_name)

plot_hap_ranges <- subsetByOverlaps(hap_g_ranges, window_g_range)
plot_haps <- data.frame(
  start = start(plot_hap_ranges@ranges),
  end = end(plot_hap_ranges@ranges),
  y = plot_hap_ranges@elementMetadata$y,
  line = plot_hap_ranges@elementMetadata$taxa
)

ggplot(plot_haps) +
  geom_rect(aes(xmin = pmax(start, 1), ymin = y - 0.25, xmax = pmin(end, 10000), ymax = y + 0.25), color = "gray") +
  geom_vline(xintercept = plot_haps$end[plot_haps$end <= 10000] + 0.5, alpha = 0.5) +
  scale_y_continuous(breaks = y_breaks, labels = y_labels, limits = c(0, length(y_labels) + 1)) +
  xlab("Position (bp)") +
  ylab("Taxa") +
  theme_minimal() +
  theme(legend.position = "none")
