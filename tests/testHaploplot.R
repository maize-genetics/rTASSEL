library(ggplot2)
library(DBI)

# Parameters
phg_path <- "/Users/tw493/Downloads/phgMaizeGzip_deltaG_depthAndNonMappedFix.db"
genic_maize_id_cutoff <- 37804L
chr <- 2L
num_ref_ranges_to_plot <- 100L
# ---

phg_db <- dbConnect(RSQLite::SQLite(), phg_path)
min_ref_range_id <- as.integer(dbGetQuery(phg_db, 'SELECT MIN(ref_range_id) FROM reference_ranges WHERE chrom == :chr',
                               params = list(chr = chr)))
max_ref_range_id <- as.integer(dbGetQuery(phg_db, 'SELECT MAX(ref_range_id) FROM reference_ranges WHERE chrom == :chr AND ref_range_id < :cutoff',
                               params = list(chr = chr, cutoff = genic_maize_id_cutoff)))
haps <- dbGetQuery(phg_db, 'SELECT haplotypes.ref_range_id, genotypes.line_name FROM haplotypes INNER JOIN genotypes ON haplotypes.gamete_grp_id=genotypes.genoid WHERE haplotypes.ref_range_id >= :min AND haplotypes.ref_range_id <= :max',
                   params = list(min = min_ref_range_id, max = max_ref_range_id))
dbDisconnect(phg_db)

filtered_haps <- haps[which(haps$ref_range_id < (min_ref_range_id + num_ref_ranges_to_plot)),]
filtered_haps$ref_range_id <- as.factor(filtered_haps$ref_range_id)
filtered_haps$line_name <- as.factor(filtered_haps$line_name)

ggplot(filtered_haps) +
  geom_point(aes(x = ref_range_id, y = line_name), shape = 15, size = 1, color = "gray") +
  ggtitle(paste0("Haploplot (Chr ",  chr, ")")) +
  xlab("Reference Range") +
  ylab("Taxa") +
  theme_minimal() +
  theme(legend.position = "none", axis.text.x = element_blank())
