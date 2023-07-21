library(data.table)
library(ggplot2)
library(GenomicRanges)


source("./methods.R")
source("../../digital_karyotype/R/utils.R")


# raw file
raw_dt <- fread("../data/breakPointSummary.txt")
raw_dt


ProcFile(raw_dt)
raw_dt
# strsplit("P3254_i438.sort.mdup.bam.RData.1", "\\.")[[1]][1]
raw_dt[, (tstrsplit(filenames, "\\.", keep = 1))]
# raw_dt[, length(unique(strand_state))]

brk_color <- hcl.colors(raw_dt[, length(unique(strand_state))])
names(brk_color) <- raw_dt[, unique(strand_state)]
brk_color


big_dt <- raw_dt[, .(
    chrom,
    big_start,
    big_end,
    strand_state,
    cell_name
)]
setnames(
    big_dt,
    c("big_start", "big_end"),
    c("start_loc", "end_loc")
)
big_dt
big_dt[1, cell_name]



str(plot_digital_karyotype)
PlotDigitalKaryotype(
    cell_sv_dt = big_dt,
    fill_arg = "strand_state",
    color_pal = brk_color,
    legend_title_arg = "Strand Switch",
    plot_both_haplotypes = F,
    kar_label = "All break points across cells",
    save_digital_karyotype = T,
    save_dir = "all_cells"
)


# binning genome
bin_size <- 10e6
ret_list <- BinGenome(bin_size = bin_size)
bin_genome <- ret_list[["bin_genome"]]
bin_gen_dt <- ret_list[["bin_gen_dt"]]


bin_genome
# bin_gen_dt <- as.data.table(bin_genome)
# bin_gen_dt[, strand := NULL]
# setnames(bin_gen_dt,
#          c("seqnames", "start", "end"),
#          c("chrom", "start_loc", "end_loc")
# )
bin_gen_dt


big_dt
big_gro <- makeGRangesFromDataFrame(big_dt,
    seqnames.field = "chrom",
    start.field = "start_loc",
    end.field = "end_loc",
    keep.extra.columns = T
)
big_gro



# overlaps <- findOverlaps(bin_genome,
#                          big_gro
# )
# overlaps
#
#
# check <- big_dt
# check[subjectHits(overlaps),      bin_id := bin_gen_dt[queryHits(overlaps), bin_id]]
# check
#
#
# big_overlap <- big_dt[subjectHits(overlaps)]
# big_overlap[, bin_id := bin_gen_dt[queryHits(overlaps), bin_id]]
# big_overlap
# big_overlap_dt <- as.data.table(big_overlap)
# big_overlap_dt



big_overlap_5mb <- PercentOverlap(bin_genome, big_gro, bin_size = bin_size)
big_overlap_1mb <- PercentOverlap(bin_genome, big_gro, bin_size = bin_size)
big_overlap_10mb <- PercentOverlap(bin_genome, big_gro, bin_size = bin_size)


big_overlap_5mb
big_overlap_1mb
big_overlap_10mb




# to keep the entry which has
# higher perc between two bin_id
rem_entry_idx <- big_overlap_10mb[perc_overlap < 100,
    .SD[order(-perc_overlap)][2, idx],
    by = .(cell_name, strand_state),
][, V1]
# sanity check
big_overlap_10mb[8:11]
big_overlap_10mb <- big_overlap_10mb[!rem_entry_idx]







# making plotting dirs
tmp_plot <- "../plots/tmp/tmp_plot.pdf"
hist_plot_dir <- "../plots/bin_hist/"
if (!dir.exists(hist_plot_dir)) {
    dir.create(hist_plot_dir)
}




hist_1mb <- ggplot(
    data = big_overlap_1mb,
    aes(
        x = bin_id,
        fill = chrom
    )
) +
    geom_histogram(bins = nrow(bin_gen_dt)) +
    ggtitle("Counts of Breakpoints in bin size: 1MB")


ggsave(paste0(hist_plot_dir, "bins_hist_1mb.pdf"),
    plot = hist_1mb,
    device = "pdf",
    width = 29.7,
    height = 21,
    units = "cm"
)



hist_5mb <- ggplot(
    data = big_overlap_5mb,
    aes(
        x = bin_id,
        fill = chrom
    )
) +
    geom_histogram(bins = nrow(bin_gen_dt)) +
    ggtitle("Counts of Breakpoints in bin size: 5MB")


ggsave(paste0(hist_plot_dir, "bins_hist_5mb.pdf"),
    plot = hist_5mb,
    device = "pdf",
    width = 29.7,
    height = 21,
    units = "cm"
)



hist_10mb <- ggplot(
    data = big_overlap_10mb,
    aes(
        x = bin_id,
        fill = chrom
    )
) +
    geom_histogram(bins = nrow(bin_gen_dt)) +
    ggtitle("Counts of Breakpoints in bin size: 10MB")

ggsave(paste0(hist_plot_dir, "bins_hist_10mb.pdf"),
    plot = hist_10mb,
    device = "pdf",
    width = 29.7,
    height = 21,
    units = "cm"
)




big_overlap_10mb
bin_freq_10mb <- data.table(table(big_overlap_10mb[, bin_id]))
setnames(
    bin_freq_10mb,
    c("V1", "N"),
    c("bin_id", "freq")
)
bin_freq_10mb
typeof(bin_freq_10mb[, bin_id])
bin_freq_10mb[, bin_id := as.integer(bin_id)]




big_overlap_10mb <- merge(bin_gen_dt,
    bin_freq_10mb,
    by = "bin_id"
)
big_overlap_10mb



# plotting gradient digKar
dig_kar_h2 <- generate_skeleton_h2(ideo_dt)
dig_kar_h2 <- generate_info_layer_gradient_h2(dig_kar_h2,
    big_overlap_10mb,
    fill_param = "freq",
)
dig_kar_h2 <- stitch_skeleton_h2(dig_kar_h2,
    haplo_label = "Frequency of Breakpoints in 10MB bins"
)


# create combined legend for the cell
combined_legend <- generate_combined_legend(
    combined_haplo_calls = big_overlap_10mb,
    fill_param = "freq",
    legend_title = "Frequency",
    gradient = T
)


dig_kar <- plot_grid(combined_legend,
    dig_kar_h2,
    ncol = 2,
    rel_widths = MASTER_PLOT_LEGEND_REL_WIDTHS_HALF,
    align = "hv",
    axis = "l"
)


save_digital_karyotype(
    dig_kar,
    "../plots/digKar_bins_freq",
    "P3254"
)
