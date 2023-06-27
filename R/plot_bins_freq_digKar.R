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



# binning genome
bin_size = 10e6
ret_list <- BinGenome(bin_size = bin_size)
bin_genome <- ret_list[["bin_genome"]]
bin_gen_dt <- ret_list[["bin_gen_dt"]]


bin_genome
bin_gen_dt



big_dt <- raw_dt[, .(chrom,
                     big_start,
                     big_end,
                     strand_state,
                     cell_name)]
setnames(big_dt,
         c("big_start", "big_end"),
         c("start_loc", "end_loc")
)
big_gro <- makeGRangesFromDataFrame(big_dt,
                                    seqnames.field = "chrom",
                                    start.field = "start_loc",
                                    end.field = "end_loc",
                                    keep.extra.columns = T
)
big_gro



big_overlap_5mb <- PercentOverlap(bin_genome, big_gro, bin_size = bin_size)



rem_entry_idx <- big_overlap_10mb[perc_overlap < 100,
                                  .SD[order(-perc_overlap)][2, idx],
                                  by = .(cell_name, strand_state),
                                  ][, V1]
big_overlap_10mb <- big_overlap_10mb[!rem_entry_idx]



bin_freq_10mb <- data.table(table(big_overlap_10mb[, bin_id]))
setnames(bin_freq_10mb,
         c("V1", "N"),
         c("bin_id", "freq")
)
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
combined_legend <- generate_combined_legend(combined_haplo_calls = big_overlap_10mb,
                                            fill_param = "freq",
                                            legend_title = "Frequency",
                                            gradient = T
)


dig_kar <- plot_grid(combined_legend,
                     dig_kar_h2,
                     ncol = 2,
                     rel_widths = MASTER_PLOT_LEGEND_REL_WIDTHS_HALF,
                     align = "hv",
                     axis = "l")


save_digital_karyotype(dig_kar,
                       "../plots/digKar_bins_freq",
                       paste0("P3254", bin_size)
)
