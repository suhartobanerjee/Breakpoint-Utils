library(data.table)
library(ggplot2)
library(GenomicRanges)


source("./methods.R")
source("../../digital_karyotype/R/utils.R",
       chdir = T
)


# raw file
raw_dt <- fread("../proc/breakpointR_concat.tsv.gz")
raw_dt



# binning genome
bin_size = 1e6
ret_list <- BinGenome(bin_size = bin_size)
bin_genome <- ret_list[["bin_genome"]]
bin_gen_dt <- ret_list[["bin_gen_dt"]]


bin_genome
bin_gen_dt



Get_Interval_Dt <- function(raw_dt, interval = "big") {

    big_dt <- raw_dt[, .(chrom,
                         get(paste0(interval, "_start")),
                         get(paste0(interval, "_end")),
                         strand_state,
                         sample,
                         cell_name)]
    setnames(big_dt,
             c("V2", "V3"),
             c("start_loc", "end_loc")
    )


    return(big_dt)
}


big_dt <- Get_Interval_Dt(raw_dt = raw_dt)
big_dt


big_gro <- makeGRangesFromDataFrame(big_dt,
                                    seqnames.field = "chrom",
                                    start.field = "start_loc",
                                    end.field = "end_loc",
                                    keep.extra.columns = T
)
big_gro



big_overlap <- PercentOverlap(bin_genome,
                              big_gro,
                              bin_size = bin_size) 
                                    
big_overlap


take_high_dups <- big_overlap[perc_overlap < 100,
                             .SD[order(-perc_overlap)][2, idx],
                             by = .(cell_name, strand_state),
                             ][, V1]
big_overlap <- big_overlap[!take_high_dups]




bin_freq <- data.table(table(big_overlap[, bin_id]))
setnames(bin_freq,
         c("V1", "N"),
         c("bin_id", "freq")
)
bin_freq[, bin_id := as.integer(bin_id)]
bin_freq



cell_contrib <- big_overlap[, length(unique(cell_name)), 
                            by = bin_id]
setnames(cell_contrib,
         "V1",
         "n_cells"
)
cell_contrib[, n_cells := as.integer(n_cells)]
big_overlap[1:12]

# merging cell_contrib and bin_freq
bin_freq_ncells <- merge(bin_freq,
                         cell_contrib,
                         by = "bin_id"
)
bin_freq_ncells



big_freq_ncells_verbose <- merge.data.table(bin_gen_dt,
                                 bin_freq_ncells,
                                 by = "bin_id",
                                 all = T
)
na_bins <- big_freq_ncells_verbose[is.na(big_freq_ncells_verbose[, freq])
                        ][,`:=`(freq = 0, n_cells = 0)
                        ][, bin_id]
big_freq_ncells_verbose[bin_id %in% na_bins,
                        `:=`(freq = 0,
                             n_cells = 0)]


big_freq_ncells_verbose[, mid_point := (start_loc + width / 2 - 1)]
big_freq_ncells_verbose[, ymin := 0]
str(big_freq_ncells_verbose)


big_freq_ncells_verbose
big_freq_ncells_verbose[chrom == "chr21",
                        max(end_loc)]
ideo_dt[chrom == "chr21",
        max(end_loc)]






# trying out line plot with freq
# tmp_plot <- "../plots/tmp/tmp_plot.pdf"
# 
# freq_plot <- ggplot(data = big_freq_ncells_verbose[chrom == "chr1"],
#                     aes(xmin = start_loc,
#                         xmax = end_loc,
#                         ymin = ymin,
#                         ymax = freq)
#                     ) +
#                                  geom_rect() +
#                                  coord_fixed(ratio = 300000,
#                                              xlim = CHR_PROP_XLIM_H2) +
#                                  facet_grid(rows = vars(chrom), switch = "y") +
#                                  theme_dig_kar
# 
# freq_plot2 <- ggplot(data = big_freq_ncells_verbose[chrom == "chr2"],
#                     aes(xmin = start_loc,
#                         xmax = end_loc,
#                         ymin = ymin,
#                         ymax = freq)
#                     ) +
#                                  geom_rect() +
#                                  coord_fixed(ratio = 300000,
#                                              xlim = CHR_PROP_XLIM_H2) +
#                                  facet_grid(rows = vars(chrom), switch = "y") +
#                                  theme_dig_kar
# 
# freq_plot3 <- ggplot(data = big_freq_ncells_verbose[chrom == "chr3"],
#                     aes(xmin = start_loc,
#                         xmax = end_loc,
#                         ymin = ymin,
#                         ymax = freq)
#                     ) +
#                                  geom_rect() +
#                                  coord_fixed(ratio = 300000,
#                                              xlim = CHR_PROP_XLIM_H2) +
#                                  facet_grid(rows = vars(chrom), switch = "y") +
#                                  theme_dig_kar
# 
# freq_plot4 <- ggplot(data = big_freq_ncells_verbose[chrom == "chr4"],
#                     aes(xmin = start_loc,
#                         xmax = end_loc,
#                         ymin = ymin,
#                         ymax = freq)
#                     ) +
#                                  geom_rect() +
#                                  coord_fixed(ratio = 300000,
#                                              xlim = CHR_PROP_XLIM_H2) +
#                                  facet_grid(rows = vars(chrom), switch = "y") +
#                                  theme_dig_kar
# ggsave("../plots/digKar_bins_freq/freq_plot.pdf")




# plotting gradient digKar
centro_dt <- Load_Centro_Ranges(centro_path = "../../digital_karyotype/proc/centromeres_concise.tsv.gz")
centro_dt <- Proc_Centro_Trapezium(centro_dt)
centro_dt


dig_kar_h2 <- generate_skeleton_h2(ideo_dt[chrom %in% all_chr[10:14]])
dig_kar_h2 <- generate_skeleton_h2(ideo_dt)
dig_kar_h2 <- add_centromere(dig_kar = dig_kar_h2,
                             centro_dt = centro_dt)
dig_kar <- stitch_skeleton_h2(dig_kar_h2)
save_digital_karyotype(dig_kar,
                       save_dir = "digKar_bins_freq",
                       cell_name = "ideo"
)

# comb <- c()
# 
# comb <- plot_grid(comb,
#                   freq_plot,
#                   dig_kar_h2[["chr1"]],
#                   ncol = 1,
#                   align = "v",
#                   axis = "l"
# )
# 
# comb <- plot_grid(comb,
#                   freq_plot2,
#                   dig_kar_h2[["chr2"]],
#                   ncol = 1,
#                   align = "v",
#                   axis = "l"
# )
# 
# comb <- plot_grid(comb,
#                   freq_plot3,
#                   dig_kar_h2[["chr3"]],
#                   ncol = 1,
#                   align = "v",
#                   axis = "l"
# )
# 
# comb <- plot_grid(comb,
#                   freq_plot4,
#                   dig_kar_h2[["chr4"]],
#                   ncol = 1,
#                   align = "v",
#                   axis = "l"
# )
# save_digital_karyotype(comb,
#                        save_dir = "digKar_bins_freq",
#                        cell_name = "chr_freq"
# )
# 
# ggsave("../plots/digKar_bins_freq/chr_freq.pdf")
# 


################################################################################
Plot_Frequency <- function(freq_dt) {

    freq_counter <- 1
    dig_kar_counter <- 2
    freq_plot_list <- c()
    for (chr in seq_along(all_chr)) {

        freq_plot_list[[freq_counter]] <- ggplot() +
                       geom_rect(data = freq_dt[chrom == all_chr[chr]],
                            aes(ymin = KAR_YMIN,
                                ymax = freq,
                                xmin = start_loc,
                                xmax = end_loc,
                                fill = n_cells)) +
#                        xlim(CHR_PROP_XLIM_H2) +
#                        coord_fixed(ratio = CHR_PROP_RATIO - 1e7,
#                                    xlim = CHR_PROP_XLIM_H2) +
                       facet_grid(rows = vars(chrom), switch = "y") +
#                        coord_cartesian(
#                                        ylim = c(min(freq_dt[, freq]),
#                                                 max(freq_dt[, freq])),
#                                        expand = F) +
                       scale_x_continuous(limits = CHR_PROP_XLIM_H2) +
                       scale_y_continuous(limits = c(min(freq_dt[, freq]),
                                                     max(freq_dt[, freq])),
                                          breaks = seq(min(freq_dt[, freq]),
                                                       max(freq_dt[, freq]),
                                                       by = 20)) +
                       scale_fill_gradientn(colors = hsv(1, seq(0, 1, length.out = length(unique(freq_dt$n_cells))), 1),
                                            limits = c(min(freq_dt$n_cells), max(freq_dt$n_cells)),
                                            breaks = c(seq(min(freq_dt[, freq]),
                                                       max(freq_dt[, freq])))) +
                       theme_dig_kar + 
                       theme(plot.margin = margin(t = 0,
                                                  r = 0,
                                                  b = 0,
                                                  l = 0,
                                                  unit = "pt"),
                             strip.text.y.left = element_blank(),
                             axis.line.y.left = element_line(),
                             axis.ticks.y = element_line(),
                             axis.text.y.left = element_text(),
                             legend.position = "none"
                       )
        

        freq_plot_list[[dig_kar_counter]] <- dig_kar_h2[[all_chr[chr]]]

        freq_counter <- freq_counter + 2
        dig_kar_counter <- dig_kar_counter + 2
    }


    return(freq_plot_list)
}



freq_plot_list <- Plot_Frequency(big_freq_ncells_verbose)
# dig_kar_h2[["chr1"]]
# freq_plot_list



list_plot <- plot_grid(plotlist = freq_plot_list,
          align = "v",
          axis = "l",
          ncol = 1,
          rel_heights = rep(c(3,1), 24)
)

comb_legend <- generate_combined_legend(combined_haplo_calls = big_freq_ncells_verbose,
                                        fill_param = "n_cells",
                                        legend_title = "n_cells",
                                        gradient = T

)

comb_plot <- plot_grid(comb_legend,
                       list_plot,
                       ncol = 2,
                       rel_widths = MASTER_PLOT_LEGEND_REL_WIDTHS_HALF,
                       align = "hv",
                       axis = "l"
)

# 
save_digital_karyotype(comb_plot,
                       save_dir = "digKar_bins_freq",
                       cell_name = "list_plot"
)
# ggsave("../plots/digKar_bins_freq/chr_freq.pdf",
#        width = 40,
#        height = 40,
#        units = "cm"
# )




# dig_kar <- stitch_skeleton_h2(freq_plot_list,
#                               haplo_label = "Breakpoint Distribution")
# save_digital_karyotype(dig_kar,
#                        save_dir = "digKar_bins_freq",
#                        cell_name = "list_plot_stitch"
# )

# for super large output
save_plot(filename = "../plots/digKar_bins_freq/chr_freq_cowplot.pdf",
          plot = dig_kar,
          base_height = 4,
          base_asp = 20,
          nrow = 48,
          limitsize = F
)


# ggsave("../plots/digKar_bins_freq/chr_freq.pdf",
#        width = 40,
#        height = 40,
# )




################################################################################
# dig_kar_h2 <- generate_info_layer_gradient_h2(dig_kar_h2,
#                                      big_overlap_10mb,
#                                      fill_param = "freq",
# )
# dig_kar_h2 <- stitch_skeleton_h2(dig_kar_h2,
#                                  haplo_label = "Frequency of Breakpoints in 10MB bins"
# )
# 
# 
# # create combined legend for the cell
# combined_legend <- generate_combined_legend(combined_haplo_calls = big_overlap_10mb,
#                                             fill_param = "freq",
#                                             legend_title = "Frequency",
#                                             gradient = T
# )
# 
# 
# dig_kar <- plot_grid(combined_legend,
#                      dig_kar_h2,
#                      ncol = 2,
#                      rel_widths = MASTER_PLOT_LEGEND_REL_WIDTHS_HALF,
#                      align = "hv",
#                      axis = "l")
# 
# 
# save_digital_karyotype(dig_kar,
#                        "../plots/digKar_bins_freq",
#                        paste0("P3254", bin_size)
# )
