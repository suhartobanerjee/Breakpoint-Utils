source("/fast/groups/ag_sanders/work/projects/suharto/breakpoint_analysis/R/methods.R")
source("/fast/groups/ag_sanders/work/projects/suharto/digital_karyotype/R/utils.R",
    chdir = T
)


################################################################################
# USER CHANGE REQUIRED
#
# 1. raw file path
raw_file_path <- "../proc/breakpointR_concat.tsv.gz"
# 2. Sample name (if required)
sample_name <- "P1530"
# 3. save_dir for the plots.
save_path <- "digKar_bins_freq"
#
################################################################################




################################################################################
# FUNCTIONS
#
# to choose which interval to consider.
# big by default
.Get_Interval_Dt <- function(raw_dt, interval = "big") {
    interval_dt <- raw_dt[, .(
        chrom,
        get(paste0(interval, "_start")),
        get(paste0(interval, "_end")),
        strand_state,
        sample,
        cell_name
    )]
    setnames(
        interval_dt,
        c("V2", "V3"),
        c("start_loc", "end_loc")
    )


    return(interval_dt)
}



PreProc_Data <- function(raw_dt, interval = "big") {
    interval_dt <- .Get_Interval_Dt(raw_dt = raw_dt)

    interval_gro <- makeGRangesFromDataFrame(interval_dt,
        seqnames.field = "chrom",
        start.field = "start_loc",
        end.field = "end_loc",
        keep.extra.columns = T
    )
    interval_overlap_dt <- PercentOverlap(bin_genome,
        interval_gro,
        bin_size = bin_size
    )



    take_high_dups <- interval_overlap_dt[perc_overlap < 100,
        .SD[order(-perc_overlap)][2, idx],
        by = .(cell_name, strand_state),
    ][, V1]
    interval_overlap_dt <- interval_overlap_dt[!take_high_dups]


    return(interval_overlap_dt)
}





Calculate_Freq <- function(interval_overlap_dt) {
    # calculating the number of bps in a bin
    bin_freq <- data.table(table(big_overlap[, bin_id]))
    setnames(
        bin_freq,
        c("V1", "N"),
        c("bin_id", "freq")
    )
    bin_freq[, bin_id := as.integer(bin_id)]


    # calc the n_cells
    cell_contrib <- big_overlap[, length(unique(cell_name)),
        by = bin_id
    ]
    setnames(
        cell_contrib,
        "V1",
        "n_cells"
    )
    cell_contrib[, n_cells := as.integer(n_cells)]


    # merging cell_contrib and bin_freq
    bin_freq_ncells <- merge(bin_freq,
        cell_contrib,
        by = "bin_id"
    )


    return(bin_freq_ncells)
}




Merge_Freq_Dt <- function(bin_freq_ncells, bin_gen_dt) {
    # merging to get bin start and end
    # replacing NA with 0
    big_freq_ncells_verbose <- merge.data.table(bin_gen_dt,
        bin_freq_ncells,
        by = "bin_id",
        all = T
    )
    na_bins <- big_freq_ncells_verbose[is.na(big_freq_ncells_verbose[, freq])][, `:=`(freq = 0, n_cells = 0)][, bin_id]
    big_freq_ncells_verbose[
        bin_id %in% na_bins,
        `:=`(
            freq = 0,
            n_cells = 0
        )
    ]


    return(big_freq_ncells_verbose)
}


# plotting the freq and ncells.
Plot_Frequency <- function(freq_dt, dig_kar_h2) {
    freq_counter <- 1
    dig_kar_counter <- 2
    freq_plot_list <- c()
    for (chr in seq_along(all_chr)) {
        freq_plot_list[[freq_counter]] <- ggplot() +
            geom_rect(
                data = freq_dt[chrom == all_chr[chr]],
                aes(
                    ymin = KAR_YMIN,
                    ymax = freq,
                    xmin = start_loc,
                    xmax = end_loc,
                    fill = n_cells
                )
            ) +
            facet_grid(rows = vars(chrom), switch = "y") +
            scale_x_continuous(limits = CHR_PROP_XLIM_H2) +
            scale_y_continuous(
                limits = c(
                    min(freq_dt[, freq]),
                    max(freq_dt[, freq])
                ),
                breaks = seq(min(freq_dt[, freq]),
                    max(freq_dt[, freq]),
                    by = 10
                )
            ) +
            scale_fill_gradientn(
                colors = hsv(1, seq(0, 1, length.out = length(unique(freq_dt$n_cells))), 1),
                limits = c(min(freq_dt$n_cells), max(freq_dt$n_cells)),
                breaks = c(seq(
                    min(freq_dt[, freq]),
                    max(freq_dt[, freq])
                ))
            ) +
            theme_dig_kar +
            theme(
                plot.margin = margin(
                    t = 0,
                    r = 0,
                    b = 0,
                    l = 0,
                    unit = "pt"
                ),
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
################################################################################


################################################################################
# raw file path here
raw_dt <- fread(raw_file_path)


# pre-process.
ProcFile(raw_dt)


# ignore this if sample is not required.
raw_dt <- raw_dt[sample == sample_name]


# binning genome, set bin size
bin_size <- 1e6
ret_list <- BinGenome(bin_size = bin_size)
bin_genome <- ret_list[["bin_genome"]]
bin_gen_dt <- ret_list[["bin_gen_dt"]]


interval_overlap_dt <- PreProc_Data(raw_dt)
bin_freq_ncells <- Calculate_Freq(interval_overlap_dt)
big_freq_ncells_verbose <- Merge_Freq_Dt(bin_freq_ncells, bin_gen_dt)
big_freq_ncells_verbose




# plotting skeleton digKar
dig_kar_h2 <- generate_skeleton_h2(ideo_dt,
    chr_name_size = 10
)
dig_kar_h2 <- add_centromere(
    dig_kar = dig_kar_h2,
    centro_dt = centro_dt
)



freq_plot_list <- Plot_Frequency(big_freq_ncells_verbose, dig_kar_h2)



# add the sample name here
dig_kar <- stitch_skeleton_h2(
    dig_kar_h2 = freq_plot_list,
    haplo_label = paste0(
        "Breakpoint Distribution: ",
        sample_name
    ),
    subplot_scale = rep(c(3, 1), 48)
)



comb_legend <- generate_combined_legend(
    combined_haplo_calls = big_freq_ncells_verbose,
    fill_param = "n_cells",
    legend_title = "n_cells",
    gradient = T
)


comb_plot <- plot_grid(comb_legend,
    dig_kar,
    ncol = 2,
    nrow = 1,
    rel_widths = MASTER_PLOT_LEGEND_REL_WIDTHS_HALF,
    align = "hv",
    axis = "l"
)


# add the sample name in cell_name
# change the name of save_dir. Plots
# will be saved to ../plots/{save_dir}
save_digital_karyotype(comb_plot,
    save_dir = save_path,
    cell_name = paste0("SV_BP_", sample_name)
)
################################################################################
