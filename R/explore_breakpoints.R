library(data.table)
library(ggplot2)


source("./methods.R")
source("../../digital_karyotype/R/utils.R")


# raw file
raw_dt <- fread("../data/breakPointSummary.txt")
raw_dt


ProcFile(raw_dt)
raw_dt
raw_dt[, length(unique(strand_state))]

brk_color <- hcl.colors(raw_dt[, length(unique(strand_state))])
names(brk_color) <- raw_dt[, unique(strand_state)]
brk_color


big_dt <- raw_dt[, .(chrom,
                     big_start,
                     big_end,
                     strand_state,
                     filenames)]
setnames(big_dt,
         c("big_start", "big_end", "filenames"),
         c("start_loc", "end_loc", "cell_name")
)
big_dt
big_dt[1, cell_name]



str(plot_digital_karyotype)
PlotDigitalKaryotype(cell_sv_dt = big_dt,
                       fill_arg = "strand_state",
                       color_pal = brk_color,
                       legend_title_arg = "Strand Switch",
                       plot_both_haplotypes = F,
                       kar_label = "All break points across cells",
                       save_digital_karyotype = T,
                       save_dir = "all_cells"
)

    
