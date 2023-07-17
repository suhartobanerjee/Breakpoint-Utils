library(data.table)
library(ggplot2)
library(GenomicRanges)
# library(foreach)
# library(parallel)


source("./methods.R")



sv_dt <- fread("../proc/sctrip_concat.tsv.gz")
brk_dt <- fread("../proc/breakpointR_concat.tsv.gz")


filter_sv_overlaps <- function(breakpoint_dt, sctrip_dt) {

    sctrip_dt <- sctrip_dt[cell_name %in% breakpoint_dt[, cell_name]]
    breakpoint_dt <- breakpoint_dt[cell_name %in% sctrip_dt[, cell_name]]

    all_cells <- breakpoint_dt[, unique(cell_name)]



    sce_big <- data.table()
    for(ncell in seq_along(all_cells)) {


        brk_big_gro <- makeGRangesFromDataFrame(breakpoint_dt[cell_name == all_cells[ncell]],
                                                start.field = "big_start",
                                                end.field = "big_end",
                                                keep.extra.columns = T
        )


        sv_gro <- makeGRangesFromDataFrame(sctrip_dt[cell_name == all_cells[ncell]],
                                           keep.extra.columns = T
        )


        hits <- findOverlaps(brk_big_gro,
                                  sv_gro)
        if(length(queryHits(hits)) > 0) {

            sce_big <- rbind(sce_big,
                             breakpoint_dt[cell_name == all_cells[ncell]
                              ][!unique(queryHits(hits))])
        }
    }



    return(sce_big)

}


sce_big <- filter_sv_overlaps(breakpoint_dt = brk_dt,
                              sctrip_dt = sv_dt
)
sce_big


corr_dt <- merge(sv_dt[, .N, by = .(sample, cell_name)],
      sce_big[, .N, by = .(sample, cell_name)],
      by = c("sample","cell_name")
)
setnames(corr_dt,
         c("N.x", "N.y"),
         c("sv", "sce")
)
corr_dt




# corr plot
corr_plot <- ggplot(data = corr_dt,
       aes(x = sv,
           y = sce,
           color = sample,
           fill = sample
       )) +
    geom_point() +
    geom_smooth(method='lm', formula= y~x) +
    ggtitle("SV vs SCE")


ggsave("../plots/sv_sce/all_together.pdf",
       corr_plot
)



