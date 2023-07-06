library(data.table)
library(ggplot2)
library(GenomicRanges)
library(foreach)
library(parallel)


source("./methods.R")



brk_17A <- fread("../proc/martina/breakPointSummary_17_A.txt")
ProcFile(brk_17A)
brk_17A


sv_calls <- fread("../proc/martina/stringent_filterTRUE_17_A.tsv")
sv_calls


sv_calls[!cell %in% brk_17A[, cell_name]]
brk_17A <- brk_17A[cell_name %in% sv_calls[, cell]]
brk_17A


sv_calls[, length(unique(cell))]
brk_17A[, length(unique(cell_name))]
all_cell <- brk_17A[, unique(cell_name)]

brk_17A[cell_name == all_cell[1]]
sv_calls[cell == all_cell[1]]




sv_calls[cell == "P1467_i473"]
brk_17A[cell_name == "P1467_i473"]


brk_gro <- makeGRangesFromDataFrame(brk_17A,
                                    start.field = "big_start",
                                    end.field = "big_end",
                                    keep.extra.columns = T
)
brk_gro


brk_gro <- makeGRangesFromDataFrame(brk_17A,
                                    start.field = "small_start",
                                    end.field = "small_end",
                                    keep.extra.columns = T
)
brk_gro

sv_gro <- makeGRangesFromDataFrame(sv_calls,
                                   keep.extra.columns = T
)
sv_gro





hits <- findOverlaps(brk_gro,
             sv_gro
)


unique(queryHits(hits))

brk_17A
sce <- brk_17A[!unique(queryHits(hits))]
nrow(sce)

foreach(c = seq_along(all_cell)) %do% {
    print(c)
}


sce_big <- NULL

filter_sv_overlaps <- function(brk_dt, sv_dt) {

    sv_dt <- sv_dt[cell %in% brk_dt[, cell_name]]
    brk_dt <- brk_dt[cell_name %in% sv_dt[, cell]]

    all_cells <- brk_dt[, unique(cell_name)]



    sce_big <- data.table()
    for(ncell in seq_along(all_cells)) {


        brk_big_gro <- makeGRangesFromDataFrame(brk_dt[cell_name == all_cells[ncell]],
                                                start.field = "big_start",
                                                end.field = "big_end",
                                                keep.extra.columns = T
        )


        sv_gro <- makeGRangesFromDataFrame(sv_dt[cell == all_cells[ncell]],
                                           keep.extra.columns = T
        )


        hits <- findOverlaps(brk_big_gro,
                                  sv_gro)
        if(length(queryHits(hits)) > 0) {

            sce_big <- rbind(sce_big,
                             brk_dt[cell_name == all_cells[ncell]
                              ][!unique(queryHits(hits))])
        }
    } |>
    system.time()



    return(sce_big)

}



files <- list.files("../proc/marti_files/")
files
break_files <- files[grep("breakPointSummary", files)]
break_files


# strsplit(strsplit(files, "\\.")[[1]], "_")[[1]][2]
# paste0(strsplit(strsplit(files, "\\.")[[1]], "_")[[1]][2],
#        "_",
#        strsplit(strsplit(files, "\\.")[[1]], "_")[[1]][3]
# )
# 
# files[grep("17_A", files)[2]]
# 
# paste0(strsplit(strsplit(break_files[3], "\\.")[[1]], "_")[[1]][2],
#        "_",
#        strsplit(strsplit(break_files[3], "\\.")[[1]], "_")[[1]][3]
# )
# 
# nfile
#     f_id <- paste0(strsplit(strsplit(file, "\\.")[[1]], "_")[[1]][2],
#                    "_",
#                    strsplit(strsplit(file, "\\.")[[1]], "_")[[1]][3]
#     )
# f_id
#     r_id <- grep(f_id, files)
# r_id
# 
#     brk_dt <- fread(paste0("../proc/martina/",
#                            files[r_id[1]]))


cluster <- makeForkCluster(16)
doParallel::registerDoParallel(cluster)

for(file in break_files) {

    f_id <- paste0(strsplit(strsplit(file, "\\.")[[1]], "_")[[1]][2],
                   "_",
                   strsplit(strsplit(file, "\\.")[[1]], "_")[[1]][3]
    )
    r_id <- grep(f_id, files)


    brk_dt <- fread(paste0("../proc/marti_files/",
                           files[r_id[1]]))
    sv_dt <- fread(paste0("../proc/marti_files/",
                          files[r_id[2]]))

    ProcFile(brk_dt)
#     brk_dt[, cell_name := substr(cell_name,
#                                  start = 0,
#                                  stop = 10)]


    ret_list <- filter_sv_overlaps(brk_dt, sv_dt)


    fwrite(ret_list,
           paste0("../proc/sv_filter/sv_filtered_big_", 
                  f_id,
                  ".tsv"
           ),
           sep = "\t"
    )

#     fwrite(ret_list[["sce_small"]],
#            paste0("../proc/sv_filter/sv_filtered_small_", 
#                   f_id,
#                   ".tsv"
#            ),
#            sep = "\t"
#     )
}

# stopCluster(cluster)
warnings()

sv_dt
brk_dt
file
all_cells


all_cell <- brk_dt[, unique(cell_name)]
all_cell


b <- brk_dt[cell_name == all_cell[1]]
s <- sv_dt[cell == all_cell[1]]

b
s

brk_dt <- brk_17A
ProcFile(brk_dt)
brk_dt
b_gro <- makeGRangesFromDataFrame(brk_dt[cell_name == all_cell[1]],
                                                start.field = "big_start",
                                                end.field = "big_end",
                                                keep.extra.columns = T
        )

sv_gro <- makeGRangesFromDataFrame(sv_dt[cell == all_cell[1]],
                                   keep.extra.columns = T
) 

b_gro
sv_gro


h <- findOverlaps(b_gro, sv_gro)
queryHits(h) |>
length()

queryHits(h)

brk_dt[cell_name == all_cell[1]]

brk_dt[cell_name == all_cell[1]
       ][!queryHits(h)]




brk_dt[cell_name == all_cell[1] &
       chrom == "chr1"]
sv_dt[cell == all_cell[1] &
      chrom == "chr1"]

files[r_id]

substr("P1530_i402_", 
       start = 0,
       stop = 10
)





# corr plot
filtered <- list.files("../proc/sv_filter/")

bind_sce <- data.table()
for(file in filtered) {
    f <- fread(paste0("../proc/sv_filter/",
                      file))
    bind_sce <- rbind(bind_sce,
                       f
    )
}

bind_sce
bind_sce <- bind_sce[chrom != "chrX" |
         chrom != "chrY"]

sv_files <- grep("TRUE", 
                 list.files("../proc/marti_files/"),
                 value = T
)
sv_files

sv_bind <- data.table()
for(file in sv_files) {
    f <- fread(paste0("../proc/marti_files/",
                      file))
    sv_bind <- rbind(sv_bind,
                       f
    )
}
sv_bind
sv_bind <- sv_bind[!(chrom == "chrX" |
        chrom == "chrY")]


setnames(bind_sce,
         "cell_name",
         "cell"
)
bind_sce


corr_dt <- merge(sv_bind[, .N, by = cell],
      bind_sce[, .N, by = cell],
      by = "cell"
)

setnames(corr_dt,
         c("N.x", "N.y"),
         c("sv", "sce")
)
corr_dt
corr_dt[, sample := substr(cell,
                           start = 0,
                           stop = 5)]
corr_dt




corr_plot <- ggplot(data = corr_dt,
       aes(x = sv,
           y = sce,
           fill = sample
       )) +
    geom_point(aes(color = factor(sample))) +
    geom_smooth(method='lm', formula= y~x)

ggsave("../plots/sv_sce/all_together.pdf")







# individually
id <- "17_D"
sv_dt <- fread(paste0("../proc/marti_files/stringent_filterTRUE_", 
                      id,
                      ".tsv"))
sce_dt <- fread(paste0("../proc/sv_filter/sv_filtered_big_",
                       id,
                       ".tsv"))

sce_dt
setnames(sce_dt,
         "cell_name",
         "cell"
)




corr_dt <- merge(sv_dt[, .N, by = cell],
      sce_dt[, .N, by = cell],
      by = "cell"
)

setnames(corr_dt,
         c("N.x", "N.y"),
         c("sv", "sce")
)
corr_dt


corr_plot <- ggplot(data = corr_dt,
       aes(x = sv,
           y = sce
       )) +
    geom_point() +
    geom_smooth(method='lm', formula= y~x) +
    ggtitle(paste0("SV vs SCE for cell: ",
                    id))


ggsave(paste0("../plots/sv_sce/",
              id,
              ".pdf"),
       corr_plot
)



