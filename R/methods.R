ProcFile <- function(raw_dt) {
    setnames(raw_dt,
             c("seqnames", 
               "CI.start", 
               "CI.end",
               "start", 
               "end", 
               "genoT"
               ),
             c("chrom",
               "big_start",
               "big_end",
               "small_start",
               "small_end",
               "strand_state"
             )
    )

    raw_dt[, `:=`(geno1 = strsplit(strand_state, "-")[[1]][1],
                  geno2 = strsplit(strand_state, "-")[[1]][2])]

    raw_dt[, sample := substr(filenames,
                                 start = 0,
                                 stop = 5)]
    raw_dt[, cell_name := substr(filenames,
                                 start = 0,
                                 stop = 10)]
    #     raw_dt[, cell_name := tstrsplit(filenames, "\\.", keep = 1)]

    setcolorder(raw_dt,
                c("chrom",
                  "big_start",
                  "big_end",
                  "small_start",
                  "small_end",
                  "geno1",
                  "geno2",
                  "sample",
                  "cell_name")
    )



    return(raw_dt)
}


BinGenome <- function(bin_size) {

    ideo_dt <- fread("../../digital_karyotype/proc/ideogram_scaffold.tsv")
    ideo_width <- ideo_dt[, width]
    names(ideo_width) <- ideo_dt[, chrom]
    ideo_width


    bin_genome <- tileGenome(ideo_width,
                             tilewidth = bin_size,
                             cut.last.tile.in.chrom = T
    )
    bin_genome$bin_id <- 1:length(bin_genome)

    bin_gen_dt <- as.data.table(bin_genome)
    bin_gen_dt[, strand := NULL]
    setnames(bin_gen_dt,
             c("seqnames", "start", "end"),
             c("chrom", "start_loc", "end_loc")
    )

    return(list(bin_genome = bin_genome,
                bin_gen_dt = bin_gen_dt
        )
    )
}


PercentOverlap <- function(bin_genome, calls_gro, bin_size) {

    overlaps <- findOverlaps(bin_genome,
                             calls_gro
    )

    calls_dt <- as.data.table(calls_gro)
    bin_gen_dt <- as.data.table(bin_genome)

    calls_dt <- calls_dt[subjectHits(overlaps)]
    calls_dt[, bin_id := bin_gen_dt[queryHits(overlaps), bin_id]]


    setnames(bin_gen_dt,
             c("seqnames", "start", "end"),
             c("chrom", "start_loc", "end_loc")
    )
    setnames(calls_dt,
             c("seqnames", "start", "end"),
             c("chrom", "start_loc", "end_loc")
    )

    for (entry in 1:nrow(calls_dt)) {
        curr_bin <- calls_dt[entry, bin_id]


        # for completely embedded
        if(calls_dt[entry, start_loc] - bin_gen_dt[bin_id == curr_bin, start_loc] > 0 &
           bin_gen_dt[bin_id == curr_bin, end_loc] - calls_dt[entry, end_loc] > 0) {
            perc_overlap_temp <- 100

        } else {

            # Left shifted
            if((bin_gen_dt[bin_id == curr_bin, end_loc]-
                calls_dt[entry, end_loc]) > 0) {

                perc_overlap_temp <- round((calls_dt[entry, end_loc]-
                                      bin_gen_dt[bin_id == curr_bin, start_loc]) / bin_size * 100, 2)
            }


            # Right shifted
            if((bin_gen_dt[bin_id == curr_bin, end_loc]-
                calls_dt[entry, end_loc]) < 0) {

                perc_overlap_temp <- round((bin_gen_dt[bin_id == curr_bin, end_loc]-
                                      calls_dt[entry, start_loc]) / bin_size * 100, 2)
            }
        }


#         if (perc_overlap_temp > 100) {
#             perc_overlap_temp  <- 100
#         }


        calls_dt[entry, perc_overlap := perc_overlap_temp]
    }


    calls_dt[, idx := 1:nrow(calls_dt)]
    setcolorder(calls_dt,
                "idx"
    )


    return(calls_dt)
}

