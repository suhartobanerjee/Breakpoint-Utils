ProcFile <- function(raw_dt) {
    setnames(raw_dt,
             c("seqnames", 
               "start", 
               "end", 
               "CI.start", 
               "CI.end",
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


    setcolorder(raw_dt,
                c("chrom",
                  "big_start",
                  "big_end",
                  "small_start",
                  "small_end",
                  "geno1",
                  "geno2",
                  "filenames")
    )



    return(raw_dt)
}
