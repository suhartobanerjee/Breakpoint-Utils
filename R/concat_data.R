library(data.table)


dir <- list.files("../proc/marti_files/")
dir


sctrip <- grep("filter", dir, value = T)
sctrip


brkpnt <- dir[!(dir %in% sctrip)]


sctrip_files <- data.table()
brkpnt_files <- data.table()
for (file in dir) {
    if (file %in% sctrip) {
        sctrip_files <- rbind(
            sctrip_files,
            fread(paste0(
                "../proc/marti_files/",
                file
            ))
        )
    }

    if (file %in% brkpnt) {
        brkpnt_files <- rbind(
            brkpnt_files,
            fread(paste0(
                "../proc/marti_files/",
                file
            ))
        )
    }
}

# sanity check
sctrip_files
sctrip_files[, unique(sample)]

source("./methods.R")
ProcFile(brkpnt_files)
brkpnt_files
brkpnt_files[, unique(sample)]


# adding sample col to match breakpoint R dt
setnames(
    sctrip_files,
    c("sample", "cell"),
    c("sample_patient", "cell_name")
)
sctrip_files[, sample := substr(cell_name,
    start = 0,
    stop = 5
)]
sctrip_files



fwrite(sctrip_files,
    "../proc/sctrip_concat.tsv.gz",
    sep = "\t"
)

fwrite(brkpnt_files,
    "../proc/breakpointR_concat.tsv.gz",
    sep = "\t"
)
