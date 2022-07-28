library(data.table)
load(here::here("data/hyp_pfc.rda"))

setDT(hyp_pfc$metadata)
setDT(hyp_pfc$counts)

hyp_pfc$metadata[, geo_name := paste(stringr::str_replace_all(`Cell type`,
                                                              pattern = " ",
                                                              replacement = "-"),
                                     `Condition 1`,
                                     stringr::str_replace_all(`Condition 2`,
                                                          pattern = " ",
                                                          replacement = "-"),
                                     seq_len(.N),
                                     sep = "_"),
                 by = c("Cell type", "Condition 1", "Condition 2")]



writeClipboard(hyp_pfc$metadata[, paste0("cp ", `Sample ID (Library ID)`,
                                         "_*_first_R", c(1,2),
                                         "_001.fastq.gz ../geo/", geo_name,
                                         "_1_R", c(1,2), ".fastq.gz"),
                                by = "geo_name"]$V1)

writeClipboard(hyp_pfc$metadata[, paste0("cp ", `Sample ID (Library ID)`,
                                         "_*_R", c(1,2),
                                         "_001.fastq.gz ../geo/", geo_name,
                                         "_2_R", c(1,2), ".fastq.gz"),
                                by = "geo_name"]$V1)

setnames(hyp_pfc$counts, hyp_pfc$metadata$`Sample ID (Library ID)`, hyp_pfc$metadata$geo_name)
fwrite(hyp_pfc$counts, file = here::here("data-raw/counts.csv.gz"))

writexl::write_xlsx(hyp_pfc$metadata, path = here::here("data-raw/metadata_for_geo.xlsx"))



writeClipboard(hyp_pfc$metadata[, paste0(geo_name,
                                         "_1_R", 1, ".fastq.gz"),
                                by = "geo_name"]$V1)

writeClipboard(hyp_pfc$metadata[, paste0(geo_name,
                                         "_2_R", 1, ".fastq.gz"),
                                by = "geo_name"]$V1)
