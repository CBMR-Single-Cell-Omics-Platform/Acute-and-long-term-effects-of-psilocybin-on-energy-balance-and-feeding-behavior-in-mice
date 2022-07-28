## code to prepare `DATASET` dataset goes here
library(data.table)
library(stringr)
library(readxl)
library(cbmr)

prepare_function <- function(metadata, pattern, path = here::here("data-raw/featureCounts/")){
  metadata <- readxl::read_excel(metadata)
  data.table::setDT(metadata)
  data.table::setnames(metadata, stringr::str_replace_all(colnames(metadata), "\\s", " "))

  featurecounts_files <- list.files(path = path,
                                    pattern = pattern, full.names = TRUE)

  names(featurecounts_files) <- stringr::str_remove(basename(featurecounts_files),
                                           "_S[[:digit:]]+_Aligned.sortedByCoord.out.bam_featureCounts.txt.gz"
  )
  names(featurecounts_files) <- stringr::str_replace(names(featurecounts_files), "-", "_")
  featurecounts_files <- featurecounts_files[metadata$`Sample ID (Library ID)`]

  counts <- cbmr::prepare_featurecounts(featurecounts_files, regex = "fq2count_raw/02_star/|_S[[:digit:]]+_Aligned.sortedByCoord.out.bam")
  data.table::setnames(counts, stringr::str_replace(colnames(counts), "-", "_"))

  list(metadata = metadata,
       counts = counts)
}

hypothalamus <- prepare_function(
  metadata = here::here("data-raw/Copy of 0135_metadata_A (003) with sample information updates 240521.xlsx"),
  pattern = "0135A.*.txt.gz"
)

pfc <- prepare_function(
  metadata = here::here("data-raw/Copy of 0135_metadata_B (002) with sample information updates 240521.xlsx"),
  pattern = "0135B.*.txt.gz"
)


hypothalamus_first <- prepare_function(
  metadata = here::here("data-raw/Copy of 0135_metadata_A (003) with sample information updates 240521.xlsx"),
  pattern = "0135A.*.txt.gz",
  path = here::here("data-raw/featureCounts_first_round/")
)

pfc_first <- prepare_function(
  metadata = here::here("data-raw/Copy of 0135_metadata_B (002) with sample information updates 240521.xlsx"),
  pattern = "0135B.*.txt.gz",
  path = here::here("data-raw/featureCounts_first_round/")
)

all(hypothalamus$counts$Geneid == hypothalamus_first$counts$Geneid)
all(pfc$counts$Geneid == pfc_first$counts$Geneid)
# See other script that there are no systematic differences between the runs
idx_hyp <- colnames(hypothalamus$counts)[7:38]
hypothalamus$counts <- cbind(
  hypothalamus$counts[, 1:6],
  hypothalamus$counts[, .SD, .SDcols = idx_hyp] +
  hypothalamus_first$counts[, .SD, .SDcols = idx_hyp]
)

idx_pfc <- colnames(pfc$counts)[7:38]
pfc$counts <- cbind(
  pfc$counts[, 1:6],
  pfc$counts[, .SD, .SDcols = idx_pfc] +
    pfc_first$counts[, .SD, .SDcols = idx_pfc]
)

### To analyse all together we merge the two experiments.
all(pfc$counts$Geneid == hypothalamus$counts$Geneid)

# The extraction pools are different between the tissues.
pfc$metadata[, `Extraction pool (A,B,C)`:=fifelse(`Extraction pool (A,B,C)` == "A", "C", "D")]

hypothalamus$metadata$`RIN value` <- 0
pfc$metadata$`RIN value` <- 0

metadatas <- list(
  hypothalamus = hypothalamus$metadata,
  pfc = pfc$metadata,
)

metadatas$hypothalamus[, tissue:="Hypothalamus"]
metadatas$pfc[, tissue:="PFC"]

common_cols <- Reduce(intersect, lapply(metadatas, colnames))
metadata <- rbindlist(lapply(metadatas, `[`, , ..common_cols))

common_cols <- Reduce(intersect, lapply(metadatas[c("hypothalamus", "pfc")], colnames))
metadata <- rbindlist(lapply(metadatas[c("hypothalamus", "pfc")], `[`, , ..common_cols))

hyp_pfc <- list(
  counts = cbind(hypothalamus$counts, pfc$counts[, 7:38]),
  metadata = metadata
)

usethis::use_data(hyp_pfc, overwrite = TRUE)
