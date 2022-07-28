library("cbmr")
library("limma")
library("edgeR")
library("ggplot2")
library("data.table")
library("stringr")
library("org.Mm.eg.db")

#### Hypothalamus
hyp <- prepare_rna_dgelist(hypothalamus$counts,
                           hypothalamus$metadata,
                           sample_col = "Sample ID (Library ID)")
hyp_first <- prepare_rna_dgelist(hypothalamus_first$counts,
                                 hypothalamus_first$metadata,
                                 sample_col = "Sample ID (Library ID)")

hyp$samples$run <- "B"
hyp_first$samples$run <- "A"

y <- cbind(hyp, hyp_first)

y <- y[, y$samples$Sample.ID..Library.ID. != "0135A_2"]

y$samples$Condition.2 <- fifelse(y$samples$Condition.2 == "3 hours post-injection",
                                  "3hours", "4weeks")

y$samples$group <- factor(paste0(y$samples$Condition.1, "_",
                                 y$samples$Condition.2),
                          levels = c("Saline_3hours",
                                     "Psilocybin_3hours",
                                     "Saline_4weeks",
                                     "Psilocybin_4weeks"))

names(y$samples)[18] <- "Extraction_pool"
design <- model.matrix(~ 0 + group + Extraction_pool, data = y$samples)

idx_expressed <- filterByExpr(y, design = design)
y <- y[idx_expressed, ]

colour_scale <- scale_color_manual(name = NULL,
                                   values = c("Saline_3hours" = "#a6cee3",
                                              "Psilocybin_3hours" = "#1f78b4",
                                              "Saline_4weeks" = "#b2df8a",
                                              "Psilocybin_4weeks" = "#33a02c"),
                                   labels = c("Saline_3hours" = "Saline, 3 hours",
                                              "Psilocybin_3hours" = "Psilocybin, 3 hours",
                                              "Saline_4weeks" = "Saline, 4 weeks",
                                              "Psilocybin_4weeks" = "Psilocybin, 4 weeks")
)
p <- ggplot_mds(y, dim_plot = c(1, 2),
                colour_by = "group",
                col_scale = colour_scale)
p

cpm_values <- cpm(y, log = TRUE)

my_grp <- as.character(c(1:31,1:31))

cpm_corrected <- removeBatchEffect(cpm_values, batch = my_grp,
                                   design = model.matrix(~ run, data = y$samples))

ggplot_mds(cpm_corrected, dim_plot = c(3, 4),
           colour_by = y$samples$group, col_scale = colour_scale)



###### VTA
pfc <- prepare_rna_dgelist(pfc$counts,
                           pfc$metadata,
                           sample_col = "Sample ID (Library ID)")
pfc_first <- prepare_rna_dgelist(pfc_first$counts,
                                 pfc_first$metadata,
                                 sample_col = "Sample ID (Library ID)")

pfc$samples$run <- "B"
pfc_first$samples$run <- "A"

y <- cbind(pfc, pfc_first)

y <- y[, y$samples$Sample.ID..Library.ID. != "0135B_6"]

y$samples$Condition.2 <- fifelse(y$samples$Condition.2 == "3 hours post-injection",
                                 "3hours", "4weeks")

y$samples$group <- factor(paste0(y$samples$Condition.1, "_",
                                 y$samples$Condition.2),
                          levels = c("Saline_3hours",
                                     "Psilocybin_3hours",
                                     "Saline_4weeks",
                                     "Psilocybin_4weeks"))

names(y$samples)[18] <- "Extraction_pool"
design <- model.matrix(~ 0 + group + Extraction_pool, data = y$samples)

idx_expressed <- filterByExpr(y, design = design)
y <- y[idx_expressed, ]

colour_scale <- scale_color_manual(name = NULL,
                                   values = c("Saline_3hours" = "#a6cee3",
                                              "Psilocybin_3hours" = "#1f78b4",
                                              "Saline_4weeks" = "#b2df8a",
                                              "Psilocybin_4weeks" = "#33a02c"),
                                   labels = c("Saline_3hours" = "Saline, 3 hours",
                                              "Psilocybin_3hours" = "Psilocybin, 3 hours",
                                              "Saline_4weeks" = "Saline, 4 weeks",
                                              "Psilocybin_4weeks" = "Psilocybin, 4 weeks")
)
p <- ggplot_mds(y, dim_plot = c(1, 2),
                colour_by = "group",
                col_scale = colour_scale)
p

cpm_values <- cpm(y, log = TRUE)

my_grp <- as.character(c(1:31,1:31))

cpm_corrected <- removeBatchEffect(cpm_values, batch = my_grp,
                                   design = model.matrix(~ run, data = y$samples))

ggplot_mds(cpm_corrected, dim_plot = c(3, 4),
           colour_by = y$samples$group, col_scale = colour_scale)
