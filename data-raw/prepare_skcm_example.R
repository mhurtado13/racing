## ------------------------------------------------------------------
## Prepare the skcm_example dataset shipped with the RaCInG package.
##
## Source files live in example_data/ and correspond to SKCM melanoma
## data that has already been preprocessed through prepare_input_files().
##   celltype_ligand.csv            → Lmatrix
##   celltype_receptor.csv          → Rmatrix
##   SKCM_TMEmod_cell_fractions.csv → Cmatrix
##   SKCM_LRpairs_weights_min.csv   → LRmatrix
##
## This script reads them through the same helpers used by
## generateInput(), subsets to 10 patients, and saves the result as
## data/skcm_example.rda.
## ------------------------------------------------------------------

library(dplyr)
library(tibble)
library(tidyr)

## ---- helpers (mirrors the package functions) ----

sortPermute <- function(stringlist) {
  ord <- order(stringlist)
  list(sortlist = stringlist[ord], perm = ord)
}

## ---- read Lmatrix ----
lig_df <- read.csv("example_data/celltype_ligand.csv", row.names = 1)
ligands   <- colnames(lig_df)
celltypes_lig <- rownames(lig_df)

Lmatrix <- lig_df %>%
  rownames_to_column("Celltype") %>%
  select(-Celltype) %>%
  mutate(across(everything(), as.numeric)) %>%
  as.matrix()

sp <- sortPermute(celltypes_lig)
celltypes <- sp$sortlist
Lmatrix   <- Lmatrix[sp$perm, , drop = FALSE]

## ---- read Rmatrix ----
rec_df <- read.csv("example_data/celltype_receptor.csv", row.names = 1)
receptors   <- colnames(rec_df)
celltypes_rec <- rownames(rec_df)

Rmatrix <- rec_df %>%
  rownames_to_column("Celltype") %>%
  select(-Celltype) %>%
  mutate(across(everything(), as.numeric)) %>%
  as.matrix()

sp_r <- sortPermute(celltypes_rec)
Rmatrix <- Rmatrix[sp_r$perm, , drop = FALSE]

## ---- read Cmatrix (cell-type fractions) ----
ctype_df <- read.csv("example_data/SKCM_TMEmod_cell_fractions.csv", row.names = 1)
celltypes_c <- colnames(ctype_df)
datasets    <- rownames(ctype_df)

Cmatrix <- ctype_df %>%
  rownames_to_column("Datasets") %>%
  select(-Datasets) %>%
  mutate(across(everything(), as.numeric)) %>%
  as.matrix()

# Normalize rows to sum to 1
row_sums <- rowSums(Cmatrix)
Cmatrix  <- Cmatrix / row_sums

# Sort cell types alphabetically
sp_c <- sortPermute(celltypes_c)
celltypes_c <- sp_c$sortlist
Cmatrix     <- Cmatrix[, sp_c$perm, drop = FALSE]

# Merge M1 + M2 → M (matches the package logic)
if (length(celltypes) != length(celltypes_c)) {
  idx_M1 <- which(celltypes_c == "M1")
  idx_M2 <- which(celltypes_c == "M2")
  if (length(idx_M1) > 0 && length(idx_M2) > 0) {
    Cmatrix[, idx_M1] <- Cmatrix[, idx_M1] + Cmatrix[, idx_M2]
    Cmatrix <- Cmatrix[, -idx_M2, drop = FALSE]
    celltypes_c[idx_M1] <- "M"
    celltypes_c <- celltypes_c[-idx_M2]
  }
}

## ---- read LRmatrix (interaction weights) ----
lr_df <- read.csv("example_data/SKCM_LRpairs_weights_min.csv", row.names = 1)
interactions <- colnames(lr_df)

distr <- lr_df %>%
  rownames_to_column("Datasets") %>%
  select(-Datasets) %>%
  mutate(across(everything(), ~ replace_na(as.numeric(.), 0))) %>%
  as.matrix()

# Split interaction names into ligand and receptor parts
parts <- strsplit(interactions, "_")
rowWordIndex <- sapply(parts, `[`, 1)
colWordIndex <- sapply(parts, `[`, 2)

row_idx <- match(rowWordIndex, ligands)
col_idx <- match(colWordIndex, receptors)

# Build 3D tensor: ligands × receptors × patients
DconnectionTensor <- array(0, dim = c(length(ligands), length(receptors), nrow(distr)))
for (i in seq_len(nrow(distr))) {
  DconnectionTensor[cbind(row_idx, col_idx, i)] <- distr[i, ]
}

# Normalize each patient slice to sum to 1
normsum <- apply(DconnectionTensor, 3, sum)
for (i in seq_along(normsum)) {
  if (normsum[i] != 0) {
    DconnectionTensor[, , i] <- DconnectionTensor[, , i] / normsum[i]
  }
}

## ---- Subset to 10 patients ----
n_example <- 10
Cmatrix_sub   <- Cmatrix[seq_len(n_example), , drop = FALSE]
LRmatrix_sub  <- DconnectionTensor[, , seq_len(n_example), drop = FALSE]
patient_names <- datasets[seq_len(n_example)]

## ---- Sign matrix (zeros = unknown) ----
Sign_matrix <- matrix(0, nrow = length(ligands), ncol = length(receptors))

## ---- Assemble the list ----
skcm_example <- list(
  Lmatrix     = Lmatrix,
  Rmatrix     = Rmatrix,
  Cmatrix     = Cmatrix_sub,
  LRmatrix    = LRmatrix_sub,
  celltypes   = celltypes,
  ligands     = ligands,
  receptors   = receptors,
  Sign_matrix = Sign_matrix
)

## ---- Save ----
usethis::use_data(skcm_example, overwrite = TRUE, compress = "xz")

cat("Done. skcm_example saved with", n_example, "patients,",
    length(ligands), "ligands,", length(receptors), "receptors,",
    length(celltypes), "cell types.\n")
