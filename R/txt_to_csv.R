# Load required library
library(data.table)

Read_Sim_Output <- function(filename) {

  # ------------------------------------------------------------
  # Purpose:
  #   Read simulation output (.out file) for any communication type
  #   (D, W, TT, CT, GSCC)
  #
  # The function parses:
  #   - Metadata
  #   - Raw counts (average + std)
  #   - Composition (average + std)
  #
  # Output format adapts automatically depending on type:
  #   D    → 2D (cell x cell)
  #   W/TT/CT → 3D (cell x cell x cell)
  #   GSCC → 1D (cell)
  # ------------------------------------------------------------

  # Read full file
  lines <- read.csv(filename, header = FALSE, stringsAsFactors = FALSE)

  idx <- 1

  # -----------------------------
  # First line: communication type
  # -----------------------------
  triangle_type <- as.character(lines[idx, 1])
  idx <- idx + 1

  # -----------------------------
  # Second line: metadata
  # Format: cancer_type, weight_type, NoPat, N, itNo, av
  # -----------------------------
  meta <- lines[idx, ]
  idx <- idx + 1

  cancer_type <- meta[1]
  weight_type <- meta[2]
  NoPat <- as.integer(meta[3])
  N <- as.integer(meta[4])
  itNo <- as.integer(meta[5])
  av <- as.numeric(meta[6])

  # -----------------------------
  # Determine dimensionality
  # -----------------------------
  if (triangle_type == "GSCC") {
    dim_type <- 1
  } else if (triangle_type == "D") {
    dim_type <- 2
  } else {
    dim_type <- 3
  }

  # Infer number of cells
  nCells <- ncol(lines[4,])

  # -----------------------------
  # Initialize storage
  # -----------------------------
  if (dim_type == 3) {
    data_a <- array(0, dim = c(nCells, nCells, nCells, NoPat))
    data_s <- array(0, dim = c(nCells, nCells, nCells, NoPat))
  } else if (dim_type == 2) {
    data_a <- array(0, dim = c(nCells, nCells, NoPat))
    data_s <- array(0, dim = c(nCells, nCells, NoPat))
  } else {
    data_a <- array(0, dim = c(nCells, NoPat))
    data_s <- array(0, dim = c(nCells, NoPat))
  }

  raw_count <- matrix(0, nrow = NoPat, ncol = 2)

  # -----------------------------
  # Read patient data
  # -----------------------------
  while (idx <= nrow(lines)) {

    p <- as.integer(lines[idx, 1]) # patient ID, the other columns correspond to N and avg
    idx <- idx + 1

    if (is.na(p)) break

    # Skip cell names line
    idx <- idx + 1

    # Raw count
    vals <- lines[idx, ]
    idx <- idx + 1

    raw_count[p, 1] <- as.numeric(vals[2])  # average
    raw_count[p, 2] <- as.numeric(vals[3])  # std

    # Skip "Composition - Average"
    idx <- idx + 1

    # -----------------------------
    # Read averages
    # -----------------------------
    if (dim_type == 3) {
      for (n in 1:(nCells^3)) {
        row <- lines[idx, ]; idx <- idx + 1
        i <- as.integer(row[1])
        j <- as.integer(row[2])
        k <- as.integer(row[3])
        val <- as.numeric(row[4])
        data_a[i, j, k, p] <- val
      }
    } else if (dim_type == 2) {
      for (n in 1:(nCells^2)) {
        row <- lines[idx, ]; idx <- idx + 1
        i <- as.integer(row[1])
        j <- as.integer(row[2])
        val <- as.numeric(row[3])
        data_a[i, j, p] <- val
      }
    } else {
      for (n in 1:nCells) {
        row <- lines[idx, ]; idx <- idx + 1
        i <- as.integer(row[1])
        val <- as.numeric(row[2])
        data_a[i, p] <- val
      }
    }

    # Skip "Composition - Std"
    idx <- idx + 1

    # -----------------------------
    # Read standard deviations
    # -----------------------------
    if (dim_type == 3) {
      for (n in 1:(nCells^3)) {
        row <- lines[idx, ]; idx <- idx + 1
        i <- as.integer(row[1])
        j <- as.integer(row[2])
        k <- as.integer(row[3])
        val <- as.numeric(row[4])
        data_s[i, j, k, p] <- val
      }
    } else if (dim_type == 2) {
      for (n in 1:(nCells^2)) {
        row <- lines[idx, ]; idx <- idx + 1
        i <- as.integer(row[1])
        j <- as.integer(row[2])
        val <- as.numeric(row[3])
        data_s[i, j, p] <- val
      }
    } else {
      for (n in 1:nCells) {
        row <- lines[idx, ]; idx <- idx + 1
        i <- as.integer(row[1])
        val <- as.numeric(row[2])
        data_s[i, p] <- val
      }
    }
  }

  # -----------------------------
  # Return output
  # -----------------------------
  return(list(
    data_avg = data_a,
    data_std = data_s,
    raw_count = raw_count
  ))
}

# Generate_normalised_count_csv: create normalized CSV for triangle/graphlet data
Generate_normalised_count_csv <- function(cancer_type,
                                          weight_type,
                                          triangle_type,      ### CHECK WHETHER WE NEED TO CONSIDER DIFFERENT TYPES OF TRIANGLE!!
                                          average = 15,
                                          noCells = 10000,
                                          folder = "Input_data_RaCInG",
                                          remove_direction = TRUE,
                                          output_folder = "Output_data_RaCInG") {

  # -----------------------------
  # Load input data
  # -----------------------------
  # generateInput loads cell types, patient info, ligand-receptor lists, etc.
  input_data <- generateInput(weight_type, cancer_type, folder = folder)
  celltypes <- input_data$celltypes
  patient_names <- rownames(input_data$Dtypes) ############################ (TO BE REMOVED!) Assuming patient names are row names in Dtypes; adjust if needed

  # Standardize cell type names
  celltypes[celltypes == "CD8+ T"] <- "CD8" #####(TO BE REMOVED!)
  nCells <- length(celltypes)  # number of cell types

  # -----------------------------
  # Read simulation outputs
  # -----------------------------
  # data_raw  : non-normalized simulation
  # data_norm : normalized (uniform) simulation
  data_raw  <- Read_Sim_Output(file.path(output_folder,
                                         sprintf("%s_%s_%s.out", cancer_type, triangle_type, average)))
  data_norm <- Read_Sim_Output(file.path(output_folder,
                                         sprintf("%s_%s_%s_norm.out", cancer_type, triangle_type, average)))

  # Extract arrays of averages and patient-level raw counts
  av      <- data_raw$data_avg
  avN     <- data_norm$data_avg
  summary <- data_raw$raw_count
  summaryN <- data_norm$raw_count

  # -----------------------------
  # Convert arrays → data frames
  # -----------------------------
  # This step converts multi-dimensional arrays (triangle, direct, or GSCC) into
  # data frames where each column corresponds to a labeled combination of cell types.
  non_unif_data <- list()
  unif_data     <- list()

  if (length(dim(av)) == 4) {  # Triangles (3D interactions) + patients
    for (i in 1:nCells) {
      for (j in 1:nCells) {
        for (k in 1:nCells) {
          label <- paste(celltypes[i], celltypes[j], celltypes[k], sep = "_")
          non_unif_data[[label]] <- av[i, j, k, ]
          unif_data[[label]] <- avN[i, j, k, ]
        }
      }
    }
  } else if (length(dim(av)) == 3) {  # Direct interactions (2D) + patients
    for (i in 1:nCells) {
      for (j in 1:nCells) {
        label <- paste(celltypes[i], celltypes[j], sep = "_")
        non_unif_data[[label]] <- av[i, j, ]
        unif_data[[label]] <- avN[i, j, ]
      }
    }
  } else {  # GSCC (1D) + patients
    for (i in 1:nCells) {
      label <- celltypes[i]
      non_unif_data[[label]] <- av[i, ]
      unif_data[[label]] <- avN[i, ]
    }
  }

  # Convert lists to data frames for easier manipulation
  df  <- as.data.frame(non_unif_data)
  dfN <- as.data.frame(unif_data)

  # -----------------------------
  # Remove direction (merge labels)
  # -----------------------------
  # Some triangle types (like W) are directional; this step merges counts for
  # interactions that are equivalent when ignoring direction.
  if (remove_direction) {
    new  <- data.frame(matrix(nrow = nrow(df), ncol = 0))
    newN <- data.frame(matrix(nrow = nrow(dfN), ncol = 0))

    for (colname in colnames(df)) {
      parts <- strsplit(colname, "_")[[1]]

      if (triangle_type == "W" && length(parts) == 3) {
        # Sort only first and last cell for W triangles
        sorted_parts <- sort(c(parts[1], parts[3]))
        parts[1] <- sorted_parts[1]
        parts[3] <- sorted_parts[2]
        new_label <- paste(parts, collapse = "_")
      } else {
        # Sort all parts for other types (fully direction-agnostic)
        new_label <- paste(sort(parts), collapse = "_")
      }

      # Accumulate counts for merged labels
      if (new_label %in% colnames(new)) {
        new[[new_label]] <- new[[new_label]] + df[[colname]]
        newN[[new_label]] <- newN[[new_label]] + dfN[[colname]]
      } else {
        new[[new_label]] <- df[[colname]]
        newN[[new_label]] <- dfN[[colname]]
      }
    }

    # Update arrays for normalization
    av  <- as.matrix(new)
    avN <- as.matrix(newN)
  }

  # -----------------------------
  # Normalization
  # -----------------------------
  # Replace zero baseline counts to avoid division by zero.
  zero_idx <- which(av == 0 & avN == 0, arr.ind = TRUE)
  if (nrow(zero_idx) > 0) {
    av[zero_idx]  <- 1
    avN[zero_idx] <- 1
  }

  # Replace any remaining zeros in baseline to 1 (to prevent division errors)
  avN[avN == 0] <- 1

  # Normalize counts: divide non-uniform by uniform counts
  Norm <- av / avN

  # -----------------------------
  # Remove patients where simulation failed
  # -----------------------------
  # Patients with zero total counts in raw or baseline simulations are excluded
  # because they indicate failed or empty simulations.
  delete_indices <- unique(c(which(summary[,1] == 0), which(summaryN[,1] == 0)))
  if (length(delete_indices) > 0) {
    Norm <- Norm[-delete_indices, , drop = FALSE]
    patients <- patient_names[-delete_indices]
  } else {
    patients <- patient_names
  }

  # -----------------------------
  # Convert normalized data to dataframe
  # -----------------------------
  df <- as.data.frame(Norm)
  rownames(df) <- patients

  # -----------------------------
  # Save CSV
  # -----------------------------
  write.csv(df,
            file = file.path(output_folder,
                             sprintf("%s_%s_%s_cells_%s_deg_data.csv", cancer_type, triangle_type, noCells, average)),
            row.names = TRUE)

  return(df)
}
