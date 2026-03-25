# Load required library
library(data.table)

# Triangle_Prop_Read: Reads triangle data from a .txt/.csv file
Triangle_Prop_Read <- function(filename) {
  # Read the entire file as a data frame
  # header = FALSE ensures no row is treated as header
  lines <- read.csv(filename, header = FALSE, stringsAsFactors = FALSE)
  
  idx <- 1  # index to keep track of current line
  
  # -----------------------------
  # First line: triangle type (e.g., "W")
  # -----------------------------
  triangle_type <- lines[idx, 1]
  idx <- idx + 1
  
  # -----------------------------
  # Second line: metadata
  # Format: cancer_type, weight_type, NoPat, N, itNo, av
  # -----------------------------
  meta <- lines[idx, ]
  idx <- idx + 1
  
  cancer_type <- meta[1]
  weight_type <- meta[2]
  NoPat <- as.integer(meta[3])  # number of patients
  N <- as.integer(meta[4])      # number of cells per graph
  itNo <- as.integer(meta[5])   # number of graphs generated
  av <- as.numeric(meta[6])     # average degree per graph
  
  # -----------------------------
  # Initialize arrays/matrices for storing the data
  # -----------------------------
  triangle_data_a <- array(0, dim = c(9, 9, 9, NoPat)) # average values
  triangle_data_s <- array(0, dim = c(9, 9, 9, NoPat)) # standard deviations
  triangle_raw_count <- matrix(0, nrow = NoPat, ncol = 2) # raw counts per patient
  
  # -----------------------------
  # Read patient data iteratively
  # -----------------------------
  while (idx <= nrow(lines)) {
    p <- lines[idx, 1]  # patient ID
    idx <- idx + 1
    
    if (is.na(p)) break  # stop if end of file reached
    
    # Skip a line (like in Python: next(reader))
    idx <- idx + 1
    
    # Next line: raw averages and std for this patient
    vals <- lines[idx, ]
    idx <- idx + 1
    triangle_raw_count[as.integer(p), 1] <- as.numeric(vals[2])
    triangle_raw_count[as.integer(p), 2] <- as.numeric(vals[3])
    
    # Skip another line
    idx <- idx + 1
    
    # -----------------------------
    # Read 9x9x9 array of averages
    # -----------------------------
    for (i in 1:9) {
      for (j in 1:9) {
        for (k in 1:9) {
          row <- lines[idx, ]
          idx <- idx + 1
          triangle_data_a[as.integer(row[1]), as.integer(row[2]), as.integer(row[3]), as.integer(p)] <- as.numeric(row[4])
        }
      }
    }
    
    # Skip a line before reading standard deviations
    idx <- idx + 1
    
    # -----------------------------
    # Read 9x9x9 array of standard deviations
    # -----------------------------
    for (i in 1:9) {
      for (j in 1:9) {
        for (k in 1:9) {
          row <- lines[idx, ]
          idx <- idx + 1
          triangle_data_s[as.integer(row[1]), as.integer(row[2]), as.integer(row[3]), as.integer(p)] <- as.numeric(row[4])
        }
      }
    }
  }
  
  # -----------------------------
  # Return a list equivalent to Python tuple
  # -----------------------------
  return(list(
    triangle_type = triangle_type,
    triangle_data_a = triangle_data_a,
    triangle_data_s = triangle_data_s,
    triangle_raw_count = triangle_raw_count,
    weight_type = weight_type,
    cancer_type = cancer_type,
    N = N,
    itNo = itNo,
    av = av
  ))
}

# Direct_Comm_Read: Reads direct communication data from a .txt/.csv file
Direct_Comm_Read <- function(filename) {
  # Read the entire file
  lines <- read.csv(filename, header = FALSE, stringsAsFactors = FALSE)
  
  idx <- 1  # keep track of current line
  
  # -----------------------------
  # First line: triangle type (e.g., "D")
  # -----------------------------
  triangle_type <- lines[idx, 1]
  idx <- idx + 1
  
  # -----------------------------
  # Second line: metadata
  # Format: cancer_type, weight_type, NoPat, N, itNo, av
  # -----------------------------
  meta <- lines[idx, ]
  idx <- idx + 1
  
  cancer_type <- meta[1]
  weight_type <- meta[2]
  NoPat <- as.integer(meta[3])  # number of patients
  N <- as.integer(meta[4])      # number of cells per graph
  itNo <- as.integer(meta[5])   # number of graphs generated
  av <- as.numeric(meta[6])     # average degree per graph
  
  # -----------------------------
  # Initialize arrays/matrices
  # triangle_data_a/s: 9 x 9 x NoPat
  # triangle_raw_count: 2 columns per patient
  # -----------------------------
  triangle_data_a <- array(0, dim = c(9, 9, NoPat))
  triangle_data_s <- array(0, dim = c(9, 9, NoPat))
  triangle_raw_count <- matrix(0, nrow = NoPat, ncol = 2)
  
  # -----------------------------
  # Read patient data iteratively
  # -----------------------------
  while (idx <= nrow(lines)) {
    p <- lines[idx, 1]  # patient ID
    idx <- idx + 1
    
    if (is.na(p)) break  # end of file
    
    # Skip one line
    idx <- idx + 1
    
    # Next line: raw averages and std for this patient
    vals <- lines[idx, ]
    idx <- idx + 1
    triangle_raw_count[as.integer(p), 1] <- as.numeric(vals[2])
    triangle_raw_count[as.integer(p), 2] <- as.numeric(vals[3])
    
    # Skip one line before reading averages
    idx <- idx + 1
    
    # -----------------------------
    # Read 9x9 array of averages
    # -----------------------------
    for (i in 1:9) {
      for (j in 1:9) {
        row <- lines[idx, ]
        idx <- idx + 1
        triangle_data_a[as.integer(row[1]), as.integer(row[2]), as.integer(p)] <- as.numeric(row[3])
      }
    }
    
    # Skip one line before reading std
    idx <- idx + 1
    
    # -----------------------------
    # Read 9x9 array of standard deviations
    # -----------------------------
    for (i in 1:9) {
      for (j in 1:9) {
        row <- lines[idx, ]
        idx <- idx + 1
        triangle_data_s[as.integer(row[1]), as.integer(row[2]), as.integer(p)] <- as.numeric(row[3])
      }
    }
  }
  
  # -----------------------------
  # Return a list equivalent to Python tuple
  # -----------------------------
  return(list(
    triangle_type = triangle_type,
    triangle_data_a = triangle_data_a,
    triangle_data_s = triangle_data_s,
    triangle_raw_count = triangle_raw_count,
    weight_type = weight_type,
    cancer_type = cancer_type,
    N = N,
    itNo = itNo,
    av = av
  ))
}

# GSCC_Read: Reads GSCC data from a .txt/.csv file
GSCC_Read <- function(filename) {
  # Read the entire file
  lines <- read.csv(filename, header = FALSE, stringsAsFactors = FALSE)
  
  idx <- 1  # current line index
  
  # -----------------------------
  # First line: triangle type (e.g., "GSCC")
  # -----------------------------
  triangle_type <- lines[idx, 1]
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
  # Initialize arrays/matrices
  # triangle_data_a/s: 9 x NoPat
  # triangle_raw_count: 2 columns per patient
  # -----------------------------
  triangle_data_a <- array(0, dim = c(9, NoPat))
  triangle_data_s <- array(0, dim = c(9, NoPat))
  triangle_raw_count <- matrix(0, nrow = NoPat, ncol = 2)
  
  # -----------------------------
  # Read patient data iteratively
  # -----------------------------
  while (idx <= nrow(lines)) {
    p <- lines[idx, 1]  # patient ID
    idx <- idx + 1
    
    if (is.na(p)) break  # end of file
    
    # Skip one line
    idx <- idx + 1
    
    # Next line: raw averages and std for this patient
    vals <- lines[idx, ]
    idx <- idx + 1
    triangle_raw_count[as.integer(p), 1] <- as.numeric(vals[2])
    triangle_raw_count[as.integer(p), 2] <- as.numeric(vals[3])
    
    # Skip one line before reading averages
    idx <- idx + 1
    
    # -----------------------------
    # Read 9 values for averages
    # -----------------------------
    for (i in 1:9) {
      row <- lines[idx, ]
      idx <- idx + 1
      triangle_data_a[as.integer(row[1]), as.integer(p)] <- as.numeric(row[2])
    }
    
    # Skip one line before reading standard deviations
    idx <- idx + 1
    
    # -----------------------------
    # Read 9 values for standard deviations
    # -----------------------------
    for (i in 1:9) {
      row <- lines[idx, ]
      idx <- idx + 1
      triangle_data_s[as.integer(row[1]), as.integer(p)] <- as.numeric(row[2])
    }
  }
  
  # -----------------------------
  # Return all data in a list
  # -----------------------------
  return(list(
    triangle_type = triangle_type,
    triangle_data_a = triangle_data_a,
    triangle_data_s = triangle_data_s,
    triangle_raw_count = triangle_raw_count,
    weight_type = weight_type,
    cancer_type = cancer_type,
    N = N,
    itNo = itNo,
    av = av
  ))
}

# Load required libraries
library(data.table)

# Generate_normalised_count_csv: create normalized CSV for triangle/graphlet data
Generate_normalised_count_csv <- function(cancer_type,
                                          weight_type,
                                          triangle_types,
                                          average = 15,
                                          noCells = 10000,
                                          folder = "Input_data_RaCInG",
                                          remove_direction = TRUE) {

  input_data <- generateInput(weight_type, cancer_type, folder = folder)
  
  CellLig <- input_data$CellLigList
  CellRec <- input_data$CellRecList
  Dtypes <- input_data$Dtypes
  Dconn <- input_data$DconnectionTensor
  celltypes <- input_data$celltypes
  lig <- input_data$ligands
  rec <- input_data$receptors
  signs <- input_data$Sign_matrix
  patient_names <- rownames(input_data$Dtypes)
  
  # Standardize cell type naming
  celltypes[celltypes == "CD8+ T"] <- "CD8"
  
  if (remove_direction) {
    # Initialize data frames for accumulating triangle types
    for (index in seq_along(triangle_types)) {
      triangle_type <- triangle_types[index]
      
      # Read non-normalized and normalized triangle data
      data_raw <- Triangle_Prop_Read(sprintf("%s_%s_%s.out", cancer_type, triangle_type, average))
      data_norm <- Triangle_Prop_Read(sprintf("%s_%s_%s_norm.out", cancer_type, triangle_type, average))
      
      av <- data_raw$triangle_data_a
      avN <- data_norm$triangle_data_a
      summary <- data_raw$triangle_raw_count
      summaryN <- data_norm$triangle_raw_count
      
      # ------------------------------------------------------------
      # Create data frames with columns as "cell1_cell2_cell3" labels
      # ------------------------------------------------------------
      non_unif_data <- list()
      unif_data <- list()
      
      for (i in 1:9) {
        for (j in 1:9) {
          for (k in 1:9) {
            label <- paste(celltypes[i], celltypes[j], celltypes[k], sep = "_")
            # Extract patient data (last dimension)
            non_unif_data[[label]] <- av[i, j, k, ]
            unif_data[[label]] <- avN[i, j, k, ]
          }
        }
      }
      
      df <- as.data.frame(non_unif_data)
      dfN <- as.data.frame(unif_data)
      
      # Initialize accumulation frames
      if (index == 1) {
        new <- data.frame(matrix(nrow = nrow(df), ncol = 0))
        newN <- data.frame(matrix(nrow = nrow(dfN), ncol = 0))
      }
      
      # ------------------------------------------------------------
      # Combine columns based on direction removal
      # ------------------------------------------------------------
      if (triangle_type == "W") {
        # Only sort first and last element of label
        for (colname in colnames(df)) {
          parts <- strsplit(colname, "_")[[1]]
          if (length(parts) > 1) {
            sorted_parts <- sort(c(parts[1], parts[3]))
            parts[1] <- sorted_parts[1]
            parts[3] <- sorted_parts[2]
            new_label <- paste(parts, collapse = "_")
          } else {
            new_label <- colname
          }
          # Accumulate or initialize
          if (new_label %in% colnames(new)) {
            new[[new_label]] <- new[[new_label]] + df[[colname]]
            newN[[new_label]] <- newN[[new_label]] + dfN[[colname]]
          } else {
            new[[new_label]] <- df[[colname]]
            newN[[new_label]] <- dfN[[colname]]
          }
        }
      } else {
        # Fully sort for other triangle types
        for (colname in colnames(df)) {
          sorted_label <- paste(sort(strsplit(colname, "_")[[1]]), collapse = "_")
          if (sorted_label %in% colnames(new)) {
            new[[sorted_label]] <- new[[sorted_label]] + df[[colname]]
            newN[[sorted_label]] <- newN[[sorted_label]] + dfN[[colname]]
          } else {
            new[[sorted_label]] <- df[[colname]]
            newN[[sorted_label]] <- dfN[[colname]]
          }
        }
      }
    }
    
    if (length(triangle_types) > 1) triangle_types <- "Tr"
    
    # ------------------------------------------------------------
    # Normalizing: divide non-uniform by uniform data
    # ------------------------------------------------------------
    av <- as.matrix(new)
    avN <- as.matrix(newN)
    
    # Set entries that are zero in both to 1
    zero_idx <- which(av == 0 & avN == 0, arr.ind = TRUE)
    if (nrow(zero_idx) > 0) {
      for (r in 1:nrow(zero_idx)) {
        av[zero_idx[r, "row"], zero_idx[r, "col"]] <- 1
        avN[zero_idx[r, "row"], zero_idx[r, "col"]] <- 1
      }
    }
    
    # Ensure no zeros in avN
    avN[avN == 0] <- 1
    
    # Normalize
    Norm <- av / avN
    
    # Remove patients where summary failed
    delete_indices <- unique(c(which(summary[,1] == 0), which(summaryN[,1] == 0)))
    if (length(delete_indices) > 0) {
      Norm <- Norm[-delete_indices, , drop = FALSE]
      patients <- patient_names[-delete_indices]
    } else {
      patients <- patient_names
    }
    
    # Convert normalized data to data frame
    df <- as.data.frame(Norm)
    rownames(df) <- patients
    
    # Save CSV
    write.csv(df, sprintf("%s_%s_%s_cells_%s_deg_data_bundle.csv",
                          cancer_type, triangle_types, noCells, average),
              row.names = TRUE)
    return(df)
  }
  
  # ------------------------------------------------------------
  # Case without remove_direction
  # (similar but simpler normalization)
  # ------------------------------------------------------------
  data_raw <- Triangle_Prop_Read(sprintf("%s_%s_%s.out", cancer_type, triangle_types, average))
  data_norm <- Triangle_Prop_Read(sprintf("%s_%s_%s_norm.out", cancer_type, triangle_types, average))
  
  av <- data_raw$triangle_data_a
  avN <- data_norm$triangle_data_a
  summary <- data_raw$triangle_raw_count
  summaryN <- data_norm$triangle_raw_count
  
  # Set zero intersections to 1
  zero_idx <- which(av == 0 & avN == 0, arr.ind = TRUE)
  if (nrow(zero_idx) > 0) {
    for (r in 1:nrow(zero_idx)) {
      av[zero_idx[r, "row"], zero_idx[r, "col"], zero_idx[r, "dim3"], zero_idx[r, "dim4"]] <- 1
      avN[zero_idx[r, "row"], zero_idx[r, "col"], zero_idx[r, "dim3"], zero_idx[r, "dim4"]] <- 1
    }
  }
  
  avN[avN == 0] <- 1
  Norm <- av / avN
  
  delete_indices <- unique(c(which(summary[,1] == 0), which(summaryN[,1] == 0)))
  if (length(delete_indices) > 0) {
    Norm <- Norm[, , , -delete_indices, drop = FALSE]
    patients <- patient_names[-delete_indices]
  } else {
    patients <- patient_names
  }
  
  # Convert normalized data to data frame
  Normdata <- list()
  for (i in 1:9) {
    for (j in 1:9) {
      for (k in 1:9) {
        label <- paste(celltypes[i], celltypes[j], celltypes[k], sep = "_")
        Normdata[[label]] <- Norm[i, j, k, ]
      }
    }
  }
  
  df <- as.data.frame(Normdata)
  rownames(df) <- patients
  
  # Save CSV
  write.csv(df, sprintf("%s_%s_%s_cells_%s_deg_data.csv",
                        cancer_type, triangle_types, noCells, average),
            row.names = TRUE)
  
  return(df)
}