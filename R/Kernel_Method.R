# ==============================================================
# Kernel.R
# Contains functions to compute the kernel and graph properties
# ==============================================================

library(Matrix)
library(pracma)  # for root finding similar to scipy.optimize.root

# --------------------------------------------------------------
# Calculate_kernel: compute the kernel of the random graph model
# --------------------------------------------------------------
Calculate_kernel <- function(liglist, reclist, Dcell, Dconn, normalize = FALSE) {
  # Ensure Dcell is a matrix: if it's a vector (single patient), convert it to 1-row matrix
  Dcell <- if (is.null(dim(Dcell))) matrix(Dcell, nrow = 1) else Dcell

  dim1 <- ncol(Dcell) # number of cell types (columns of Dcell)
  dim2 <- nrow(Dcell) # number of patients (rows of Dcell)

  # Ensure Dconn is always a 3D array: [cell types x cell types x patients]
  Dconn <- array(Dconn, dim = c(nrow(Dconn), ncol(Dconn), dim2))

  # Initialize commexpec: 3D array to store expected communications
  # Dimensions: [ligand cell type x receptor cell type x patient]
  commexpec <- array(0, dim = c(dim1, dim1, dim2))

  # Normalize connection strengths if requested
  if (normalize) {
    # Count the number of positive connections for each patient (3rd dimension)
    # Then invert the counts to get the weight for each positive connection
    normvec <- 1 / apply(Dconn > 0, 3, sum)

    # Loop over each patient
    for (i in 1:dim2) {
      # Extract the 2D connection matrix for this patient
      copy <- Dconn[,,i]

      # Set all positive connections to the normalized value
      # This ensures that the sum of outgoing connections for this patient is 1
      copy[copy > 0] <- normvec[i]

      # Update the original Dconn array
      Dconn[,,i] <- copy
    }
  }

  for (patient in 1:dim2) {
    cat("Patient:", patient, "\n")
    for (i in 1:dim1) {
      for (j in 1:dim1) {
        for (k in 1:dim1) {
          for (l in 1:dim1) {
            weightlig <- sum(Dcell[patient, ] * liglist[, i])
            weightrec <- sum(Dcell[patient, ] * reclist[, j])

            if (weightlig != 0 && weightrec != 0) {
              commexpec[k, l, patient] <- commexpec[k, l, patient] +
                Dconn[i, j, patient] * (liglist[k, i] / weightlig) * (reclist[l, j] / weightrec)
            }
          }
        }
      }
    }
  }

  return(commexpec)
}

# --------------------------------------------------------------
# Save kernel to .RData
# --------------------------------------------------------------
saveKernel <- function(weight, cancer, folder, test = FALSE) {
  input <- generateInput(weight, cancer, folder)
  liglist <- input[[1]]
  reclist <- input[[2]]
  Dcell <- input[[3]]
  Dconn <- input[[4]]

  if (test) {
    Dcell <- Dcell[1:2, ]
    Dconn <- Dconn[,,1:2]
  }

  kernel <- Calculate_kernel(liglist, reclist, Dcell, Dconn, normalize = FALSE)
  unifKernel <- Calculate_kernel(liglist, reclist, Dcell, Dconn, normalize = TRUE)

  save(kernel, unifKernel, file = paste0("kernel_", cancer, ".RData"))
}

# --------------------------------------------------------------
# Calculate direct communication
# --------------------------------------------------------------
calculateDirect <- function(cancer, bundle = TRUE, norm = TRUE, test = FALSE) {
  load(file = paste0("kernel_", cancer, ".RData"))
  input <- generateInput("min", cancer)
  Dcell <- input[[3]]
  cells <- input[[5]]
  cells[cells == "CD8+ T"] <- "CD8"
  names <- get_patient_names(cancer)

  if (test) {
    Dcell <- Dcell[1:2, ]
    names <- names[1:2]
  }

  if (norm) {
    outN <- unifKernel
  }

  df_list <- list()
  for (i in 1:length(cells)) {
    for (j in if (bundle) i:length(cells) else 1:length(cells)) {
      colname <- paste0("Dir_", cells[i], "_", cells[j])
      if (norm) {
        value <- if (bundle) {
          (kernel[i, j, ] + kernel[j, i, ]) / (outN[i, j, ] + outN[j, i, ])
        } else {
          kernel[i, j, ] / outN[i, j, ]
        }
      } else {
        value <- if (bundle) {
          Dcell[, i] * Dcell[, j] * (kernel[i, j, ] + kernel[j, i, ])
        } else {
          Dcell[, i] * Dcell[, j] * kernel[i, j, ]
        }
      }
      df_list[[colname]] <- value
    }
  }

  df <- as.data.frame(df_list)
  rownames(df) <- names
  write.csv(df, file = paste0(cancer, "_min_weight_Dir", if (bundle) "_bundle", if (norm) "_norm", ".csv"))
  return(df)
}

# --------------------------------------------------------------
# Calculate wedge values (W)
# --------------------------------------------------------------
calculateWedges <- function(cancer, bundle = TRUE, norm = TRUE, test = FALSE) {
  load(file = paste0("kernel_", cancer, ".RData"))
  input <- generateInput("min", cancer)
  Dcell <- input[[3]]
  cells <- input[[5]]
  cells[cells == "CD8+ T"] <- "CD8"
  names <- get_patient_names(cancer)

  if (test) {
    Dcell <- Dcell[1:2, ]
    names <- names[1:2]
  }

  df_list <- list()
  for (i in 1:length(cells)) {
    for (j in 1:length(cells)) {
      for (k in if (bundle) j:length(cells) else 1:length(cells)) {
        colname <- paste0("W_", cells[i], "_", cells[j], "_", cells[k])
        if (norm) {
          value <- if (bundle) {
            (kernel[i,j,] * kernel[j,k,] + kernel[k,j,]*kernel[j,i,]) /
              (unifKernel[i,j,] * unifKernel[j,k,] + unifKernel[k,j,]*unifKernel[j,i,])
          } else {
            (kernel[i,j,] * kernel[j,k,]) / (unifKernel[i,j,] * unifKernel[j,k,])
          }
        } else {
          value <- if (bundle) {
            Dcell[,i] * Dcell[,j] * Dcell[,k] * (kernel[i,j,] * kernel[j,k,] + kernel[k,j,]*kernel[j,i,])
          } else {
            Dcell[,i] * Dcell[,j] * Dcell[,k] * (kernel[i,j,] * kernel[j,k,])
          }
        }
        df_list[[colname]] <- value
      }
    }
  }

  df <- as.data.frame(df_list)
  rownames(df) <- names
  write.csv(df, file = paste0(cancer, "_min_weight_W", if (bundle) "_bundle", if (norm) "_norm", ".csv"))
  return(df)
}

# --------------------------------------------------------------
# Calculate GSCC analytically
# --------------------------------------------------------------
getGSCCAnalytically <- function(cancer, lab = 1, norm = TRUE, test = FALSE) {
  load(file = paste0("kernel_", cancer, ".RData"))
  input <- generateInput("min", cancer)
  Dcell <- input[[3]]
  cells <- input[[5]]
  cells[cells == "CD8+ T"] <- "CD8"
  names <- get_patient_names(cancer)

  if (test) {
    Dcell <- Dcell[1:2, ]
    names <- names[1:2]
  }

  sens <- length(cells)
  GSCCsizes <- matrix(0, nrow = length(names), ncol = sens+1)
  GSCCsizesN <- matrix(0, nrow = length(names), ncol = sens+1)

  for (k in 1:length(names)) {
    q <- Dcell[k, ]
    muP <- matrix(0, nrow = sens, ncol = sens)
    muM <- matrix(0, nrow = sens, ncol = sens)
    for (i in 1:sens) {
      for (j in 1:sens) {
        muP[i,j] <- lab * kernel[j,i,k] * q[j]
        muM[i,j] <- lab * kernel[i,j,k] * q[j]
      }
    }
    solP <- nleqslv::nleqslv(rep(1, sens), poiBPFunc, muP = muP, sens = sens)
    solM <- nleqslv::nleqslv(rep(1, sens), poiBPFunc, muP = muM, sens = sens)
    x <- solP$x
    y <- solM$x
    GSCCsizes[k, 1:sens] <- x * y * q
    GSCCsizes[k, sens+1] <- sum(x * y * q)
  }

  # Normalize if required
  if (norm) {
    normvals <- GSCCsizes / GSCCsizesN
    normvals[is.nan(normvals)] <- 1
    df <- as.data.frame(normvals)
    colnames(df) <- c(paste0("GSCC_", cells), "GSCC")
    rownames(df) <- names
    write.csv(df, file = paste0(cancer, "_min_weight_GSCC_norm.csv"))
  } else {
    df <- as.data.frame(GSCCsizes)
    colnames(df) <- c(paste0("GSCC_", cells), "GSCC")
    rownames(df) <- names
    write.csv(df, file = paste0(cancer, "_min_weight_GSCC.csv"))
  }

  return(df)
}
