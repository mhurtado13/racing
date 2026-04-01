#' Compute the RaCInG kernel for one or more patients
#'
#' @param liglist Cell-by-ligand compatibility matrix.
#' @param reclist Cell-by-receptor compatibility matrix.
#' @param Cmatrix Patient-by-cell-type abundance matrix.
#' @param LRmatrix Ligand-by-receptor-by-patient interaction tensor.
#' @param normalize Logical; if `TRUE`, also compute a uniformized baseline kernel.
#'
#' @return Either a 3D kernel array or a list with `kernel` and `kernel_norm`.
#' @export
compute_kernel <- function(liglist, reclist, Cmatrix, LRmatrix, normalize = FALSE) {
  # liglist:   [num_genes × num_ligands] binary/weighted matrix mapping genes to ligands
  # reclist:   [num_genes × num_receptors] binary/weighted matrix mapping genes to receptors
  # Cmatrix:     [num_patients × num_genes] cell type abundances for each patient
  # LRmatrix:     [num_ligands × num_receptors × num_patients] ligand–receptor interaction strengths for each patient
  # normalize: if TRUE, normalize LRmatrix so that each patient has the same total interaction strength
  #
  # Return:
  #  - if normalize == FALSE: kernel array [sender x receiver x patient]
  #  - if normalize == TRUE: list(kernel = unnormalized_kernel, kernel_norm = normalized_kernel)

  # Ensure matrix format
  Cmatrix <- if (is.null(dim(Cmatrix))) matrix(Cmatrix, nrow = 1) else Cmatrix
  n_patients <- nrow(Cmatrix)
  n_celltypes <- ncol(Cmatrix)

  # Prepare normalized LRmatrix copy if requested (do not overwrite original LRmatrix)
  if (normalize) {
    LRmatrix_norm <- LRmatrix
    normvec <- 1 / apply(LRmatrix_norm > 0, 3, sum)  # inverse of positive interactions per patient
    for (p in 1:n_patients) {
      copy <- LRmatrix_norm[, , p]
      copy[copy > 0] <- normvec[p]           # scale positive interactions
      LRmatrix_norm[, , p] <- copy
    }
  }

  # Output: [sender x receiver x patient]
  kernel <- array(0, dim = c(n_celltypes, n_celltypes, n_patients))
  if (normalize) kernel_norm <- array(0, dim = c(n_celltypes, n_celltypes, n_patients))

  for (p in 1:n_patients) {

    Cvec <- Cmatrix[p, ]  # cell abundances

    ############################# LIGANDS #############################

    # ---- Step 1: ligand weights ----
    # Row-wise multiplication:
    # each row (cell type A) is scaled by its abundance C(A) in patient p,
    # producing abundance-weighted ligand expression C(A) * L(A,i) per patient
    lig_weight_raw <- sweep(liglist, 1, Cvec, "*")

    # ---- Step 1.5: Calculate normalization factors for ligands ----
    # For each ligand i, sum contributions from all cell types:
    # gives total C(A) * L(A,i) used for normalization
    lig_norm <- colSums(lig_weight_raw) # How much total of this ligand exists across all cell types, weighted by abundance?

    # avoid division by zero
    lig_norm[lig_norm == 0] <- 1

    # ---- Step 2: normalize ligand weights ----
    # dim: [celltype x ligand]
    lig_weight <- sweep(lig_weight_raw, 2, lig_norm, "/") #turns raw values into a probability distribution over sender cell types

    ############################## RECEPTORS #############################

    # ---- Step 2: receptor weights ----
    rec_weight_raw <- sweep(reclist, 1, Cvec, "*") ## Abundance-weighted receptor expression C(A) * R(A,j) per patient
    rec_norm <- colSums(rec_weight_raw) ## Normalization factor for each receptor: total C(A) * R(A,j) across all cell types
    rec_norm[rec_norm == 0] <- 1 ## avoid division by zero

    rec_weight <- sweep(rec_weight_raw, 2, rec_norm, "/") # normalized receptor weights: probability distribution over receiver cell types for each receptor

    # ---- Step 3: aggregate over ligand–receptor pairs ----
    # Matrix multiplication performs the full kernel computation:
    #
    # lig_weight:   [cell A × ligand i]
    #   → contribution of each sender cell type A to ligand i
    #
    # LRmatrix:        [ligand i × receptor j]
    #   → interaction strength between ligand i and receptor j
    #
    # t(rec_weight): [receptor j × cell B]
    #   → contribution of receptor j to receiver cell type B
    #
    # Multiplication order:
    #   (lig_weight %*% LRmatrix)        → [cell A × receptor j]
    #   (... %*% t(rec_weight))       → [cell A × cell B]
    #
    # Final result:
    #   kernel[A, B] = sum over all ligand–receptor pairs (i,j) of
    #                  (sender contribution) × (LR interaction) × (receiver contribution)
    #
    # Interpretation:
    #   Expected interaction strength from cell type A → cell type B
    #   aggregated over all ligand–receptor pairs

    kernel[, , p] <- lig_weight %*% LRmatrix[, , p] %*% t(rec_weight)

    # If requested, compute normalized kernel for this patient using LRmatrix_norm
    if (normalize) {
      kernel_norm[, , p] <- lig_weight %*% LRmatrix_norm[, , p] %*% t(rec_weight)
    }
  }

  if (normalize) {
    return(list(kernel = kernel, kernel_norm = kernel_norm))
  } else {
    return(kernel)
  }
}

#' Calculate direct communication features from a kernel
#'
#' @param kernel Kernel array returned by [compute_kernel()].
#' @param unifKernel Optional normalized baseline kernel.
#' @param cells Character vector of cell-type names.
#' @param bundle Logical; if `TRUE`, combine reciprocal directions.
#'
#' @return A patient-by-feature data frame of direct communication scores.
#' @export
calculateDirect <- function(kernel, unifKernel = NULL, cells, bundle = TRUE) {
  #---------------------------------------------------------------
  # Function: calculateDirectSimple
  # Purpose:  Compute direct communication scores between cell types
  #           from precomputed kernel arrays.
  #
  # Inputs:
  #   kernel      - 3D array [sender cell x receiver cell x patient], raw scores
  #   unifKernel  - optional 3D array, same dimensions as kernel, normalized scores
  #   cells       - character vector of cell type names (length = dim(kernel)[1])
  #   bundle      - logical; if TRUE, sum reciprocal interactions (A->B + B->A)
  #
  # Output:
  #   Data frame: rows = patients, columns = direct interaction scores for each cell type pair
  #---------------------------------------------------------------

  n_celltypes <- length(cells)      # number of cell types
  n_patients <- dim(kernel)[3]      # number of patients

  df_list <- list()                 # list to store each column before converting to data frame

  #---------------------------------------------------------------
  # Loop over all cell type pairs (sender i, receiver j)
  # If bundle = TRUE, we only consider upper triangle to avoid double counting
  #---------------------------------------------------------------
  for (i in 1:n_celltypes) {
    for (j in if (bundle) i:n_celltypes else 1:n_celltypes) {

      colname <- paste0("Dir_", cells[i], "_", cells[j])  # column name for this pair

      #-----------------------------------------------------------
      # Compute direct communication value
      # - If unifKernel is provided, normalize by the corresponding uniform kernel
      # - If bundle = TRUE, sum the reciprocal interaction (A->B + B->A)
      #-----------------------------------------------------------
      if (!is.null(unifKernel)) {
        # normalized direct score
        value <- if (bundle) {
          (kernel[i, j, ] + kernel[j, i, ]) / (unifKernel[i, j, ] + unifKernel[j, i, ])
        } else {
          kernel[i, j, ] / unifKernel[i, j, ]
        }
      } else {
        # raw direct score (already abundance-weighted in kernel)
        value <- if (bundle) {
          kernel[i, j, ] + kernel[j, i, ]
        } else {
          kernel[i, j, ]
        }
      }

      df_list[[colname]] <- value   # add the vector of patient scores as a new column
    }
  }

  #---------------------------------------------------------------
  # Convert list of columns into a data frame
  # Rows are patients, columns are cell type pair direct interactions
  #---------------------------------------------------------------
  df <- as.data.frame(df_list)
  rownames(df) <- paste0("Patient_", 1:n_patients)

  return(df)
}

#' Calculate wedge features from a kernel
#'
#' @param kernel Kernel array returned by [compute_kernel()].
#' @param unifKernel Optional normalized baseline kernel.
#' @param cells Character vector of cell-type names.
#' @param bundle Logical; if `TRUE`, combine directionally equivalent wedges.
#'
#' @return A patient-by-feature data frame of wedge scores.
#' @export
calculateWedges <- function(kernel, unifKernel = NULL, cells, bundle = TRUE) {
  #---------------------------------------------------------------
  # Function: calculateWedgesSimple
  # Purpose:  Compute "wedge" values W for triplets of cell types
  #           from precomputed kernel arrays.
  #
  # Inputs:
  #   kernel      - 3D array [cell A x cell B x patient], direct interaction scores
  #   unifKernel  - optional 3D array, same dimensions as kernel, normalized scores
  #   cells       - character vector of cell type names (length = dim(kernel)[1])
  #   bundle      - logical; if TRUE, sum reciprocal interactions in wedge calculation
  #
  # Output:
  #   Data frame: rows = patients, columns = wedge scores for each cell triplet
  #---------------------------------------------------------------

  n_celltypes <- length(cells)
  n_patients <- dim(kernel)[3]

  df_list <- list()   # list to collect columns before converting to data frame

  #---------------------------------------------------------------
  # Loop over all triplets of cell types (i, j, k)
  # i = sender, j = intermediate, k = receiver
  #---------------------------------------------------------------
  for (i in 1:n_celltypes) {
    for (j in 1:n_celltypes) {
      for (k in if (bundle) j:n_celltypes else 1:n_celltypes) {

        colname <- paste0("W_", cells[i], "_", cells[j], "_", cells[k])

        #-----------------------------------------------------------
        # Compute wedge score
        # - If unifKernel provided, normalize by corresponding uniform kernel
        # - If bundle = TRUE, sum reciprocal contributions (i<->k via j)
        #-----------------------------------------------------------
        if (!is.null(unifKernel)) {
          value <- if (bundle) {
            (kernel[i,j,] * kernel[j,k,] + kernel[k,j,] * kernel[j,i,]) /
              (unifKernel[i,j,] * unifKernel[j,k,] + unifKernel[k,j,] * unifKernel[j,i,])
          } else {
            (kernel[i,j,] * kernel[j,k,]) /
              (unifKernel[i,j,] * unifKernel[j,k,])
          }
        } else {
          # raw wedge (kernel already includes cell abundances)
          value <- if (bundle) {
            kernel[i,j,] * kernel[j,k,] + kernel[k,j,] * kernel[j,i,]
          } else {
            kernel[i,j,] * kernel[j,k,]
          }
        }

        df_list[[colname]] <- value   # store column for this triplet
      }
    }
  }

  #---------------------------------------------------------------
  # Convert list of columns into a data frame
  # Rows = patients, Columns = wedges for each triplet
  #---------------------------------------------------------------
  df <- as.data.frame(df_list)
  rownames(df) <- paste0("Patient_", 1:n_patients)

  return(df)
}

#' Legacy GSCC helper
#'
#' This helper depended on project-specific files from the original development
#' workflow. The packaged interface now recommends using [computeGSCC()] directly.
#'
#' @param cancer Legacy dataset identifier.
#' @param lab Scaling factor.
#' @param norm Logical; kept for backward compatibility.
#' @param test Logical; kept for backward compatibility.
#'
#' @return This function stops with a message directing users to [computeGSCC()].
#' @export
getGSCCAnalytically <- function(cancer, lab = 1, norm = TRUE, test = FALSE) {
  stop(
    "`getGSCCAnalytically()` is a legacy helper. Use `computeGSCC()` with an explicit kernel and abundance matrix instead.",
    call. = FALSE
  )
}

#' Compute GSCC features from kernel matrices
#'
#' @param kernel Kernel array from [compute_kernel()].
#' @param Dcell Patient-by-cell-type abundance matrix.
#' @param cell_names Character vector of cell-type labels.
#' @param patient_names Character vector of patient names.
#' @param unifKernel Optional normalized baseline kernel.
#' @param norm Logical; if `TRUE`, divide by the baseline GSCC values.
#' @param lab Scaling factor for interaction strengths.
#'
#' @return A data frame with GSCC feature values per patient.
#' @export
computeGSCC <- function(kernel, Dcell, cell_names, patient_names,
                        unifKernel = NULL, norm = FALSE, lab = 1) {
  # kernel: [sender x receiver x patient] matrix from compute_kernel
  # Dcell: [patient x cell type] cell abundances (can be used to scale if needed)
  # cell_names: vector of cell type names
  # patient_names: vector of patient names
  # unifKernel: optional uniform kernel (same dimensions as kernel) for normalization
  # norm: boolean, whether to normalize by uniform kernel
  # lab: scaling factor for interaction strengths

  n_patients <- nrow(Dcell)
  n_cells <- length(cell_names)

  # Matrices to store GSCC per patient and per cell type
  GSCCsizes <- matrix(0, nrow = n_patients, ncol = n_cells + 1)  # last column = total GSCC
  GSCCsizesN <- if (!is.null(unifKernel)) matrix(0, nrow = n_patients, ncol = n_cells + 1) else NULL

  # Loop over patients
  for (p in 1:n_patients) {
    q <- Dcell[p, ]  # cell abundances vector for this patient

    # Construct mu matrices for Poisson branching process
    muP <- matrix(0, nrow = n_cells, ncol = n_cells)
    muM <- matrix(0, nrow = n_cells, ncol = n_cells)
    for (i in 1:n_cells) {
      for (j in 1:n_cells) {
        muP[i, j] <- lab * kernel[j, i, p] * q[j]  # incoming contribution
        muM[i, j] <- lab * kernel[i, j, p] * q[j]  # outgoing contribution
      }
    }

    # Solve the branching process equations using nleqslv
    solP <- nleqslv::nleqslv(rep(1, n_cells), poiBPFunc, M = muP, sens = n_cells)
    solM <- nleqslv::nleqslv(rep(1, n_cells), poiBPFunc, M = muM, sens = n_cells)

    x <- solP$x
    y <- solM$x

    # GSCC per cell type and total
    GSCCsizes[p, 1:n_cells] <- x * y * q
    GSCCsizes[p, n_cells + 1] <- sum(x * y * q)

    # If uniform kernel is provided, calculate same for normalization
    if (!is.null(unifKernel)) {
      muP_N <- matrix(0, nrow = n_cells, ncol = n_cells)
      muM_N <- matrix(0, nrow = n_cells, ncol = n_cells)
      for (i in 1:n_cells) {
        for (j in 1:n_cells) {
          muP_N[i, j] <- lab * unifKernel[j, i, p] * q[j]
          muM_N[i, j] <- lab * unifKernel[i, j, p] * q[j]
        }
      }
      solP_N <- nleqslv::nleqslv(rep(1, n_cells), poiBPFunc, M = muP_N, sens = n_cells)
      solM_N <- nleqslv::nleqslv(rep(1, n_cells), poiBPFunc, M = muM_N, sens = n_cells)
      xN <- solP_N$x
      yN <- solM_N$x
      GSCCsizesN[p, 1:n_cells] <- xN * yN * q
      GSCCsizesN[p, n_cells + 1] <- sum(xN * yN * q)
    }
  }

  # Normalize by uniform kernel if requested
  if (norm && !is.null(unifKernel)) {
    normvals <- GSCCsizes / GSCCsizesN
    normvals[is.nan(normvals)] <- 1  # handle division by zero
    df <- as.data.frame(normvals)
  } else {
    df <- as.data.frame(GSCCsizes)
  }

  # Set column and row names
  colnames(df) <- c(paste0("GSCC_", cell_names), "GSCC_total")
  rownames(df) <- patient_names

  return(df)
}

#' Compute triangle features from kernel matrices
#'
#' @param kernel Kernel array from [compute_kernel()].
#' @param cell_names Character vector of cell-type names.
#' @param patient_names Character vector of patient names.
#' @param unifKernel Optional normalized baseline kernel.
#' @param norm Logical; if `TRUE`, divide by the baseline triangle scores.
#' @param bundle Logical; if `TRUE`, aggregate directionally equivalent triangles.
#'
#' @return A patient-by-feature data frame of triangle scores.
#' @export
computeTriangles <- function(kernel, cell_names, patient_names,
                             unifKernel = NULL, norm = FALSE, bundle = TRUE) {
  # kernel: 3D array [sender x receiver x patient] from compute_kernel
  #         Already includes cell abundances, so no extra weighting needed
  # cell_names: vector of cell type names
  # patient_names: vector of patient IDs/names
  # unifKernel: optional 3D array same shape as kernel, for normalization
  # norm: if TRUE, divide triangle scores by the corresponding uniform kernel
  # bundle: if TRUE, sum all permutations of three nodes to remove directionality

  n_patients <- dim(kernel)[3]   # number of patients
  n_cells <- length(cell_names)   # number of cell types

  df_list <- list()  # will store triangle values per patient

  # Loop over all triples of cell types
  for (i in 1:n_cells) {
    for (j in if (bundle) i:n_cells else 1:n_cells) {  # avoid double counting if bundling
      for (k in if (bundle) j:n_cells else 1:n_cells) {

        # Create a column name that describes the triangle
        colname <- if (bundle) {
          paste0("Tr_", cell_names[i], "_", cell_names[j], "_", cell_names[k])
        } else {
          paste0("TT_", cell_names[i], "_", cell_names[j], "_", cell_names[k])
        }

        # Initialize vector to store triangle value per patient
        num <- rep(0, n_patients)

        # Determine which permutations of the triple to consider
        if (bundle) {
          # Include all 8 permutations to remove directionality
          perms <- list(
            c(i,j,k), c(i,j,k), c(i,k,j), c(i,k,j),
            c(j,i,k), c(j,k,i), c(j,k,i), c(k,j,i)
          )
        } else {
          # Only one permutation for directed triangles (trust triangles)
          perms <- list(c(i,j,k))
        }

        # Compute the triangle score for each patient
        # For each permutation, multiply the kernel values along the triangle path
        for (p in perms) {
          # kernel[p[1],p[2],] = contribution from cell type p[1] to p[2]
          # multiply the three edges to get triangle contribution
          num <- num + kernel[p[1], p[2], ] * kernel[p[2], p[3], ] * kernel[p[3], p[1], ]
        }

        # Normalize by uniform kernel if requested
        if (norm && !is.null(unifKernel)) {
          denom <- rep(0, n_patients)  # store denominator (triangle in uniform kernel)
          for (p in perms) {
            denom <- denom + unifKernel[p[1], p[2], ] * unifKernel[p[2], p[3], ] * unifKernel[p[3], p[1], ]
          }
          num <- num / denom   # divide numerator by denominator
          num[is.nan(num)] <- 1  # handle divisions by zero
        }

        # Store the triangle values for this triple
        df_list[[colname]] <- num
      }
    }
  }

  # Convert the list into a dataframe: rows = patients, columns = triangles
  df <- as.data.frame(df_list)
  rownames(df) <- patient_names

  return(df)
}

#' Derive communication features from a kernel
#'
#' @param kernel Kernel array from [compute_kernel()].
#' @param unifKernel Optional normalized baseline kernel.
#' @param celltypes Character vector of cell-type labels.
#' @param communication_type Feature family to compute (`"D"`, `"W"`, `"TT"`, or `"GSCC"`).
#' @param bundle Logical; if `TRUE`, merge directionally equivalent features where appropriate.
#' @param patient_names Optional patient labels.
#' @param Dcell Optional abundance matrix required for `"GSCC"`.
#' @param norm Logical; if `TRUE`, compute normalized features when a baseline is supplied.
#' @param patient_idx Optional patient index subset.
#'
#' @return A data frame of feature values for the selected patients.
#' @export
compute_kernel_features <- function(kernel, unifKernel = NULL, celltypes, communication_type = "D", bundle = TRUE,
                                    patient_names = NULL, Dcell = NULL, norm = FALSE, patient_idx = NULL) {

  # Subset patients if requested
  if (!is.null(patient_idx)) {
    kernel <- kernel[,,patient_idx, drop = FALSE]
    if (!is.null(unifKernel)) unifKernel <- unifKernel[,,patient_idx, drop = FALSE]
    if (!is.null(patient_names)) patient_names <- patient_names[patient_idx]
    if (!is.null(Dcell)) Dcell <- Dcell[patient_idx, , drop = FALSE]
  }

  comm <- toupper(communication_type)

  if (comm == "D") {
    df <- calculateDirect(kernel = kernel, unifKernel = unifKernel, cells = celltypes, bundle = bundle)

  } else if (comm == "W") {
    df <- calculateWedges(kernel = kernel, unifKernel = unifKernel, cells = celltypes, bundle = bundle)

  } else if (comm == "TT") {
    # trust triangles (directed): use bundle = FALSE in computeTriangles
    df <- computeTriangles(kernel = kernel, cell_names = celltypes,
                           patient_names = if (!is.null(patient_names)) patient_names else paste0("Patient_", seq_len(dim(kernel)[3])),
                           unifKernel = unifKernel, norm = norm, bundle = FALSE)

  } else if (comm == "GSCC") {
    if (is.null(Dcell)) stop("Dcell (cell abundance matrix) is required for GSCC computation.")
    # computeGSCC expects patient_names and Dcell rows == patients
    df <- computeGSCC(kernel = kernel, Dcell = Dcell,
                      cell_names = celltypes,
                      patient_names = if (!is.null(patient_names)) patient_names else paste0("Patient_", seq_len(nrow(Dcell))),
                      unifKernel = unifKernel, norm = norm)

  } else {
    stop("Unsupported communication_type. Use D, W, TT, or GSCC.")
  }

  return(df)

}

#' Run the full kernel-based RaCInG workflow
#'
#' @param counts Gene-by-sample count matrix. Required when `input_data` is not
#'   supplied; ignored otherwise.
#' @param output_folder Directory used to write and read intermediate input files.
#' @param deconv Optional deconvolution matrix.
#' @param cc_network Optional ligand-receptor prior network.
#' @param fun_LR Function used to combine ligand and receptor expression values.
#' @param cell_expr_profile Optional cell-type expression profile matrix.
#' @param source,target Column names to use as ligand and receptor identifiers in `cc_network`.
#' @param signed Logical; if `TRUE`, also try to load a sign matrix.
#' @param deconv_method Deconvolution method used when `deconv` is not supplied.
#' @param cbsx.name,cbsx.token Optional credentials for the deconvolution workflow.
#' @param file_name File stem used for intermediate input files.
#' @param nPatients Number of patients to process, or `"all"`.
#' @param communication_type Feature family to compute.
#' @param norm Logical; if `TRUE`, compute a normalized baseline kernel.
#' @param pt_idx Optional single patient index to process.
#' @param remove_direction Logical; if `TRUE`, merge directionally equivalent features.
#' @param input_data Optional named list of pre-computed input matrices as returned
#'   by [prepare_input_files()].
#'   Must contain `Lmatrix`, `Rmatrix`, `Cmatrix`, `LRmatrix`, `celltypes`,
#'   `ligands`, and `receptors`.
#'   When supplied, the `counts` argument and all preprocessing parameters
#'   (`deconv`, `cc_network`, etc.) are ignored.
#'
#' @return A list with the kernel arrays and the derived feature matrix.
#' @export
compute_racing_kernel = function(counts = NULL, output_folder = "~/Documents/racing/vignettes/", deconv = NULL, cc_network = NULL, fun_LR = min, 
                                 cell_expr_profile = NULL, source = "source_genesymbol", target = "target_genesymbol", signed = FALSE,
                                 deconv_method = "Quantiseq", cbsx.name = NULL, cbsx.token = NULL, file_name = NULL, nPatients = "all", 
                                 communication_type = "W", norm = TRUE, pt_idx = NULL, remove_direction = TRUE,
                                 input_data = NULL) {

  if (is.null(file_name)) {
    file_name <- "RaCInG_input"
  }

  if (!is.null(input_data)) {
    cat("Using pre-computed input matrices; skipping input generation.\n")
    res <- input_data
  } else {
    if (is.null(counts)) {
      stop("Either 'counts' or 'input_data' must be provided.", call. = FALSE)
    }
    res <- prepare_input_files(counts, output_folder = output_folder, deconv = deconv, cc_network = cc_network, fun_LR = fun_LR, 
                               cell_expr_profile = cell_expr_profile, source = source, target = target,
                               deconv_method = deconv_method, cbsx.name = cbsx.name, cbsx.token = cbsx.token, file_name = file_name,
                               signed = signed)
  }

  Lmatrix   <- res$Lmatrix
  Rmatrix   <- res$Rmatrix
  Cmatrix    <- res$Cmatrix
  LRmatrix   <- res$LRmatrix
  cellTypes <- res$celltypes
  ligs      <- res$ligands
  recs      <- res$receptors

  if(nPatients == "all"){
    nPatients = nrow(Cmatrix)
  }                        

  if(!is.null(pt_idx)){
    cat("Because patient index is provided, nPatients argument will be ignored.\n")
  }else{
    cat("Computing kernel for ", ifelse(nPatients == "all", "all", nPatients), " patients\n")
    Cmatrix  <- Cmatrix[1:nPatients, , drop = FALSE]
    LRmatrix <- LRmatrix[,,1:nPatients, drop = FALSE]
  }

  # -----------------------------
  # Calculate kernel
  # -----------------------------
  cat("Calculating kernel...\n")
  
  res <- compute_kernel(Lmatrix, Rmatrix, Cmatrix, LRmatrix, normalize = norm)
  kernel_raw <- if (isTRUE(norm)) res$kernel else res
  kernel_norm <- if (isTRUE(norm)) res$kernel_norm else NULL

  # -----------------------------
  # Calculate features
  # -----------------------------
  cat("Calculating features...\n")
  patient_names <- if (!is.null(counts)) colnames(counts) else rownames(Cmatrix)
  features <- compute_kernel_features(
    kernel = kernel_raw,
    unifKernel = kernel_norm,
    celltypes = cellTypes,
    communication_type = communication_type,
    bundle = remove_direction,
    patient_names = patient_names,
    Dcell = Cmatrix,
    norm = norm,
    patient_idx = pt_idx
  )

  return(list(kernel = kernel_raw, kernel_norm = kernel_norm, features = features))

}