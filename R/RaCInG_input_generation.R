# ================================
# Required libraries
# ================================
library(data.table)
library(dplyr)

# ================================
# 1. sortPermute
# ================================
sortPermute <- function(stringlist) {
  ord <- order(stringlist)
  sortlist <- stringlist[ord]
  perm <- ord
  return(list(sortlist = sortlist, perm = perm))
}

# ================================
# 2. createCellLigList
# ================================
createCellLigList <- function(filename) {
  if (!file.exists(filename)) {
    print(filename)
    stop("ERROR: The file with cell-ligand interactions does not exist...")
  }
  
  df <- read.csv(filename, row.names = 1)
  
  # Store the column names (ligands) and row names (cell types)
  ligands <- colnames(df)
  celltypes <- rownames(df)
  
  # Convert the data frame into a numeric matrix while removing rownames
  CellLigList = df %>%
    tibble::rownames_to_column("Celltype") %>%
    dplyr::select(-Celltype) %>%
    dplyr::mutate(dplyr::across(dplyr::everything(), as.numeric)) %>%
    as.matrix() 
  
  # Sort the cell types alphabetically and get the permutation indices
  sp <- sortPermute(celltypes)
  celltypes <- sp$sortlist # sorted cell type names
  perm <- sp$perm # permutation to reorder rows
  
  CellLigList <- CellLigList[perm, , drop = FALSE]
  
  return(list(CellLigList = CellLigList,
              ligands = ligands,
              celltypes = celltypes))
}

# ================================
# 3. createCellRecList
# ================================
createCellRecList <- function(filename) {
  if (!file.exists(filename)) {
    stop("ERROR: The file with cell-receptor interactions does not exist...")
  }
  
  df <- read.csv(filename, row.names = 1)
  
  # Store the column names (receptors) and row names (cell types)
  receptors <- colnames(df)
  celltypes <- rownames(df)

  # Convert the data frame into a numeric matrix while removing rownames
  CellRecList = df %>%
    tibble::rownames_to_column("Celltype") %>%  # move rownames into a column
    dplyr::select(-Celltype) %>%                # remove the rowname column
    dplyr::mutate(dplyr::across(dplyr::everything(), as.numeric)) %>%  # convert all columns to numeric
    as.matrix()   

  # Sort the cell types alphabetically and get the permutation indices
  sp <- sortPermute(celltypes)
  celltypes <- sp$sortlist
  perm <- sp$perm
  
  CellRecList <- CellRecList[perm, , drop = FALSE]
  
  return(list(CellRecList = CellRecList,
              receptors = receptors,
              celltypes = celltypes))
}

# ================================
# 4. createCellTypeDistr
# ================================
createCellTypeDistr <- function(cells, filename) {
  if (!file.exists(filename)) {
    stop("Cell type distribution file does not exist...")
  }
  
  # Read CSV file, using the first column as rownames (datasets)
  df <- read.csv(filename, row.names = 1)
  
  # Extract column names as cell type names and row names as dataset names
  celltypes <- colnames(df)
  datasets <- rownames(df)
  
  # Convert the data frame to a numeric matrix and remove rownames
  Dtypes = df %>%
    tibble::rownames_to_column("Datasets") %>%
    dplyr::select(-Datasets) %>%
    dplyr::mutate(dplyr::across(dplyr::everything(), as.numeric)) %>%
    as.matrix()
  
  # Normalize rows so each dataset sums to 1 (relative fractions) --> ############### Estimates from deconv should be always proportions (this should not be necessary unless there is MCP or XCell)
  row_sums <- rowSums(Dtypes)
  normalized <- Dtypes / row_sums
  
  # Sort cell types alphabetically
  sp <- sortPermute(celltypes)
  celltypes <- sp$sortlist
  perm <- sp$perm
  Dtypes <- normalized[, perm, drop = FALSE]
  
  # Merge M1 and M2 if needed --> This will break for other deconvolution results (need to be fix)
  if (length(cells) != length(celltypes)) {
    idx_M1 <- which(celltypes == "M1")
    idx_M2 <- which(celltypes == "M2")
    
    if (length(idx_M1) > 0 && length(idx_M2) > 0) {
      Dtypes[, idx_M1] <- Dtypes[, idx_M1] + Dtypes[, idx_M2]
      Dtypes <- Dtypes[, -idx_M2, drop = FALSE]
      celltypes[idx_M1] <- "M"
      celltypes <- celltypes[-idx_M2]
    }
  }
  
  return(list(Dtypes = Dtypes,
              celltypes = celltypes,
              datasets = datasets))
}

# ================================
# 5. createInteractionDistr
# ================================
createInteractionDistr <- function(filename, ligands, receptors) {
  if (!file.exists(filename)) {
    stop("Interaction distribution file does not exist...")
  }
  
  df <- read.csv(filename, row.names = 1)
  
  # Store column names (interaction names like LIG_REC) and row names (datasets)
  interactions <- colnames(df)
  datasets <- rownames(df)
  
  # Convert data frame to numeric matrix, replacing NA values with 0
  distr <- df %>%
    tibble::rownames_to_column("Datasets") %>%
    dplyr::select(-Datasets) %>%
    dplyr::mutate(dplyr::across(
      dplyr::everything(),                      # apply to all columns
      ~ tidyr::replace_na(as.numeric(.), 0)     # convert to numeric and replace NA with 0
    )) %>%
    as.matrix() 
  
  # Split interaction names "LIG_REC" into ligand and receptor parts
  parts <- strsplit(interactions, "_")
  rowWordIndex <- sapply(parts, `[`, 1)
  colWordIndex <- sapply(parts, `[`, 2)
  
  # Convert ligand/receptor names into numeric indices matching the provided ligand/receptor lists
  row <- match(rowWordIndex, ligands)
  col <- match(colWordIndex, receptors)
  
  # Create a 3D tensor to hold ligand-receptor probabilities
  # Dimensions: ligands x receptors x datasets
  DconnectionTensor <- array(0, dim = c(length(ligands), length(receptors), nrow(distr)))
  
  # Fill the tensor: for each dataset, assign probability values to the correct ligand-receptor pairs
  for (i in 1:nrow(distr)) {
    DconnectionTensor[cbind(row, col, i)] <- distr[i, ]
  }
  
  # Normalize each dataset slice so that the sum of all probabilities equals 1
  # This ensures that each dataset represents a valid probability distribution over all ligand-receptor pairs
  normsum <- apply(DconnectionTensor, 3, sum)
  
  for (i in 1:length(normsum)) {
    if (normsum[i] != 0) { # avoids division by 0
      DconnectionTensor[, , i] <- DconnectionTensor[, , i] / normsum[i]
    }
  }
  
  # Return the normalized 3D tensor
  return(DconnectionTensor)
}

# ================================
# 6. Read_Lig_Rec_Interaction
# ================================
Read_Lig_Rec_Interaction <- function(filename) {
  if (!file.exists(filename)) {
    stop("Interaction sign file does not exist...")
  }
  
  df <- fread(filename, header = TRUE)
  
  receptor_names <- colnames(df)[-1]
  ligand_names <- df[[1]]
  
  sign_matrix <- as.matrix(df[, -1])
  storage.mode(sign_matrix) <- "numeric"
  
  return(list(sign_matrix = sign_matrix,
              ligand_names = ligand_names,
              receptor_names = receptor_names))
}

# ================================
# 7. generateInput
# ================================
generateInput <- function(weight_type, cancer_name,
                          read_signs = FALSE,
                          folder = "Example input") {
  
  # -----------------------------
  # Read cell-to-ligand compatibility information
  # The function returns:
  #   CellLigList: numeric matrix indicating which cell type can interact with which ligand
  #   ligands: vector of ligand names
  #   celltypes: vector of cell type names (alphabetically sorted)
  # -----------------------------
  lig_data <- createCellLigList(file.path(folder, "celltype_ligand.csv"))
  
  # -----------------------------
  # Read cell-to-receptor compatibility information
  # The function returns:
  #   CellRecList: numeric matrix indicating which cell type can interact with which receptor
  #   receptors: vector of receptor names
  #   celltypes: vector of cell type names (alphabetically sorted)
  # -----------------------------
  rec_data <- createCellRecList(file.path(folder, "celltype_receptor.csv"))
  
  # -----------------------------
  # Read cell type distribution per dataset/patient
  # The function returns:
  #   Dtypes: matrix with relative abundance of each cell type per dataset (rows normalized)
  #   celltypes: vector of cell types
  #   datasets: vector of dataset/patient names
  # -----------------------------
  dist_data <- createCellTypeDistr(
    cells = lig_data$celltypes,
    filename = file.path(folder, paste0(cancer_name, "_TMEmod_cell_fractions.csv"))
  )
  
  # -----------------------------
  # Read ligand-receptor interaction weights
  # The function returns:
  #   DconnectionTensor: 3D array (ligands x receptors x datasets) of interaction probabilities
  # -----------------------------
  DconnectionTensor <- createInteractionDistr(
    filename = file.path(folder, paste0(cancer_name, "_LRpairs_weights_", weight_type, ".csv")),
    ligands = lig_data$ligands,
    receptors = rec_data$receptors
  )
  
  # -----------------------------
  # Optionally read the sign of interactions (+1 stimulating, -1 inhibiting, 0 unknown)
  # If read_signs = FALSE, create a matrix of zeros (unknown)
  # -----------------------------
  if (read_signs) { ### missing to check this
    sign_data <- Read_Lig_Rec_Interaction(
      file.path(current_path, folder,
                paste0(cancer_name, "_LRpairs_sign_interaction.csv"))
    )
    Sign_matrix <- sign_data$sign_matrix
  } else {
    Sign_matrix <- matrix(0,
                          nrow = dim(DconnectionTensor)[1],
                          ncol = dim(DconnectionTensor)[2])
  }
  
  # -----------------------------
  # Return all input data in a list for use in the network generation model
  # -----------------------------
  return(list(
    CellLigList = lig_data$CellLigList,       # cell-ligand compatibility matrix
    CellRecList = rec_data$CellRecList,       # cell-receptor compatibility matrix
    Dtypes = dist_data$Dtypes,                # normalized cell type distributions
    DconnectionTensor = DconnectionTensor,    # ligand-receptor interaction tensor
    celltypes = lig_data$celltypes,           # vector of cell type names
    ligands = lig_data$ligands,               # vector of ligand names
    receptors = rec_data$receptors,           # vector of receptor names
    Sign_matrix = Sign_matrix                  # matrix of interaction signs
  ))
}

# ================================
# 8. get_patient_names
# ================================
get_patient_names <- function(cancer_type,
                              folder = "Example input") {
  
  current_path <- getwd()
  
  lig_data <- createCellLigList(file.path(current_path, folder, "celltype_ligand.csv"))
  
  dist_data <- createCellTypeDistr(
    lig_data$celltypes,
    file.path(current_path, folder,
              paste0(cancer_type, "_TMEmod_cell_fractions.csv"))
  )
  
  return(dist_data$datasets)
}