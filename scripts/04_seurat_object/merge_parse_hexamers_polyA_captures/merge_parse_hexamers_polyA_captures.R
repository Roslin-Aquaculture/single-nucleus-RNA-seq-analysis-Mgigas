library(dplyr)
library(Matrix)
 
 
merge_hexamer_polyA_columns <- function(mtx = "",
                                        kit = "",
                                        version = "",
                                        bc_directory = ""
                                        ){
  #Find the correct barcode files - see info at bottom of this script.
  print("Loading barcodes")
  bc_paths <- handle_version(kit, version, bc_directory)
  print("Matching column names - polyA and hexamers")
  #Pair up the hexamer and polyA barcodes present in matrix
  matched_barcodes <- get_matching_barcodes(mtx, bc_paths)

  #Sum the matched hexamer and polyA columns
  summed_matrix <- sum_matching_columns(mtx, matched_barcodes)
  
  return(summed_matrix)
}  


sum_matching_columns <- function(mtx, matched_barcodes) {
  
  # Initialize a list to store columns for the new sparse matrix
  output_columns <- list()
  
  # Extract column names once for faster lookup
  mtx_columns <- colnames(mtx)
  
  # Get the total number of pairs for progress tracking
  total_pairs <- length(matched_barcodes)
  progress_intervals <- seq(0.1, 1.0, by = 0.1) * total_pairs
  next_progress <- progress_intervals[1]
  progress_index <- 1
  
  # Record the start time
  start_time <- Sys.time()
  
  # Loop through each pair of matched barcodes
  for (i in seq_along(matched_barcodes)) {
    pair <- matched_barcodes[[i]]
    column_1 <- pair[[1]]  # First column in the pair
    column_2 <- pair[[2]]  # Second column in the pair
    column_name <- pair[[3]]  # Name to assign to the output column
    
    # Check which columns exist in mtx
    col1_exists <- column_1 %in% mtx_columns
    col2_exists <- column_2 %in% mtx_columns
    
    # Process based on which columns are present
    if (col1_exists && col2_exists) {
      # Both columns exist: sum them and store with the specified name
      combined_column <- mtx[, column_1] + mtx[, column_2]
      output_columns[[column_name]] <- combined_column
    } else if (col1_exists) {
      # Only column 1 exists: store it with the original name
      output_columns[[column_1]] <- mtx[, column_1]
    } else if (col2_exists) {
      # Only column 2 exists: store it with the original name
      output_columns[[column_2]] <- mtx[, column_2]
    }
    
    # Update progress every 10%
    if (i >= next_progress) {
      # Calculate elapsed time in minutes
      elapsed_time <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
      # Estimate total time and remaining time in minutes
      estimated_total_time <- (elapsed_time / (progress_index * 0.1))
      remaining_time <- estimated_total_time - elapsed_time
      
      # Print only the remaining time in minutes
      cat(sprintf("Progress: %d%% completed. Estimated remaining time: %.2f minutes.\n",
                  progress_index * 10,
                  remaining_time))
      
      # Update to the next progress checkpoint
      progress_index <- progress_index + 1
      next_progress <- ifelse(progress_index <= length(progress_intervals), progress_intervals[progress_index], Inf)
    }
  }
  
  # Convert list to a sparse matrix
  output_matrix <- do.call(cbind, output_columns)
  output_sparse_matrix <- Matrix(output_matrix, sparse = TRUE)
  
  return(output_sparse_matrix)
}

#Function to generate list of polyA/hexamer barcode pairs.
get_matching_barcodes <- function(mtx, bc_paths){
  
  # Load the Parse barcode CSVs
  bc1_data <- read.csv(bc_paths[1])
  bc2_data <- read.csv(bc_paths[2])
  bc3_data <- read.csv(bc_paths[3])
  
  # Get list of column names i.e. the full cell barcode
  column_names <- colnames(mtx)
  
  # Initialize list to store all paired barcodes and a vector to track processed barcodes
  matched_barcodes <- list()
  processed_barcodes <- character()  # Character vector to track processed barcodes

  # Count number of hexamer and polyA captures - for interest
  polyA.count <- 0
  hexamer.count <- 0
  
  # Loop through each barcode
  for (column_name in column_names){
    
    # Third barcode in column name corresponds to first plate barcode
    all_bc_sequences <- strsplit(column_name, "_")[[1]]
    plate1_sequence <- all_bc_sequences[3]
    plate2_sequence <- all_bc_sequences[2]
    plate3_sequence <- all_bc_sequences[1]
    
    # Identify the wells and capture type
    first_well_identity <- bc1_data$well[bc1_data$sequence == plate1_sequence]
    second_well_identity <- bc2_data$well[bc2_data$sequence == plate2_sequence]
    third_well_identity <- bc3_data$well[bc3_data$sequence == plate3_sequence]
    capture_type <- bc1_data$stype[bc1_data$sequence == plate1_sequence]
    
    #Make a nice cell name to label the columns with later. BC1_BC2_BC3 order
    cell_name = paste(first_well_identity,second_well_identity,third_well_identity, sep = "_")
    
    #Add to polyA and hexamer counts
    if (capture_type == "R"){
      hexamer.count = hexamer.count + 1
    } else{
      polyA.count = polyA.count + 1
    }
    
    # Skip if this barcode has already been processed - doing now so counts above are correct
    if (column_name %in% processed_barcodes) {
      next
    }
    
    # Get the corresponding hexamer or PolyA sequence
    if (capture_type == "T"){
      matching_bc1_sequence <- bc1_data$sequence[bc1_data$well == first_well_identity & bc1_data$stype == "R"]
    } else if (capture_type == "R"){
      matching_bc1_sequence <- bc1_data$sequence[bc1_data$well == first_well_identity & bc1_data$stype == "T"]
    } else {
      stop("There is a bc1 with no corresponding 'other capture' sequence - something has gone wrong!")
    }
    
    # Construct full barcodes
    full_other_barcode <- paste(all_bc_sequences[1], all_bc_sequences[2], matching_bc1_sequence, sep = "_")
    full_current_barcode <- column_name
    
    # Append the matched barcodes as a named list entry
    matched_barcodes[[length(matched_barcodes) + 1]] <- list(
      polyA = if (capture_type == "T") full_current_barcode else full_other_barcode,
      hexamer = if (capture_type == "R") full_current_barcode else full_other_barcode,
      cell_name = cell_name
    )
    
    # Add both barcodes to processed_barcodes to avoid re-processing
    processed_barcodes <- c(processed_barcodes, full_current_barcode, full_other_barcode)
  }
  
  # Print counts
  print(paste0("Total polyA barcodes: ", polyA.count))
  print(paste0("Total hexamer barcodes: ", hexamer.count))
  
  return(matched_barcodes)
}

#Function to identify which barcode files to use depending on kit
handle_version <- function(kit, version,bc_directory) {
  # Define valid kits and versions
  valid_kits <- c("WT_mini", "WT", "WT_mega")
  valid_versions <- c("v1", "v2", "v3")
  
  # Check for valid input
  if (!(kit %in% valid_kits)) {
    stop("Invalid kit: please choose from 'WT_mini', 'WT', or 'WT_mega'")
  }
  if (!(version %in% valid_versions)) {
    stop("Invalid version: please choose from 'v1', 'v2', or 'v3'")
  }
  
  # Set file names based on kit and version
  if (kit == "WT_mini") {
    if (version == "v1" || version == "v2") {
      bc1_file <- "n24_v4"
      bc2_file <- "v1"
      bc3_file <- "v1"
    } else if (version == "v3") {
      bc1_file <- "n26_R1_v3_4"
      bc2_file <- "v1"
      bc3_file <- "R3_v3"
    }
  } else if (kit == "WT") {
    if (version == "v1") {
      bc1_file <- "v2"
      bc2_file <- "v1"
      bc3_file <- "v1"
    } else if (version == "v2") {
      bc1_file <- "n99_v5"
      bc2_file <- "v1"
      bc3_file <- "v1"
    } else if (version == "v3") {
      bc1_file <- "n107_R1_v3_4"
      bc2_file <- "v1"
      bc3_file <- "R3_v3"
    }
  } else if (kit == "WT_mega") {
    if (version == "v1" || version == "v2") {
      bc1_file <- "n198_v5"
      bc2_file <- "v1"
      bc3_file <- "v1"
    } else if (version == "v3") {
      bc1_file <- "n218_R1_v3_4"
      bc2_file <- "v1"
      bc3_file <- "R3_v3"
    }
  }
  
  # Create full file names
  bc1_file <- paste0(bc_directory,"/bc_data_", bc1_file, ".csv")
  bc2_file <- paste0(bc_directory,"/bc_data_", bc2_file, ".csv")
  bc3_file <- paste0(bc_directory,"/bc_data_", bc3_file, ".csv")
  
  # Return the file names as a vector
  return(c(bc1_file, bc2_file, bc3_file))
}

#Following is taken from ParseBiosciences-Pipeline.1.3.1 (kits.py)
#     kit         chem  rows    cols    plates  bc1           bc2     bc3    ktype
#     WT_mini     v1    1       12      1       n24_v4        v1      v1     normal
#     WT          v1    4       12      1       v2            v1      v1     normal
#     WT_mega     v1    8       12      1       n198_v5       v1      v1     normal
#     WT_mini     v2    1       12      1       n24_v4        v1      v1     normal
#     WT          v2    4       12      1       n99_v5        v1      v1     normal
#     WT_mega     v2    8       12      1       n198_v5       v1      v1     normal
#     WT_mini     v3    1       12      1       n26_R1_v3_4   v1      R3_v3  normal
#     WT          v3    4       12      1       n107_R1_v3_4  v1      R3_v3  normal
#     WT_mega     v3    8       12      1       n218_R1_v3_4  v1      R3_v3  normal