#' Build and save pre-computed lookup tables for all species
#'
#' This script pre-builds gene ID lookup tables for all supported species
#' to avoid the overhead of building them during GMT processing.
#'
#' @param data_dir Directory containing the mapping files
#' @param output_file Output RData file path
#' @return Invisible NULL
#' @export
prebuild_gene_lookup_tables <- function(data_dir = "data", output_file = "data/gene_lookup_tables.RData") {

  library(data.table)

  cat("开始预构建基因查找表...\n")

  # Define all supported species and mapping types
  species_list <- c("10090", "10116", "6239", "7227", "7955", "9598", "9606")
  mapping_types <- c("symbol", "entrez", "ensembl")

  # Initialize master lookup table
  master_lookup_tables <- list()

  for (species in species_list) {
    cat(sprintf("\n处理物种: %s\n", species))

    species_mappings <- list()

    for (type in mapping_types) {
      mapping_file <- file.path(data_dir, sprintf("%s_%s.tsv.gz", species, type))

      if (file.exists(mapping_file)) {
        cat(sprintf("  加载 %s 映射文件...\n", type))

        # Load mapping data
        mapping_data <- load_local_mapping_file(mapping_file, type)

        # Build lookup table for this mapping type
        lookup_table <- list()

        # Build gene_id -> gesel_id mapping
        for (i in 1:nrow(mapping_data)) {
          # Split identifiers in case there are multiple per row
          identifiers <- strsplit(as.character(mapping_data[i, 1]), "[\t ]")[[1]]
          identifiers <- identifiers[identifiers != ""]  # Remove empty strings

          gesel_id <- i - 1  # 0-based indexing

          # Add each identifier to lookup table
          for (gene_id in identifiers) {
            lookup_table[[gene_id]] <- gesel_id
          }
        }

        species_mappings[[type]] <- lookup_table
        cat(sprintf("    ✓ 构建了 %d 个基因ID的查找表\n", length(lookup_table)))

      } else {
        warning(sprintf("映射文件不存在: %s", mapping_file))
        species_mappings[[type]] <- NULL
      }
    }

    master_lookup_tables[[species]] <- species_mappings
  }

  # Save to RData file
  cat(sprintf("\n保存查找表到: %s\n", output_file))
  save(master_lookup_tables, file = output_file, compress = "gzip")
  cat("✓ 查找表保存完成\n")

  # Print summary
  cat("\n=== 查找表摘要 ===\n")
  for (species in names(master_lookup_tables)) {
    cat(sprintf("物种 %s:\n", species))
    for (type in names(master_lookup_tables[[species]])) {
      lookup_table <- master_lookup_tables[[species]][[type]]
      if (!is.null(lookup_table)) {
        cat(sprintf("  %s: %d 个基因ID\n", type, length(lookup_table)))
      } else {
        cat(sprintf("  %s: 未找到\n", type))
      }
    }
  }

  invisible(NULL)
}

#' Load pre-built lookup tables
#'
#' @param rdata_file Path to the RData file containing pre-built lookup tables
#' @return List of lookup tables
#' @export
load_gene_lookup_tables <- function(rdata_file = "data/gene_lookup_tables.RData") {
  if (!file.exists(rdata_file)) {
    stop("查找表文件不存在: ", rdata_file)
  }

  load(rdata_file)
  return(master_lookup_tables)
}

#' Get lookup table for specific species
#'
#' @param species Species identifier
#' @param rdata_file Path to the RData file
#' @return Lookup tables for the species
#' @export
get_species_lookup_table <- function(species, rdata_file = "data/gene_lookup_tables.RData") {
  master_tables <- load_gene_lookup_tables(rdata_file)

  if (!species %in% names(master_tables)) {
    stop("不支持的物种: ", species)
  }

  return(master_tables[[species]])
}

# Internal helper function (same as in process_gmt.R)
load_local_mapping_file <- function(file_path, type) {
  # Decompress and read using system command
  temp_file <- tempfile()
  system(sprintf("gunzip -c '%s' > '%s'", file_path, temp_file))
  # Read with flexible column handling
  data <- data.table::fread(temp_file, header = FALSE, sep = "\t", quote = "", fill = Inf)
  file.remove(temp_file)

  # Validate format
  validate_mapping_format(data, type)

  return(data)
}

# Internal helper function (same as in process_gmt.R)
validate_mapping_format <- function(data, type) {
  if (nrow(data) == 0) {
    stop(type, "映射文件为空")
  }

  # Check if it has enough columns
  if (ncol(data) < 1) {
    stop(type, "映射文件格式错误：至少需要1列")
  }

  cat("✓", type, "映射文件格式验证通过，包含", nrow(data), "个基因\n")
}
