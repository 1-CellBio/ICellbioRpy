#' Preprocess GMT file with custom gene mapping support
#'
#' This function processes GMT (Gene Matrix Transposed) files and converts gene identifiers
#' to gesel format using custom mapping files. It supports multiple species and can
#' automatically download missing mapping files if needed.
#'
#' @param gmt_file Path to the GMT file to process
#' @param species Species identifier (default: "9606" for human)
#' @param gene_mapping_files Named list of custom mapping files with names: symbol, entrez, ensembl
#' @param collection_name Name for the gene set collection (auto-generated if NULL)
#' @param collection_desc Description of the gene set collection
#' @param output_dir Directory to save processed files (default: "processed_gesel")
#' @param auto_download_missing Whether to auto-download missing mapping files (default: TRUE)
#'
#' @return Invisible NULL. Creates processed files in the output directory.
#' @export
#'
#' @import data.table
#' @importFrom utils download.file
#'
#' @examples
#' \dontrun{
#' # Using custom mapping files
#' preprocess_gmt_with_custom_mapping(
#'   gmt_file = "path/to/genesets.gmt",
#'   species = "9606",
#'   gene_mapping_files = list(
#'     symbol = "path/to/symbol_mapping.tsv",
#'     entrez = "path/to/entrez_mapping.tsv"
#'   )
#' )
#'
#' # Using local data files (recommended)
#' preprocess_gmt_custom(
#'   gmt_file = "path/to/genesets.gmt",
#'   species = "9606",
#'   symbol_file = system.file("extdata", "9606_symbol.tsv.gz", package = "ICellbioRpy")
#' )
#' }
preprocess_gmt_with_custom_mapping <- function(
    gmt_file,
    species = "9606",
    gene_mapping_files = NULL,
    collection_name = NULL,
    collection_desc = "",
    output_dir = "processed_gesel",
    auto_download_missing = TRUE
) {

  # Create output directory
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  if (is.null(collection_name)) {
    collection_name <- tools::file_path_sans_ext(basename(gmt_file))
  }

  cat("开始处理 GMT 文件:", gmt_file, "\n")
  cat("物种:", species, "\n")

  # 1. Read GMT file
  gmt_data <- data.table::fread(gmt_file, header = FALSE, sep = "\t", quote = "", fill = Inf)
  cat("读取了", nrow(gmt_data), "个基因集\n")

  # 2. Parse GMT data
  parsed_gmt <- parse_gmt_data(gmt_data)

  # 3. Check if pre-built lookup tables are available
  prebuilt_available <- exists("master_lookup_tables", inherits = TRUE)
  if (prebuilt_available) {
    lookup_tables_obj <- get("master_lookup_tables", inherits = TRUE)
    if (species %in% names(lookup_tables_obj)) {
      cat("✓ 发现预构建查找表，将跳过映射文件加载\n")
    } else {
      prebuilt_available <- FALSE
      cat("预构建查找表中没有找到物种", species, "\n")
    }
  }

  # 4. Load gene mappings only if pre-built tables are not available
  gene_mappings <- NULL
  if (!prebuilt_available) {
    cat("加载基因映射...\n")
    gene_mappings <- load_gene_mappings(species, gene_mapping_files, auto_download_missing)
  } else {
    cat("使用预构建查找表，跳过映射文件加载...\n")
  }

  # 5. Convert gene IDs
  cat("转换基因ID...\n")
  converted_lists <- convert_gene_lists_to_gesel(parsed_gmt$gene_lists, gene_mappings, species)

  # 5. Generate gesel format files
  generate_gesel_files(
    converted_lists,
    parsed_gmt$set_names,
    parsed_gmt$set_descriptions,
    species,
    collection_name,
    collection_desc,
    output_dir
  )

  cat("处理完成！\n")
  invisible(NULL)
}

#' Parse GMT data
#'
#' @param gmt_data Data.table containing GMT file contents
#' @return List with set_names, set_descriptions, and gene_lists
#' @keywords internal
parse_gmt_data <- function(gmt_data) {
  set_names <- gmt_data[[1]]

  # 修复：逐行处理描述字段，避免ifelse的NA传播问题
  set_descriptions <- character(nrow(gmt_data))
  if (ncol(gmt_data) >= 2) {
    for (i in 1:nrow(gmt_data)) {
      desc <- gmt_data[i, 2, with = FALSE]
      set_descriptions[i] <- if (is.na(desc) || desc == "") "" else as.character(desc)
    }
  } else {
    set_descriptions <- rep("", nrow(gmt_data))
  }

  gene_lists <- list()
  for (i in 1:nrow(gmt_data)) {
    genes <- c()
    for (col_idx in 3:ncol(gmt_data)) {
      gene <- trimws(as.character(gmt_data[i, col_idx, with = FALSE]))
      if (!is.na(gene) && gene != "") {
        genes <- c(genes, gene)
      }
    }
    gene_lists[[i]] <- unique(genes[genes != ""])
  }

  # Filter empty gene sets
  valid_idx <- sapply(gene_lists, length) > 0
  list(
    set_names = set_names[valid_idx],
    set_descriptions = set_descriptions[valid_idx],
    gene_lists = gene_lists[valid_idx]
  )
}

#' Load gene mapping files (legacy function - now mainly for compatibility)
#'
#' @param species Species identifier
#' @param custom_files Named list of custom mapping files
#' @param auto_download Whether to auto-download missing files
#' @return List of loaded mapping data.tables
#' @keywords internal
load_gene_mappings <- function(species, custom_files = NULL, auto_download = TRUE) {
  # Since we now use pre-built lookup tables, this function is mainly for compatibility
  # It will still work if custom files are provided, but will warn about the inefficiency

  if (!is.null(custom_files)) {
    warning("自定义映射文件已提供，但建议使用预构建查找表以获得最佳性能")
  }

  mappings <- list()
  types <- c("symbol", "entrez", "ensembl")

  for (type in types) {
    file_loaded <- FALSE

    # Priority: use custom files
    if (!is.null(custom_files) && !is.null(custom_files[[type]])) {
      cat("使用自定义", type, "映射文件:", custom_files[[type]], "\n")
      tryCatch({
        mapping_data <- load_custom_mapping_file(custom_files[[type]], type)
        mappings[[type]] <- mapping_data
        file_loaded <- TRUE
        cat("✓ 成功加载自定义", type, "映射文件\n")
      }, error = function(e) {
        warning("无法加载自定义", type, "映射文件:", e$message)
      })
    }

    # Note: We no longer try to load local .tsv.gz files since they are excluded from the package
    # Users should use pre-built lookup tables instead

    if (!file_loaded) {
      warning("无法获取", type, "映射文件，将跳过此类型的映射。建议使用预构建查找表。")
    }
  }

  return(mappings)
}

#' Load custom mapping file
#'
#' @param file_path Path to the mapping file
#' @param type Type of mapping (symbol, entrez, ensembl)
#' @return Data.table with mapping data
#' @keywords internal
load_custom_mapping_file <- function(file_path, type) {
  # Check if file exists
  if (!file.exists(file_path)) {
    stop("映射文件不存在: ", file_path)
  }

  # Check if it's compressed
  if (grepl("\\.gz$", file_path)) {
  # Decompress and read
  if (requireNamespace("R.utils", quietly = TRUE)) {
    temp_file <- tempfile()
    R.utils::gunzip(file_path, temp_file, remove = FALSE)
    data <- data.table::fread(temp_file, header = FALSE, sep = "\t", quote = "")
    file.remove(temp_file)
    } else {
      # Use system command to decompress
      temp_file <- tempfile()
      system(sprintf("gunzip -c '%s' > '%s'", file_path, temp_file))
      data <- data.table::fread(temp_file, header = FALSE, sep = "\t", quote = "")
      file.remove(temp_file)
    }
  } else {
    # Read directly
    data <- data.table::fread(file_path, header = FALSE, sep = "\t", quote = "")
  }

  # Validate format
  validate_mapping_format(data, type)

  return(data)
}

#' Load local mapping file from package data
#'
#' @param file_path Path to the local mapping file
#' @param type Type of mapping (symbol, entrez, ensembl)
#' @return Data.table with mapping data
#' @keywords internal
load_local_mapping_file <- function(file_path, type) {
  # Decompress and read
  if (requireNamespace("R.utils", quietly = TRUE)) {
    temp_file <- tempfile()
    R.utils::gunzip(file_path, temp_file, remove = FALSE)
    # Read with flexible column handling
    data <- data.table::fread(temp_file, header = FALSE, sep = "\t", quote = "", fill = TRUE)
    file.remove(temp_file)
  } else {
    # Use system command to decompress
    temp_file <- tempfile()
    system(sprintf("gunzip -c '%s' > '%s'", file_path, temp_file))
    # Read with flexible column handling
    data <- data.table::fread(temp_file, header = FALSE, sep = "\t", quote = "", fill = TRUE)
    file.remove(temp_file)
  }

  # Validate format
  validate_mapping_format(data, type)

  return(data)
}

#' Validate mapping file format
#'
#' @param data Data.table to validate
#' @param type Type of mapping
#' @keywords internal
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

#' Download gene mapping file
#'
#' @param species Species identifier
#' @param type Type of mapping (symbol, entrez, ensembl)
#' @return Data.table with downloaded mapping data
#' @keywords internal
download_gene_mapping <- function(species, type) {
  base_url <- "https://github.com/LTLA/gesel-feedstock/releases/download/indices-v0.2.2"
  file_name <- sprintf("%s_%s.tsv.gz", species, type)
  url <- file.path(base_url, file_name)

  temp_file <- tempfile()
  tryCatch({
    utils::download.file(url, temp_file, quiet = TRUE)

    # Decompress
    if (requireNamespace("R.utils", quietly = TRUE)) {
      unzipped_file <- tempfile()
      R.utils::gunzip(temp_file, unzipped_file, remove = FALSE)
      data <- data.table::fread(unzipped_file, header = FALSE, sep = "\t", quote = "")
      file.remove(unzipped_file)
    } else {
      unzipped_file <- tempfile()
      system(sprintf("gunzip -c '%s' > '%s'", temp_file, unzipped_file))
      data <- data.table::fread(unzipped_file, header = FALSE, sep = "\t", quote = "")
      file.remove(unzipped_file)
    }

    file.remove(temp_file)
    return(data)

  }, error = function(e) {
    file.remove(temp_file)
    stop("下载失败: ", e$message)
  })
}

#' Convert gene lists to gesel IDs
#'
#' @param gene_lists List of gene lists
#' @param mappings List of mapping data.tables
#' @param species Species identifier for loading pre-built lookup tables
#' @return List of converted gene lists with gesel IDs
#' @keywords internal
convert_gene_lists_to_gesel <- function(gene_lists, mappings, species = NULL) {
  # Try to use pre-built lookup tables first
  lookup_tables <- NULL

  if (!is.null(species)) {
    # Check if master_lookup_tables is available (loaded by package lazy-load)
    if (exists("master_lookup_tables", inherits = TRUE)) {
      cat("Loading pre-built lookup tables...\n")
      tryCatch({
        lookup_tables_obj <- get("master_lookup_tables", inherits = TRUE)
        if (species %in% names(lookup_tables_obj)) {
          lookup_tables <- lookup_tables_obj[[species]]
          cat("✓ Successfully loaded pre-built lookup tables for species:", species, "\n")
        } else {
          cat("Pre-built lookup tables found but species", species, "not available, falling back to dynamic build\n")
        }
      }, error = function(e) {
        cat("Failed to access pre-built lookup tables:", e$message, "\n")
      })
    } else {
      cat("Pre-built lookup tables not available, falling back to dynamic build\n")
    }
  }

  # If pre-built tables not available, build them dynamically
  if (is.null(lookup_tables)) {
    if (is.null(mappings)) {
      stop("预构建查找表不可用且未提供映射文件，无法进行基因ID转换")
    }
    cat("Building lookup tables...\n")
    lookup_tables <- build_lookup_tables(mappings)
  }

  converted_lists <- list()

  for (i in 1:length(gene_lists)) {
    genes <- gene_lists[[i]]
    if (length(genes) == 0) {
      converted_lists[[i]] <- c()
      next
    }

    gesel_ids <- c()
    for (gene in genes) {
      id <- find_gesel_id_from_lookup(gene, lookup_tables)
      if (!is.null(id)) {
        gesel_ids <- c(gesel_ids, id)
      }
    }

    converted_lists[[i]] <- unique(gesel_ids)
    cat(sprintf("\r处理基因集 %d/%d: %d个基因 -> %d个gesel ID",
                i, length(gene_lists), length(genes), length(gesel_ids)))
  }
  cat("\n")

  return(converted_lists)
}

#' Build lookup tables for efficient gene ID searching
#'
#' @param mappings List of mapping data.tables
#' @return List of lookup tables (gene_id -> gesel_id)
#' @keywords internal
build_lookup_tables <- function(mappings) {
  lookup_tables <- list()

  for (type in names(mappings)) {
    mapping <- mappings[[type]]
    if (is.null(mapping)) next

    lookup_table <- list()

    # Build gene_id -> gesel_id mapping
    for (i in 1:nrow(mapping)) {
      # Split identifiers in case there are multiple per row
      identifiers <- strsplit(as.character(mapping[i, 1]), "[\t ]")[[1]]
      identifiers <- identifiers[identifiers != ""]  # Remove empty strings

      gesel_id <- i - 1  # 0-based indexing

      # Add each identifier to lookup table
      for (gene_id in identifiers) {
        lookup_table[[gene_id]] <- gesel_id
      }
    }

    lookup_tables[[type]] <- lookup_table
    cat(sprintf("Built lookup table for %s: %d gene IDs\n", type, length(lookup_table)))
  }

  return(lookup_tables)
}

#' Find gesel ID from lookup tables (efficient version)
#'
#' @param gene Gene identifier to search for
#' @param lookup_tables Pre-built lookup tables
#' @return Gesel ID (0-based) or NULL if not found
#' @keywords internal
find_gesel_id_from_lookup <- function(gene, lookup_tables) {
  # Try all mapping types
  for (type in names(lookup_tables)) {
    lookup_table <- lookup_tables[[type]]
    if (!is.null(lookup_table[[gene]])) {
      return(lookup_table[[gene]])
    }
  }

  return(NULL)
}

#' Find gesel ID from mappings (original slow version, kept for compatibility)
#'
#' @param gene Gene identifier to search for
#' @param mappings List of mapping data.tables
#' @return Gesel ID (0-based) or NULL if not found
#' @keywords internal
find_gesel_id_from_mappings <- function(gene, mappings) {
  # Try all mapping types
  for (type in names(mappings)) {
    mapping <- mappings[[type]]
    if (is.null(mapping)) next

    # Search in current mapping
    for (i in 1:nrow(mapping)) {
      identifiers <- strsplit(as.character(mapping[i, 1]), "[\t ]")[[1]]
      identifiers <- identifiers[identifiers != ""]  # Remove empty strings
      if (gene %in% identifiers) {
        return(i - 1)  # gesel IDs are 0-based
      }
    }
  }

  return(NULL)
}

#' 生成gesel格式的完整文件集
#'
#' @param converted_lists List of converted gene lists
#' @param set_names Vector of gene set names
#' @param set_descriptions Vector of gene set descriptions
#' @param species Species identifier
#' @param collection_name Name of the collection
#' @param collection_desc Description of the collection
#' @param output_dir Output directory
#' @keywords internal
generate_gesel_files <- function(converted_lists, set_names, set_descriptions,
                                species, collection_name, collection_desc, output_dir) {

  # 1. 生成collections.tsv.gz
  cat("生成 collections.tsv.gz...\n")
  generate_collections_file(
    collection_name, collection_desc, species,
    length(set_names), output_dir
  )

  # 2. 生成sets.tsv.gz
  cat("生成 sets.tsv.gz...\n")
  generate_sets_file(
    set_names, set_descriptions,
    sapply(converted_lists, length), species, output_dir
  )

  # 3. 生成set2gene.tsv.gz
  cat("生成 set2gene.tsv.gz...\n")
  generate_set2gene_file(converted_lists, species, output_dir)

  # 4. 生成gene2set.tsv.gz
  cat("生成 gene2set.tsv.gz...\n")
  generate_gene2set_file(converted_lists, species, output_dir)

  # 5. 生成summary信息
  cat("生成处理摘要...\n")
  generate_summary_report(
    set_names, converted_lists, species,
    collection_name, output_dir
  )
}

#' Simplified user interface for GMT preprocessing with custom mapping
#'
#' @param gmt_file Path to GMT file
#' @param species Species identifier (default: "9606")
#' @param symbol_file Optional path to symbol mapping file
#' @param entrez_file Optional path to Entrez ID mapping file
#' @param ensembl_file Optional path to Ensembl ID mapping file
#' @param collection_name Optional name for the collection
#' @param output_dir Output directory (default: "processed_gesel")
#'
#' @return Invisible NULL
#' @export
#'
#' @examples
#' \dontrun{
#' # Using local package data files
#' preprocess_gmt_custom(
#'   gmt_file = "path/to/genesets.gmt",
#'   species = "9606",
#'   symbol_file = system.file("data", "9606_symbol.tsv.gz", package = "ICellbioRpy"),
#'   entrez_file = system.file("data", "9606_entrez.tsv.gz", package = "ICellbioRpy")
#' )
#' }
preprocess_gmt_custom <- function(
    gmt_file,
    species = "9606",
    symbol_file = NULL,
    entrez_file = NULL,
    ensembl_file = NULL,
    collection_name = NULL,
    output_dir = "processed_gesel"
) {

  # Build mapping files list
  gene_mapping_files <- list()
  if (!is.null(symbol_file)) gene_mapping_files$symbol <- symbol_file
  if (!is.null(entrez_file)) gene_mapping_files$entrez <- entrez_file
  if (!is.null(ensembl_file)) gene_mapping_files$ensembl <- ensembl_file

  # Call main function
  preprocess_gmt_with_custom_mapping(
    gmt_file = gmt_file,
    species = species,
    gene_mapping_files = gene_mapping_files,
    collection_name = collection_name,
    output_dir = output_dir,
    auto_download_missing = length(gene_mapping_files) == 0  # Auto-download if no custom files provided
  )
}

#' 生成collections文件
#'
#' @param name Collection name
#' @param desc Collection description
#' @param species Species identifier
#' @param size Number of gene sets
#' @param output_dir Output directory
#' @keywords internal
generate_collections_file <- function(name, desc, species, size, output_dir) {
  # Ensure output directory exists
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  collections_data <- data.frame(
    title = name,
    description = desc,
    species = species,
    maintainer = "custom",
    source = "local",
    size = size
  )

  file_path <- file.path(output_dir, sprintf("%s_collections.tsv", species))
  write.table(collections_data, file_path, sep = "\t",
              row.names = FALSE, col.names = FALSE, quote = FALSE)

  # 压缩 (使用system调用gzip)
  system(sprintf("gzip '%s'", file_path))
  cat("✓ 生成:", basename(file_path), ".gz\n")
}

#' 生成sets文件
#'
#' @param names Gene set names
#' @param descriptions Gene set descriptions
#' @param sizes Gene set sizes
#' @param species Species identifier
#' @param output_dir Output directory
#' @keywords internal
generate_sets_file <- function(names, descriptions, sizes, species, output_dir) {
  sets_data <- data.frame(
    name = names,
    description = descriptions,
    size = sizes
  )

  file_path <- file.path(output_dir, sprintf("%s_sets.tsv", species))
  write.table(sets_data, file_path, sep = "\t",
              row.names = FALSE, col.names = FALSE, quote = FALSE)

  system(sprintf("gzip '%s'", file_path))
  cat("✓ 生成:", basename(file_path), ".gz\n")
}

#' 生成set2gene文件（差分编码）
#'
#' @param converted_lists List of converted gene lists
#' @param species Species identifier
#' @param output_dir Output directory
#' @keywords internal
generate_set2gene_file <- function(converted_lists, species, output_dir) {
  # 为每个基因集生成差分编码的基因ID列表
  set2gene_lines <- sapply(converted_lists, function(genes) {
    if (length(genes) == 0) return("")

    # 排序并进行差分编码
    sorted_genes <- sort(unique(genes))
    encode_delta_line(sorted_genes)
  })

  # 移除空行
  valid_lines <- set2gene_lines[set2gene_lines != ""]

  file_path <- file.path(output_dir, sprintf("%s_set2gene.tsv", species))
  writeLines(valid_lines, file_path)

  R.utils::gzip(file_path)
  cat("✓ 生成:", basename(file_path), ".gz (", length(valid_lines), "个基因集)\n")
}

#' 生成gene2set文件
#'
#' @param converted_lists List of converted gene lists
#' @param species Species identifier
#' @param output_dir Output directory
#' @keywords internal
generate_gene2set_file <- function(converted_lists, species, output_dir) {
  # 找到所有唯一的基因ID
  all_genes <- unique(unlist(converted_lists))
  all_genes <- sort(all_genes)

  cat("处理", length(all_genes), "个唯一基因...\n")

  # 使用环境而不是大向量来构建映射，这样更内存高效
  gene_to_sets <- new.env(hash = TRUE, size = length(all_genes))

  for (set_idx in 1:length(converted_lists)) {
    genes <- converted_lists[[set_idx]]
    for (gene in genes) {
      gene_key <- as.character(gene)
      if (is.null(gene_to_sets[[gene_key]])) {
        gene_to_sets[[gene_key]] <- c()
      }
      gene_to_sets[[gene_key]] <- c(gene_to_sets[[gene_key]], set_idx - 1)  # 0-based
    }
  }

  # 生成差分编码行
  gene2set_lines <- character(length(all_genes))
  for (i in 1:length(all_genes)) {
    gene_idx <- all_genes[i]
    gene_key <- as.character(gene_idx)
    sets <- gene_to_sets[[gene_key]]

    if (!is.null(sets) && length(sets) > 0) {
      gene2set_lines[i] <- encode_delta_line(sort(unique(sets)))
    } else {
      gene2set_lines[i] <- ""
    }
  }

  # 只写入非空行
  valid_lines <- gene2set_lines[gene2set_lines != ""]

  file_path <- file.path(output_dir, sprintf("%s_gene2set.tsv", species))
  writeLines(valid_lines, file_path)

  R.utils::gzip(file_path)
  cat("✓ 生成:", basename(file_path), ".gz (", length(valid_lines), "个基因有映射)\n")
}

#' 差分编码函数
#'
#' @param arr Array of integers to encode
#' @return Delta-encoded string
#' @keywords internal
encode_delta_line <- function(arr) {
  if (length(arr) == 0) return("")
  arr <- sort(unique(arr))
  if (length(arr) == 1) return(as.character(arr[1]))

  deltas <- c(arr[1], diff(arr))
  paste(deltas, collapse = "\t")
}

#' 生成处理摘要报告
#'
#' @param set_names Gene set names
#' @param converted_lists Converted gene lists
#' @param species Species identifier
#' @param collection_name Collection name
#' @param output_dir Output directory
#' @keywords internal
generate_summary_report <- function(set_names, converted_lists, species, collection_name, output_dir) {
  total_sets <- length(set_names)
  total_genes_original <- sum(sapply(converted_lists, length))
  unique_genes <- length(unique(unlist(converted_lists)))

  # 计算基因集大小分布
  set_sizes <- sapply(converted_lists, length)

  report <- sprintf("
基因集预处理报告
==================

集合信息:
- 集合名称: %s
- 物种: %s
- 处理时间: %s

统计信息:
- 基因集总数: %d
- 原始基因总数: %d
- 唯一基因数: %d
- 平均基因集大小: %.1f

基因集大小分布:
- 最小: %d 基因
- 最大: %d 基因
- 中位数: %d 基因

输出文件:
- %s_collections.tsv.gz
- %s_sets.tsv.gz
- %s_set2gene.tsv.gz
- %s_gene2set.tsv.gz

所有文件已生成并压缩为gzip格式。
",
    collection_name,
    species,
    format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    total_sets,
    total_genes_original,
    unique_genes,
    mean(set_sizes),
    min(set_sizes),
    max(set_sizes),
    median(set_sizes),
    species, species, species, species
  )

  report_file <- file.path(output_dir, "processing_report.txt")
  writeLines(report, report_file)
  cat("✓ 生成处理报告:", basename(report_file), "\n")
}
