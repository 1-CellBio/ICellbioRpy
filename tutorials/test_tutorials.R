# 教程功能测试脚本
# 这个脚本用于验证教程中的核心功能是否正常工作

cat("=== ICellbioRpy教程功能测试 ===\n")

# 加载包
suppressPackageStartupMessages({
  library(ICellbioRpy)
})

cat("✓ ICellbioRpy包加载成功\n")

# 测试1: Python环境配置
cat("\n--- 测试1: Python环境配置 ---\n")
tryCatch({
  configure_python_env(verbose = FALSE)
  cat("✓ Python环境配置成功\n")
}, error = function(e) {
  cat("✗ Python环境配置失败:", e$message, "\n")
})

# 测试2: 检查anndata可用性
cat("\n--- 测试2: anndata可用性检查 ---\n")
tryCatch({
  check_anndata_available()
  cat("✓ anndata可用性检查通过\n")
}, error = function(e) {
  cat("✗ anndata不可用:", e$message, "\n")
})

# 测试3: 函数可用性检查
cat("\n--- 测试3: 主要函数可用性 ---\n")

# 检查函数是否存在
functions_to_check <- c(
  "read1Cellbio",
  "iCellbio2H5ad", 
  "seurat_to_h5ad",
  "read_10x_mtx_to_h5ad",
  "configure_python_env",
  "check_anndata_available"
)

for (func_name in functions_to_check) {
  if (exists(func_name)) {
    cat("✓", func_name, "函数可用\n")
  } else {
    cat("✗", func_name, "函数不存在\n")
  }
}

# 测试4: 依赖包检查
cat("\n--- 测试4: 依赖包检查 ---\n")

required_packages <- c("Matrix", "data.table", "reticulate", "hdf5r", "jsonlite")
for (pkg in required_packages) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    cat("✓", pkg, "包可用\n")
  } else {
    cat("✗", pkg, "包未安装\n")
  }
}

# 测试5: 示例数据创建（用于测试）
cat("\n--- 测试5: 示例数据创建 ---\n")

tryCatch({
  # 创建示例CSV文件
  sample_csv <- tempfile(fileext = ".csv")
  
  # 创建虚拟文件路径（仅用于格式测试）
  sample_info <- data.frame(
    Sample_id = c("test1", "test2"),
    mtx_fns = c("/path/to/test1.mtx", "/path/to/test2.mtx"),
    features_fns = c("/path/to/test1_features.tsv", "/path/to/test2_features.tsv"),
    barcodes_fns = c("/path/to/test1_barcodes.tsv", "/path/to/test2_barcodes.tsv"),
    stringsAsFactors = FALSE
  )
  
  write.csv(sample_info, sample_csv, row.names = FALSE)
  
  # 验证CSV读取（不验证文件存在性）
  if (file.exists(sample_csv)) {
    csv_data <- read.csv(sample_csv, stringsAsFactors = FALSE)
    required_cols <- c("Sample_id", "mtx_fns", "features_fns", "barcodes_fns")
    
    if (all(required_cols %in% colnames(csv_data))) {
      cat("✓ CSV格式测试通过\n")
    } else {
      cat("✗ CSV格式测试失败\n")
    }
  }
  
  # 清理临时文件
  unlink(sample_csv)
  
}, error = function(e) {
  cat("✗ 示例数据创建失败:", e$message, "\n")
})

# 测试6: 文档可用性
cat("\n--- 测试6: 文档可用性检查 ---\n")

# 检查帮助文档
help_topics <- c(
  "read1Cellbio",
  "iCellbio2H5ad", 
  "seurat_to_h5ad",
  "read_10x_mtx_to_h5ad"
)

for (topic in help_topics) {
  tryCatch({
    help_file <- help(topic, package = "ICellbioRpy")
    if (length(help_file) > 0) {
      cat("✓", topic, "帮助文档可用\n")
    } else {
      cat("?", topic, "帮助文档可能不完整\n")
    }
  }, error = function(e) {
    cat("✗", topic, "帮助文档不可用\n")
  })
}

# 总结
cat("\n=== 测试完成 ===\n")
cat("如果所有测试都显示 ✓，说明教程环境准备就绪\n")
cat("如果有 ✗ 标记，请参考对应教程的环境配置部分\n")

# 环境信息
cat("\n=== 环境信息 ===\n")
cat("R版本:", R.version.string, "\n")
cat("ICellbioRpy版本:", as.character(packageVersion("ICellbioRpy")), "\n")
cat("系统:", Sys.info()["sysname"], Sys.info()["release"], "\n")

# Python信息（如果可用）
tryCatch({
  py_config <- reticulate::py_config()
  cat("Python版本:", py_config$version, "\n")
  cat("Python路径:", py_config$python, "\n")
}, error = function(e) {
  cat("Python信息不可用\n")
})

cat("\n测试脚本执行完毕！\n")