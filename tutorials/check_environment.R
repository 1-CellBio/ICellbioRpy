# ICellbioRpy环境检测脚本
# 这个脚本帮助初学者快速检查是否需要进行环境安装

cat("🔍 ICellbioRpy环境检测器\n")
cat("============================\n")

# 检测操作系统
os_info <- Sys.info()["sysname"]
cat("检测到的操作系统:", os_info, "\n\n")

# 初始化检查结果
need_setup <- FALSE
issues <- c()

# 1. 检查R版本
cat("1. 检查R版本...\n")
r_version <- R.version.string
r_version_num <- as.numeric(R.version$major) + as.numeric(R.version$minor) / 10

if (!is.na(r_version_num) && r_version_num >= 4.0) {
  cat("   ✅ R版本:", r_version, "\n")
} else {
  cat("   ❌ R版本过低:", r_version, "（需要 >= 4.0.0）\n")
  need_setup <- TRUE
  issues <- c(issues, "R版本需要更新")
}

# 2. 检查是否在RStudio中运行
cat("\n2. 检查运行环境...\n")
if (Sys.getenv("RSTUDIO") == "1") {
  cat("   ✅ 正在RStudio中运行\n")
} else {
  cat("   ⚠️  建议使用RStudio进行学习\n")
}

# 3. 检查关键R包
cat("\n3. 检查必需的R包...\n")
required_packages <- c("devtools", "reticulate", "Matrix", "data.table", "hdf5r", "jsonlite")

missing_packages <- c()
for (pkg in required_packages) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    cat("   ✅", pkg, "\n")
  } else {
    cat("   ❌", pkg, "未安装\n")
    missing_packages <- c(missing_packages, pkg)
    need_setup <- TRUE
  }
}

if (length(missing_packages) > 0) {
  issues <- c(issues, paste("缺少R包:", paste(missing_packages, collapse = ", ")))
}

# 4. 检查ICellbioRpy包
cat("\n4. 检查ICellbioRpy包...\n")
if (requireNamespace("ICellbioRpy", quietly = TRUE)) {
  cat("   ✅ ICellbioRpy已安装\n")
  
  # 尝试加载并检查主要函数
  tryCatch({
    library(ICellbioRpy, quietly = TRUE)
    main_functions <- c("read1Cellbio", "iCellbio2H5ad", "configure_python_env")
    
    missing_functions <- c()
    for (func in main_functions) {
      if (exists(func)) {
        cat("   ✅", func, "函数可用\n")
      } else {
        cat("   ❌", func, "函数不可用\n")
        missing_functions <- c(missing_functions, func)
      }
    }
    
    if (length(missing_functions) > 0) {
      issues <- c(issues, "ICellbioRpy函数不完整")
      need_setup <- TRUE
    }
    
  }, error = function(e) {
    cat("   ❌ ICellbioRpy加载失败:", e$message, "\n")
    issues <- c(issues, "ICellbioRpy无法正常加载")
    need_setup <- TRUE
  })
  
} else {
  cat("   ❌ ICellbioRpy未安装\n")
  need_setup <- TRUE
  issues <- c(issues, "ICellbioRpy包未安装")
}

# 5. 检查Python环境
cat("\n5. 检查Python环境...\n")

# 检查reticulate是否可用
if (requireNamespace("reticulate", quietly = TRUE)) {
  tryCatch({
    library(reticulate, quietly = TRUE)
    
    # 检查Python是否可用
    py_available <- py_available()
    if (py_available) {
      py_config <- py_config()
      cat("   ✅ Python可用，版本:", py_config$version, "\n")
      cat("   路径:", py_config$python, "\n")
      
      # 检查conda
      conda_available <- tryCatch({
        conda_list()
        TRUE
      }, error = function(e) FALSE)
      
      if (conda_available) {
        cat("   ✅ conda环境管理器可用\n")
        
        # 检查是否有1cellbio环境
        envs <- conda_list()
        if ("1cellbio" %in% envs$name) {
          cat("   ✅ 1cellbio conda环境已存在\n")
        } else {
          cat("   ⚠️  1cellbio conda环境不存在\n")
          issues <- c(issues, "需要创建1cellbio conda环境")
        }
      } else {
        cat("   ❌ conda不可用\n")
        need_setup <- TRUE
        issues <- c(issues, "conda环境管理器不可用")
      }
      
    } else {
      cat("   ❌ Python不可用\n")
      need_setup <- TRUE
      issues <- c(issues, "Python环境不可用")
    }
    
  }, error = function(e) {
    cat("   ❌ Python环境检查失败:", e$message, "\n")
    need_setup <- TRUE
    issues <- c(issues, "Python环境配置有问题")
  })
} else {
  cat("   ❌ reticulate包未安装\n")
  need_setup <- TRUE
  issues <- c(issues, "reticulate包未安装")
}

# 6. 检查anndata（如果Python可用）
if (requireNamespace("reticulate", quietly = TRUE) && py_available()) {
  cat("\n6. 检查anndata包...\n")
  tryCatch({
    library(ICellbioRpy, quietly = TRUE)
    check_anndata_available()
    cat("   ✅ anndata包可用\n")
  }, error = function(e) {
    cat("   ❌ anndata不可用:", e$message, "\n")
    need_setup <- TRUE
    issues <- c(issues, "anndata包不可用")
  })
}

# 生成总结报告
cat("\n", rep("=", 50), "\n")
cat("📊 环境检测总结\n")
cat(rep("=", 50), "\n")

if (!need_setup) {
  cat("🎉 恭喜！您的环境配置完整，可以直接开始使用教程！\n\n")
  cat("推荐的学习路径：\n")
  cat("1. 教程1: 1CellBio转H5AD\n")
  cat("2. 教程2: 1CellBio转Seurat\n")
  cat("3. 教程3: Seurat转H5AD\n")
  cat("4. 教程4: 10X MTX整合\n")
} else {
  cat("⚠️  检测到以下问题需要解决：\n")
  for (i in seq_along(issues)) {
    cat("   ", i, ".", issues[i], "\n")
  }
  
  cat("\n🔧 建议的解决方案：\n")
  
  if (os_info == "Windows") {
    cat("📘 请完成：[教程0a: Windows环境完整安装指南](tutorial_0_environment_setup_windows.Rmd)\n")
  } else if (os_info == "Darwin") {
    cat("📘 请完成：[教程0b: macOS环境完整安装指南](tutorial_0_environment_setup_macos.Rmd)\n")
  } else {
    cat("📘 请根据您的系统选择合适的环境安装教程\n")
  }
  
  cat("\n完成环境安装后，请重新运行此脚本进行验证。\n")
}

cat("\n💡 小贴士：\n")
cat("- 如果您是完全初学者，强烈建议完成环境安装教程\n")
cat("- 环境安装教程包含了详细的故障排除指南\n")
cat("- 遇到问题时，请查看对应教程的问题排查部分\n")

cat("\n🆘 需要帮助？\n")
cat("- 查看tutorials/README.md获取完整教程列表\n")
cat("- 访问GitHub项目页面寻求支持\n")
cat("- 联系维护者获取帮助\n")

cat("\n检测完成！\n")