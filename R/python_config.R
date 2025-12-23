#' Configure Python environment for anndata (REQUIRED)
#'
#' This function configures the Python environment to use existing anndata installation.
#' IMPORTANT: This function MUST be called before using any conversion functions.
#'
#' @param python_path Optional path to Python executable. If NULL, uses current environment.
#' @param conda_env Optional conda environment name. If NULL, uses current environment.
#' @param verbose Logical, whether to print configuration messages (default: TRUE)
#' @param timeout Numeric, timeout in seconds for Python configuration (default: 30)
#' @param force_conda Logical, whether to force conda environment usage (default: TRUE)
#'
#' @return Invisible TRUE if configuration successful, throws error if failed
#'
#' @examples
#' \dontrun{
#' # REQUIRED: Configure before any operations
#' configure_python_env(conda_env = "1cellbio", verbose = TRUE)
#' 
#' # Use specific conda environment (recommended)
#' configure_python_env(conda_env = "your_env_name")
#' 
#' # Use specific Python path (advanced)
#' configure_python_env(python_path = "/usr/local/bin/python3")
#' }
#'
#' @export
configure_python_env <- function(python_path = NULL, conda_env = NULL, verbose = TRUE, 
                                timeout = 30, force_conda = TRUE) {
  
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop("Package 'reticulate' is required but not installed. Please install it with: install.packages('reticulate')")
  }
  
  # Force environment settings to prevent automatic configuration
  Sys.setenv(RETICULATE_AUTOCONFIGURE = "FALSE")
  Sys.setenv(RETICULATE_MINICONDA_ENABLED = "FALSE")
  Sys.setenv(RETICULATE_PYTHON_FALLBACK = "FALSE")
  
  if (verbose) cat(icb_i18n("正在配置Python环境...\n", "Configuring Python environment...\n"))
  
  # Step 1: Configure Python environment with timeout protection
  config_result <- tryCatch({
    # Wrap in timeout to prevent hanging
    R.utils::withTimeout({
      if (!is.null(conda_env)) {
        if (verbose) cat(icb_i18n(
          paste0("使用指定的conda环境: ", conda_env, "\n"),
          paste0("Using specified conda environment: ", conda_env, "\n")
        ))
        # Force conda environment configuration
        reticulate::use_condaenv(conda_env, required = TRUE)
      } else if (!is.null(python_path)) {
        if (verbose) cat(icb_i18n(
          paste0("使用指定的Python路径: ", python_path, "\n"),
          paste0("Using specified Python path: ", python_path, "\n")
        ))
        reticulate::use_python(python_path, required = TRUE)
      } else if (force_conda) {
        # If no specific environment given but force_conda is TRUE, try common environments
        common_envs <- c("1cellbio", "base", "scanpy", "anndata")
        env_found <- FALSE
        for (env in common_envs) {
          tryCatch({
            reticulate::use_condaenv(env, required = FALSE)
            if (verbose) cat(icb_i18n(
              paste0("尝试使用conda环境: ", env, "\n"),
              paste0("Trying conda environment: ", env, "\n")
            ))
            env_found <- TRUE
            break
          }, error = function(e) NULL)
        }
        if (!env_found) {
          stop(icb_i18n(
            "未找到可用的conda环境。请指定环境：configure_python_env(conda_env = \"1cellbio\")",
            "No available conda environment found. Please specify: configure_python_env(conda_env = \"1cellbio\")"
          ))
        }
      } else {
        if (verbose) cat(icb_i18n("使用当前Python环境...\n", "Using current Python environment...\n"))
      }
      TRUE
    }, timeout = timeout)
  }, TimeoutException = function(e) {
    stop(icb_i18n(
      paste0("Python环境配置超时（", timeout, "秒）。请检查conda环境是否存在。"),
      paste0("Python environment configuration timed out (", timeout, " seconds). Please check if conda environment exists.")
    ))
  }, error = function(e) {
    stop(icb_i18n(
      paste0("Python环境配置失败: ", e$message),
      paste0("Python environment configuration failed: ", e$message)
    ))
  })
  
  # Step 2: Quick validation with timeout
  validation_result <- tryCatch({
    R.utils::withTimeout({
      # Quick check without full py_config() which can hang
      py_available <- reticulate::py_available(initialize = TRUE)
      if (!py_available) {
        stop(icb_i18n("Python不可用", "Python not available"))
      }
      
      # Quick anndata availability check
      anndata_available <- reticulate::py_module_available("anndata")
      if (!anndata_available) {
        stop(icb_i18n(
          paste0("在指定环境中未找到anndata模块。\n",
                "请安装: pip install anndata\n",
                "或指定已安装anndata的环境: configure_python_env(conda_env = \"1cellbio\")"),
          paste0("anndata module not found in specified environment.\n",
                "Please install: pip install anndata\n", 
                "Or specify environment with anndata: configure_python_env(conda_env = \"1cellbio\")")
        ))
      }
      
      # If verbose, get more details but with error protection
      if (verbose) {
        tryCatch({
          # Try to get Python info but don't let it hang
          py_info <- reticulate::py_discover_config()
          if (!is.null(py_info$python)) {
            cat(icb_i18n("Python配置:\n", "Python configuration:\n"))
            cat("  - Python:", py_info$python, "\n")
            if (!is.null(py_info$version)) {
              cat("  - Version:", py_info$version, "\n")
            }
          }
          
          # Try to get anndata version
          anndata <- reticulate::import("anndata", delay_load = TRUE)
          anndata_version <- anndata$`__version__`
          cat("  - anndata version:", anndata_version, "\n")
        }, error = function(e) {
          if (verbose) cat(icb_i18n(
            "  - 详细信息获取失败，但基本配置成功\n",
            "  - Failed to get details, but basic configuration successful\n"
          ))
        })
      }
      
      TRUE
    }, timeout = 10)  # Shorter timeout for validation
  }, TimeoutException = function(e) {
    stop(icb_i18n(
      "Python环境验证超时。环境可能配置不正确。",
      "Python environment validation timed out. Environment may be misconfigured."
    ))
  }, error = function(e) {
    stop(icb_i18n(
      paste0("Python环境验证失败: ", e$message),
      paste0("Python environment validation failed: ", e$message)
    ))
  })
  
  if (verbose) cat("✓ ", icb_i18n("Python环境配置成功\n", "Python environment configured successfully\n"))
  return(invisible(TRUE))
}

#' Check if Python environment is properly configured for anndata
#'
#' @param verbose Logical, whether to print status messages (default: TRUE)
#' @return Logical, TRUE if anndata is available
#'
#' @export
check_anndata_available <- function(verbose = TRUE) {
  
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    if (verbose) cat(icb_i18n("reticulate包不可用\n", "reticulate package not available\n"))
    return(FALSE)
  }
  
  # Check if Python is configured
  if (!reticulate::py_available()) {
    if (verbose) cat(icb_i18n("Python不可用\n", "Python not available\n"))
    return(FALSE)
  }
  
  # Check if anndata module is available
  anndata_available <- reticulate::py_module_available("anndata")
  
  if (verbose) {
    if (anndata_available) {
      tryCatch({
        anndata <- reticulate::import("anndata")
        anndata_version <- anndata$`__version__`
        cat("✓ ", icb_i18n("anndata可用 (版本:", "anndata available (version:"), anndata_version, ")\n")
      }, error = function(e) {
        cat("✓ ", icb_i18n("找到anndata模块但版本检查失败\n", "anndata module found but version check failed\n"))
      })
    } else {
      cat("✗ ", icb_i18n("anndata模块不可用\n", "anndata module not available\n"))
    }
  }
  
  return(anndata_available)
}

#' Check if Python environment is configured
#' 
#' @param verbose Logical, whether to print status messages (default: TRUE)
#' @return Logical, TRUE if Python environment is properly configured
#' @keywords internal
is_python_configured <- function(verbose = TRUE) {
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    return(FALSE)
  }
  
  # Check if Python is available and configured
  py_available <- tryCatch({
    reticulate::py_available(initialize = FALSE)
  }, error = function(e) FALSE)
  
  if (!py_available) {
    return(FALSE)
  }
  
  # Check if anndata is available
  return(check_anndata_available(verbose = verbose))
}

#' Automatically detect conda executable path
#'
#' @return Character path to conda executable, or NULL if not found
#' @keywords internal
detect_conda_path <- function() {
  # First check PATH for conda
  conda_from_path <- Sys.which("conda")
  if (conda_from_path != "" && file.exists(conda_from_path)) {
    return(conda_from_path)
  }

  # Check common conda installation locations
  common_paths <- c(
    "~/anaconda3/bin/conda",
    "~/miniconda3/bin/conda",
    "/usr/local/bin/conda",
    "/opt/anaconda3/bin/conda",
    "/opt/miniconda3/bin/conda",
    "~/.conda/bin/conda"
  )

  for (path in common_paths) {
    expanded <- path.expand(path)
    if (file.exists(expanded)) {
      return(expanded)
    }
  }

  return(NULL)
}

#' List available conda environments with anndata
#'
#' @param verbose Logical, whether to print status messages (default: TRUE)
#' @return Data frame of conda environments with anndata availability
#' @keywords internal
list_conda_envs_with_anndata <- function(verbose = TRUE) {
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    if (verbose) cat(icb_i18n("reticulate包不可用\n", "reticulate package not available\n"))
    return(NULL)
  }

  # Auto-detect and set conda path for reticulate
  conda_path <- detect_conda_path()
  if (!is.null(conda_path)) {
    Sys.setenv(RETICULATE_CONDA = conda_path)
    if (verbose) cat(icb_i18n(
      paste0("检测到conda: ", conda_path, "\n"),
      paste0("Detected conda: ", conda_path, "\n")
    ))
  }

  # Check if conda is available
  conda_available <- tryCatch({
    conda_envs <- reticulate::conda_list()
    !is.null(conda_envs) && nrow(conda_envs) > 0
  }, error = function(e) FALSE)
  
  if (!conda_available) {
    if (verbose) {
      cat(icb_i18n(
        "⚠️  未检测到conda环境。\n请安装Anaconda或Miniconda：https://docs.conda.io/en/latest/miniconda.html\n",
        "⚠️  No conda environments detected.\nPlease install Anaconda or Miniconda: https://docs.conda.io/en/latest/miniconda.html\n"
      ))
    }
    return(NULL)
  }
  
  # Get conda environments
  conda_envs <- reticulate::conda_list()
  
  if (verbose) {
    cat(icb_i18n("正在检查conda环境中的anndata可用性...\n", "Checking anndata availability in conda environments...\n"))
  }
  
  # Check each environment for anndata
  env_status <- data.frame(
    name = conda_envs$name,
    python = conda_envs$python,
    has_anndata = FALSE,
    anndata_version = NA_character_,
    stringsAsFactors = FALSE
  )
  
  for (i in seq_len(nrow(conda_envs))) {
    env_name <- conda_envs$name[i]
    python_path <- conda_envs$python[i]
    
    if (verbose) cat("  ", icb_i18n("检查", "Checking"), env_name, "...")
    
    # Check if anndata is available in this environment
    # Create a temporary Python script
    temp_script <- tempfile(fileext = ".py")
    writeLines("import anndata; print(anndata.__version__)", temp_script)
    
    # Run the script and capture result (suppress warnings)
    result <- tryCatch({
      suppressWarnings(system2(python_path, temp_script, stdout = TRUE, stderr = TRUE))
    }, error = function(e) {
      result <- character(0)  # Empty character vector
      attr(result, "status") <- 999  # Error status
      result
    })
    
    # Clean up
    unlink(temp_script)
    
    # Check if command succeeded
    status <- attr(result, "status")
    if (is.null(status)) status <- 0  # NULL status means success in R
    

    
    if (status == 0 && length(result) > 0) {
      version_str <- trimws(result[1])
      # Check if result looks like a version number
      if (nchar(version_str) > 0 && grepl("^[0-9]", version_str)) {
        env_status$has_anndata[i] <- TRUE
        env_status$anndata_version[i] <- version_str
        if (verbose) cat(" ✓ anndata", version_str, "\n")
      } else {
        if (verbose) cat(" ✗ (unexpected output:", version_str, ")\n")
      }
    } else {
      if (verbose) cat(" ✗ (status:", status, ")\n")
    }
  }
  
  return(env_status)
}

#' Smart Python environment auto-configuration
#' 
#' @param verbose Logical, whether to print status messages (default: TRUE)
#' @param interactive Logical, whether to prompt user for environment selection (default: TRUE)
#' @return Logical, TRUE if configuration successful
#' @export
smart_python_config <- function(verbose = TRUE, interactive = TRUE) {
  
  # Step 1: Check if already configured
  if (is_python_configured(verbose = FALSE)) {
    if (verbose) cat("✓ ", icb_i18n("Python环境已配置并可用\n", "Python environment already configured and available\n"))
    return(TRUE)
  }
  
  if (verbose) {
    cat(icb_i18n(
      "🔍 Python环境未配置，正在自动检测可用环境...\n",
      "🔍 Python environment not configured, auto-detecting available environments...\n"
    ))
  }
  
  # Step 2: List conda environments with anndata
  env_status <- list_conda_envs_with_anndata(verbose = verbose)
  
  if (is.null(env_status)) {
    # No conda available
    stop(icb_i18n(
      "❌ 未检测到conda。请安装Anaconda或Miniconda后重试。\n安装指南：https://docs.conda.io/en/latest/miniconda.html",
      "❌ Conda not detected. Please install Anaconda or Miniconda and try again.\nInstallation guide: https://docs.conda.io/en/latest/miniconda.html"
    ))
  }
  
  # Step 3: Find environments with anndata
  envs_with_anndata <- env_status[env_status$has_anndata, ]
  
  if (nrow(envs_with_anndata) == 0) {
    # No environments have anndata
    stop(icb_i18n(
      paste0("❌ 在所有conda环境中都未找到anndata。\n",
            "请在现有环境中安装anndata：\n",
            "  conda activate your_env\n",
            "  pip install anndata\n",
            "或创建新环境：\n",
            "  conda create -n 1cellbio python=3.9\n",
            "  conda activate 1cellbio\n",
            "  pip install anndata pandas numpy"),
      paste0("❌ No anndata found in any conda environment.\n",
            "Please install anndata in an existing environment:\n",
            "  conda activate your_env\n",
            "  pip install anndata\n",
            "Or create a new environment:\n",
            "  conda create -n 1cellbio python=3.9\n",
            "  conda activate 1cellbio\n",
            "  pip install anndata pandas numpy")
    ))
  }
  
  # Step 4: Auto-select or prompt user
  if (nrow(envs_with_anndata) == 1) {
    # Only one environment with anndata, use it automatically
    selected_env <- envs_with_anndata$name[1]
    if (verbose) {
      cat(icb_i18n(
        paste0("✓ 自动选择唯一可用环境: ", selected_env, " (anndata ", envs_with_anndata$anndata_version[1], ")\n"),
        paste0("✓ Auto-selecting only available environment: ", selected_env, " (anndata ", envs_with_anndata$anndata_version[1], ")\n")
      ))
    }
  } else {
    # Multiple environments available
    if (verbose) {
      cat(icb_i18n("📋 发现多个包含anndata的环境:\n", "📋 Found multiple environments with anndata:\n"))
      for (i in seq_len(nrow(envs_with_anndata))) {
        cat(sprintf("  %d. %s (anndata %s)\n", 
                   i, envs_with_anndata$name[i], envs_with_anndata$anndata_version[i]))
      }
    }
    
    if (interactive && interactive()) {
      # Interactive mode: prompt user
      cat(icb_i18n("请选择要使用的环境 (1 -", "Please select environment to use (1 -"), nrow(envs_with_anndata), "): ")
      choice <- as.integer(readline())
      
      if (is.na(choice) || choice < 1 || choice > nrow(envs_with_anndata)) {
        stop(icb_i18n("无效选择", "Invalid selection"))
      }
      
      selected_env <- envs_with_anndata$name[choice]
    } else {
      # Non-interactive mode: use first environment
      selected_env <- envs_with_anndata$name[1]
      if (verbose) {
        cat(icb_i18n(
          paste0("自动选择第一个环境: ", selected_env, "\n"),
          paste0("Auto-selecting first environment: ", selected_env, "\n")
        ))
      }
    }
  }
  
  # Step 5: Configure the selected environment
  if (verbose) {
    cat(icb_i18n(
      paste0("🔧 正在配置环境: ", selected_env, "\n"),
      paste0("🔧 Configuring environment: ", selected_env, "\n")
    ))
  }
  
  configure_python_env(conda_env = selected_env, verbose = verbose)
  
  return(TRUE)
}