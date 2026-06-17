# 1. CHECKPOINT SYSTEM ####
# Enables crash recovery: save progress -> reload on restart
# Pattern: check_checkpoint() -> if not found, compute -> save_checkpoint()


#' Save a named R object as a checkpoint .rds file (atomic write)
#' Uses temp-file + rename to prevent corrupt files on mid-save crash
save_checkpoint <- function(obj, name, base = output_path) {
  dir <- file.path(base, "RData")
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  tmp <- paste0(dir, "/",CHECKPOINT_PREFIX, name, "_TMP.rds")
  final <- paste0(dir, "/", CHECKPOINT_PREFIX, name, ".rds")
  saveRDS(obj, file = tmp)
  file.rename(tmp, final)
  message(sprintf(
    "[CHECKPOINT] Saved: %s  (%s)",
    final, format(Sys.time(), "%H:%M:%S")
  ))
}

#' Check whether a checkpoint file exists
check_checkpoint <- function(name, base = output_path) {
  file.exists(paste0(base,"/","RData/", CHECKPOINT_PREFIX, name, ".rds"))
}

#' Load a checkpoint file and return the object
load_checkpoint <- function(name, base = output_path) {
  dir <- file.path(base, "RData")
  path <- paste0(dir, "/", CHECKPOINT_PREFIX, name, ".rds")
  if (!file.exists(path)) stop("[CHECKPOINT] Not found: ", path)
  message(sprintf(
    "[CHECKPOINT] Loaded: %s  (%s)",
    path, format(Sys.time(), "%H:%M:%S")
  ))
  readRDS(path)
}

# 2. Plot Helpers ####

#' Save a ggplot to TIFF and PDF; increment global image_number counter
#' @param plotname  Base name of the file (no extension)
#' @param plot      A ggplot object
#' @param width     Width in inches
#' @param height    Height in inches
#' @param subdir    Subdirectory under figures (e.g. "QC", "DE")
save_plot <- function(plotname, plot, width = 8, height = 6, subdir = "") {
  dir <- file.path(output_path, "figures", subdir)
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  filename <- paste0(dir,"/",sprintf("%03d", image_number), "_", plotname)
  tryCatch(
    {
      if (inherits(plot, c("ggplot", "patchwork", "Heatmap"))) {
        ggsave(paste0(filename, tiff_extension), plot, width = width, height = height, units = "in", dpi = 300, bg = "white")
        ggsave(paste0(filename, pdf_extension), plot, width = width, height = height, units = "in", bg = "white")
      } else if (is.function(plot)) {
        # for base-graphics functions that draw directly to the device (e.g. DimHeatmap)
        tiff(paste0(filename, tiff_extension), width = width, height = height, units = "in", res = 300)
        plot()
        dev.off()
        pdf(paste0(filename, pdf_extension), width = width, height = height)
        plot()
        dev.off()
      } else {
        warning(sprintf("[SAVE_PLOT] '%s' — unsupported plot type '%s', skipping.", plotname, class(plot)[1]))
        return(invisible(NULL))
      }
      image_number <<- image_number + 1
    },
    error = function(e) {
      warning(sprintf("[SAVE_PLOT] Failed to save '%s': %s", filename, e$message))
    }
  )
}



# Free RAM #####

#' Get current RAM usage in MB (rounded)
#' Note: This is a simple approximation based on gc() output and may not reflect total RAM usage of the R session.
#' For more accurate monitoring, consider using the 'pryr' package or system tools.
#' @returns RAM usage in MB (numeric)
#' @examples
ram_mb <- function() round(sum(gc(verbose = FALSE)[, 2]), 0)

# Logging #####

#' Initialize a timestamped log file and return logging utilities
#' @param script_name  Label used in the log filename (e.g. "05_sct_normalization")
#' @param log_dir      Directory for log files (default: "logs")
#' @returns A list with `log_msg(...)` and `close_log()`.
#'   In the calling chunk: `on.exit(logger$close_log(), add = TRUE)`
setup_logging <- function(script_name, log_dir = "logs") {
  dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
  log_file <- file.path(
    log_dir,
    paste0(script_name, "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log")
  )
  log_con <- file(log_file, open = "wt")

  log_msg <- function(...) {
    line <- paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", paste0(...))
    message(line)
    writeLines(line, log_con)
  }

  close_log <- function() {
    if (tryCatch(isOpen(log_con), error = function(e) FALSE)) {
      log_msg("[EXIT] Script finished.")
      close(log_con)
    }
  }

  list(log_msg = log_msg, close_log = close_log, log_file = log_file)
}
# Usage example:
#   logger  <- setup_logging("05_sct_normalization")
#log_msg <- logger$log_msg
#on.exit(logger$close_log(), add = TRUE)

#log_msg("Starting SCT normalization — RAM: ", ram_mb(), " MB")


# Save Table #####
#' Save a data frame as a TSV file with a timestamped filename
#' @param df        Data frame to save
#' @param name      Base name of the file (no extension)
#' @param subdir    Subdirectory under tables/ (e.g. "QC/",
#' "DE/")
#' @returns The path of the saved file (invisibly)
#' @examples
#' save_table(my_data, "my_results", subdir = "DE/")
#' This function saves a data frame as a TSV file in the specified subdirectory under "tables/". The filename includes a timestamp to ensure uniqueness. The function returns the path of the saved file invisibly.
#' Note: The output directory is determined by the global variable `output_path`, which should be defined in the main script.
#' The function creates the target directory if it does not exist and handles any errors that may occur during the file writing process, providing informative messages to the user.
#' Example usage:
#' ```R
#' # Assuming `my_data` is a data frame you want to save
#' save_table(my_data, "my_results", subdir = "DE/")
#' ```
#' 

save_table <- function(df, name, subdir = "") {
  dir <- file.path(output_path, "tables", subdir)
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  filename <- paste0(dir, timestamp, "_", name, ".tsv")
  tryCatch(
    {
      write.table(df, file = filename, sep = "\t", row.names = FALSE, quote = FALSE)
      message(sprintf("[SAVE_TABLE] Saved: %s", filename))
      invisible(filename)
    },
    error = function(e) {
      warning(sprintf("[SAVE_TABLE] Failed to save '%s': %s", filename, e$message))
      invisible(NULL)
    }
  )
}
