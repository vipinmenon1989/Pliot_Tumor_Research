# Structured Logging Utilities for MPNST Phase 1 Analysis
# Enforces standard format, records runtime, memory, and environment.

log_info <- function(message, dataset = "global", stage = "general") {
  timestamp <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
  cat(sprintf("[%s] [INFO] [Dataset: %s] [Stage: %s] %s\n", timestamp, dataset, stage, message))
}

log_warn <- function(message, dataset = "global", stage = "general") {
  timestamp <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
  cat(sprintf("[%s] [WARNING] [Dataset: %s] [Stage: %s] %s\n", timestamp, dataset, stage, message), file = stderr())
}

log_error <- function(message, dataset = "global", stage = "general") {
  timestamp <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
  cat(sprintf("[%s] [ERROR] [Dataset: %s] [Stage: %s] %s\n", timestamp, dataset, stage, message), file = stderr())
}

# Records system resources and session details
log_system_usage <- function(dataset = "global", stage = "general") {
  mem_info <- tryCatch({
    # Attempt to read memory usage from system if on Linux
    gc_res <- gc()
    peak_mem_mb <- sum(gc_res[, 6]) # Max Mb used
    sprintf("Peak GC memory: %.2f MB", peak_mem_mb)
  }, error = function(e) "Memory info unavailable")

  log_info(sprintf("System usage status - %s", mem_info), dataset, stage)
}

# Enforce a global warning handler that logs but does not suppress
setup_strict_logging <- function() {
  options(warn = 1) # Print warnings as they occur
}
