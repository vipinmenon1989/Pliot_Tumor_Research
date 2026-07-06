# Provenance and Reproducibility Tracking Utilities
# Captures file digests, package versions, and configuration state.

library(digest)
library(jsonlite)

calculate_file_checksum <- function(filepath) {
  if (!file.exists(filepath)) {
    stop(sprintf("File does not exist: %s", filepath))
  }
  digest(filepath, file = TRUE, algo = "md5")
}

# Logs and creates a provenance metadata record for a transformation step
record_provenance <- function(step_name, inputs, outputs, parameters, dataset = "global") {
  timestamp <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")

  # Calculate input checksums
  input_checksums <- list()
  for (name in names(inputs)) {
    path <- inputs[[name]]
    if (file.exists(path)) {
      input_checksums[[name]] <- list(path = path, checksum = calculate_file_checksum(path))
    } else {
      input_checksums[[name]] <- list(path = path, checksum = "NOT_FOUND_ON_DISK")
    }
  }

  # Calculate output checksums
  output_checksums <- list()
  for (name in names(outputs)) {
    path <- outputs[[name]]
    if (file.exists(path)) {
      output_checksums[[name]] <- list(path = path, checksum = calculate_file_checksum(path))
    } else {
      output_checksums[[name]] <- list(path = path, checksum = "NOT_YET_CREATED")
    }
  }

  # Package environments
  git_commit <- tryCatch({
    trimws(system("git rev-parse HEAD", intern = TRUE))
  }, error = function(e) "NO_GIT_COMMIT")

  git_status <- tryCatch({
    paste(system("git status --porcelain", intern = TRUE), collapse = "\n")
  }, error = function(e) "NO_GIT_STATUS")

  prov_record <- list(
    timestamp = timestamp,
    step = step_name,
    dataset = dataset,
    git_commit = git_commit,
    git_status = git_status,
    inputs = input_checksums,
    outputs = output_checksums,
    parameters = parameters,
    session_info = capture.output(sessionInfo())
  )

  return(prov_record)
}

save_provenance_json <- function(prov_record, output_path) {
  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
  writeLines(toJSON(prov_record, auto_unbox = TRUE, pretty = TRUE), output_path)
}
