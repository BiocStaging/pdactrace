sb_lock_manifest <- function(paths, root = ".") {
  if (!requireNamespace("digest", quietly = TRUE)) {
    stop("Install digest to create a prospective lock manifest.", call. = FALSE)
  }
  root <- normalizePath(root, winslash = "/", mustWork = TRUE)
  paths <- sort(unique(normalizePath(paths, winslash = "/", mustWork = TRUE)))
  prefix <- paste0(root, .Platform$file.sep)
  if (any(!startsWith(paths, prefix))) {
    stop("Every locked file must be inside root.", call. = FALSE)
  }
  relative <- substring(paths, nchar(prefix) + 1L)
  data.frame(
    path = relative,
    bytes = unname(file.info(paths)$size),
    sha256 = vapply(paths, function(path) {
      digest::digest(path, algo = "sha256", serialize = FALSE, file = TRUE)
    }, character(1)),
    stringsAsFactors = FALSE
  )
}

sb_verify_lock_manifest <- function(manifest, root = ".") {
  if (is.character(manifest) && length(manifest) == 1L) {
    manifest <- utils::read.csv(manifest, check.names = FALSE, stringsAsFactors = FALSE)
  }
  required <- c("path", "bytes", "sha256")
  if (!is.data.frame(manifest) || !all(required %in% names(manifest))) {
    stop("manifest must contain path, bytes and sha256 columns.", call. = FALSE)
  }
  root <- normalizePath(root, winslash = "/", mustWork = TRUE)
  paths <- file.path(root, manifest$path)
  exists <- file.exists(paths)
  current_bytes <- rep(NA_real_, length(paths))
  current_sha256 <- rep(NA_character_, length(paths))
  current_bytes[exists] <- unname(file.info(paths[exists])$size)
  if (any(exists)) {
    if (!requireNamespace("digest", quietly = TRUE)) {
      stop("Install digest to verify a prospective lock manifest.", call. = FALSE)
    }
    current_sha256[exists] <- vapply(paths[exists], function(path) {
      digest::digest(path, algo = "sha256", serialize = FALSE, file = TRUE)
    }, character(1))
  }
  status <- ifelse(
    !exists,
    "MISSING",
    ifelse(
      current_bytes != manifest$bytes,
      "SIZE_MISMATCH",
      ifelse(current_sha256 != manifest$sha256, "HASH_MISMATCH", "PASS")
    )
  )
  data.frame(
    path = manifest$path,
    status = status,
    expected_sha256 = manifest$sha256,
    observed_sha256 = current_sha256,
    stringsAsFactors = FALSE
  )
}
