#' Download the two large upstream phase CSVs into a local cache
#'
#' The bundled `inst/extdata` carries the six small downstream phase
#' tables (phase2c / phase29 / phase42 / phase60 / phase77 / phase80)
#' but **not** the two large upstream fits — `phase33_deseq2_coef_12template.csv`
#' (~2.6 MB RNA fit) and `phase34_protein_pooled_12template.csv`
#' (~0.9 MB protein fit) — to keep the software tarball small. Both
#' are gene-level summary statistics and live in `data-raw/` of the
#' public source repository, `github.com/jibeomko/pdactrace`. This
#' helper fetches them on demand via raw GitHub URLs and caches them
#' locally with [BiocFileCache::BiocFileCache()].
#'
#' A git checkout does not need this helper: `data-raw/build_reference.R`
#' and `data-raw/build_protein_betas.R` read both files from `data-raw/`
#' directly.
#'
#' @param target Optional character.  Either `"phase33"`, `"phase34"`,
#'   or `"both"` (default).  Selects which CSV(s) to fetch.
#' @param ref Git ref (branch / tag / commit) on
#'   `jibeomko/pdactrace` to pull from.  Default `"main"`.
#' @param cache A [BiocFileCache::BiocFileCache] object.  Default
#'   uses the user's standard cache.
#' @param verbose Logical.  If `TRUE` (default), prints a progress
#'   message per file.
#' @return Invisibly, a named character vector mapping the requested
#'   `target` value(s) to the absolute cached file path(s).
#' @examples
#' args(download_phase_csvs)
#'
#' \donttest{
#'   # Pull both upstream CSVs once; subsequent calls hit the cache:
#'   paths <- download_phase_csvs("both")
#'   # ... and feed them to the build chain:
#'   p33 <- data.table::fread(paths[["phase33"]])
#' }
#' @export
download_phase_csvs <- function(target  = c("both", "phase33", "phase34"),
                                  ref     = "main",
                                  cache   = NULL,
                                  verbose = TRUE) {
  target <- match.arg(target)
  for (pkg in c("BiocFileCache", "tools")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(sprintf("`%s` is required for download_phase_csvs(). ", pkg),
           "Install via BiocManager::install(\"", pkg, "\").",
           call. = FALSE)
    }
  }
  if (is.null(cache)) {
    cache <- BiocFileCache::BiocFileCache(ask = FALSE)
  }

  base_url <- paste0(
    "https://raw.githubusercontent.com/jibeomko/pdactrace/",
    ref,
    "/data-raw/")

  files <- list(
    phase33 = "phase33_deseq2_coef_12template.csv",
    phase34 = "phase34_protein_pooled_12template.csv")
  wanted <- if (target == "both") names(files) else target

  out <- character(0L)
  for (key in wanted) {
    fname <- files[[key]]
    url   <- paste0(base_url, fname)
    rname <- paste0("pdactrace_", key, "_", ref)

    hit <- BiocFileCache::bfcquery(cache, rname, exact = TRUE)
    if (nrow(hit) > 0L) {
      fp <- BiocFileCache::bfcrpath(cache, rids = hit$rid[1L])
      if (isTRUE(verbose)) {
        message(sprintf("  cached  %s  (%s)",
                         fname,
                         format(structure(file.size(fp),
                                            class = "object_size"),
                                units = "auto")))
      }
    } else {
      if (isTRUE(verbose)) {
        message(sprintf("  fetch   %s ...", fname))
      }
      fp <- BiocFileCache::bfcadd(cache, rname = rname, fpath = url,
                                    download = TRUE)
    }
    names(fp) <- key
    out <- c(out, fp)
  }
  invisible(out)
}
