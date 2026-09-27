#' Assign exact non-dominated Pareto fronts
#'
#' @param components Numeric matrix; larger values are preferred.
#' @return Integer front labels, where 1 is non-dominated.
#' @export
sb_pareto_fronts <- function(components) {
  x <- as.matrix(components)
  storage.mode(x) <- "double"
  if (!nrow(x)) return(integer())
  x[!is.finite(x)] <- 0.5
  remaining <- seq_len(nrow(x))
  fronts <- integer(nrow(x))
  front <- 1L
  while (length(remaining)) {
    dominated <- vapply(remaining, function(i) {
      candidates <- x[remaining, , drop = FALSE]
      all_ge <- rowSums(sweep(candidates, 2, x[i, ], `>=`)) == ncol(x)
      any_gt <- rowSums(sweep(candidates, 2, x[i, ], `>`)) > 0
      any(all_ge & any_gt)
    }, logical(1))
    current <- remaining[!dominated]
    if (!length(current)) {
      current <- remaining[[which.max(rowSums(x[remaining, , drop = FALSE]))]]
    }
    fronts[current] <- front
    remaining <- setdiff(remaining, current)
    front <- front + 1L
  }
  fronts
}

