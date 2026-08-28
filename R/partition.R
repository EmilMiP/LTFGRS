#' Generate row indices for partitioning a matrix-like object
#'
#' @param nrow          Integer. Total number of rows in the object to be
#'                      partitioned.
#' @param n_partitions  Integer. Number of (roughly) equal partitions.
#'                      Mutually exclusive with `partition_size`.
#' @param partition_size Integer. Target number of rows per partition.
#'                       Mutually exclusive with `n_partitions`.
#'
#' @return A named list of integer vectors, each containing the row indices
#'         for one partition.
#'
#' @details
#' This function returns indices only and does not hold any data in memory.
#' Intended for use in loops or iterators where partitions are materialised
#' one at a time. The remainder from unequal division is retained as a smaller
#' final partition rather than distributed or discarded.
#'
#' When \code{n_partitions} is used, partition sizes differ by at most one row,
#' with larger partitions appearing first. When \code{partition_size} is used,
#' all partitions except possibly the last will have exactly
#' \code{partition_size} rows.
#'
#' @seealso \code{partition_rows()} for a version that returns subsetted
#'          objects directly.
#' @export
#'
#' @examples
#' idx <- partition_indices(nrow = 20, n_partitions = 3)
#'
partition_indices <- function(nrow, n_partitions = NULL, partition_size = NULL) {

  # --- Input validation -------------------------------------------------------

  if (!is.numeric(nrow) || length(nrow) != 1L ||
      nrow < 1L || nrow != floor(nrow)) {
    stop("`nrow` must be a single positive integer.")
  }

  nrow <- as.integer(nrow)

  if (is.null(n_partitions) == is.null(partition_size)) {
    stop("Provide exactly one of `n_partitions` or `partition_size`.")
  }

  # --- Compute groups ---------------------------------------------------------

  if (!is.null(n_partitions)) {

    if (!is.numeric(n_partitions) || length(n_partitions) != 1L ||
        n_partitions < 1L || n_partitions != floor(n_partitions)) {
      stop("`n_partitions` must be a single positive integer.")
    }

    n_partitions <- min(as.integer(n_partitions), nrow)
    groups       <- cut(seq_len(nrow), breaks = n_partitions, labels = FALSE)

  } else {

    if (!is.numeric(partition_size) || length(partition_size) != 1L ||
        partition_size < 1L || partition_size != floor(partition_size)) {
      stop("`partition_size` must be a single positive integer.")
    }

    partition_size <- as.integer(partition_size)
    groups         <- ceiling(seq_len(nrow) / partition_size)

  }

  # --- Build index list -------------------------------------------------------

  n_parts   <- max(groups)
  pad_width <- nchar(as.character(n_parts))
  nms       <- paste0("partition_", formatC(seq_len(n_parts),
                                            width = pad_width,
                                            flag  = "0"))

  indices <- lapply(seq_len(n_parts), function(i) which(groups == i))
  stats::setNames(indices, nms)
}
