#' Default Distance Bands
#'
#' Calculates distance quantiles from a reproducible sample of support points.
#' This is the default bandwidth calculation used by [seg_measure()]. Sampling uses
#' a fixed local seed and population weights when supplied. The caller's
#' random-number state is preserved. When all points fit within `n`, every
#' point is used. Distance quantiles are unweighted, and coincident-point
#' distances are omitted.
#'
#' @param x A numeric two-column coordinate matrix, a data frame with columns
#'   selected by `coords`, or an [sf::sf] object. Polygon inputs use the
#'   representative locations selected by `polygon_points`.
#' @param coords Two coordinate column names or positions for a data frame.
#'   Omit for sf input or a two-column coordinate matrix.
#' @param population Optional population-weight column name or numeric vector.
#'   Defaults to one per support point. For aggregated variables, supply the
#'   relevant count totals. At least one weight must be positive.
#' @param n Maximum number of support points to sample; an integer of at least 2.
#' @param probs Requested distance quantiles, in [0, 1]. The defaults are
#'   `c(0.01, 0.025, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5)`.
#' @param polygon_points Representative polygon locations: `"point_on_surface"`
#'   (default) or `"centroid"`. See [seg_measure()].
#' @return A numeric vector of distance bands, in coordinate units.
#' @seealso [seg_measure()]
#' @examples
#' points <- data.frame(x = 1:10, y = 0, residents = 1:10)
#' seg_default_bands(points, coords = c("x", "y"), population = "residents")
#' @export
seg_default_bands <- function(x, coords = NULL, population = NULL, n = 500,
                          probs = c(0.01, 0.025, 0.05, 0.1, 0.15, 0.2,
                                    0.25, 0.3, 0.4, 0.5),
                          polygon_points = c("point_on_surface", "centroid")) {
  if (!is.numeric(n) || length(n) != 1L || !is.finite(n) || n < 2 || n != floor(n))
    stop("'n' must be an integer of at least 2", call. = FALSE)
  if (!is.numeric(probs) || !length(probs) || any(!is.finite(probs)) ||
      any(probs < 0 | probs > 1))
    stop("'probs' must contain probabilities in [0, 1]", call. = FALSE)
  checked <- .seg_geometry(x, NULL, coords, match.arg(polygon_points))
  weights <- .seg_mass(population, checked$data, nrow(checked$coords))
  .seg_distance_bands(checked$coords, weights, n, probs)
}

# Shared by seg_measure() and seg_default_bands() so both use the same band rule
# without revalidating coordinates and weights. Uses a reproducible sample
# when needed, then returns quantiles of the positive pairwise distances.
.seg_distance_bands <- function(
  coords, weights, n = 500,
  probs = c(0.01, 0.025, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5)
) {
  total <- sum(weights)
  if (!is.finite(total) || total <= 0)
    stop("Population total must be positive and finite", call. = FALSE)

  n_points <- nrow(coords)
  if (n_points > n) {
    sample_n <- min(n, sum(weights > 0))
    selected <- withr::with_seed(1, sample.int(n_points, sample_n, prob = weights))
    coords <- coords[selected, , drop = FALSE]
  }

  d <- as.numeric(dist(coords))
  if (any(!is.finite(d)))
    stop("Default bands require finite pairwise distances", call. = FALSE)
  d <- d[d > 0]
  if (length(d) == 0) {
    stop(
      "failed to calculate default bands from the input coordinates",
      call. = FALSE
    )
  }

  as.numeric(quantile(d, probs = probs, names = FALSE, type = 8))
}
