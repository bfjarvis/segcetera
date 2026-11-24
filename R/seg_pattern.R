#' Generate a synthetic segregation pattern
#'
#' Simulates agents with lognormal resources and associated status and group
#' labels, along with an equal number of candidate residential locations.
#' Agents and locations are matched one-to-one by rank score.
#'
#' @param N Positive integer number of agents and retained locations.
#' @param side Positive side length of the square, centered at the origin.
#' @param cycles Nonnegative number of cosine cycles across each side.
#' @param center Nonnegative weight favoring locations nearer the center.
#' @param resource.gini Target Gini coefficient of the resource distribution,
#'   in `[0, 1)`.
#' @param resource.mean Positive mean of the resource distribution.
#' @param status,group Named vectors of nonnegative relative category weights,
#'   in assignment order (lower to higher score). At least one weight must
#'   be positive. Use `NULL` to omit a category assignment.
#' @param noise.agents Nonnegative standard deviation of noise added to the
#'   scores used to assign status, group, and the final agent rank.
#' @param noise.location Nonnegative standard deviation of noise added when
#'   selecting which candidate locations to retain.
#' @param coef.score Numeric weights for the standardized attributes in the
#'   final agent score. Supply any subset of the names `resource`, `status`,
#'   and `group`. Omitted coefficients retain defaults of 1, 0, and 0,
#'   respectively. Unnamed entries are ignored.
#' @param center positive number indicating how much locations and groups
#'   should be centralized. Overrides separate agent and location centers.
#' @param center.locations positive number indicating how much locations
#'   should be centralized.
#' @param center.agents positive number indicating how much particular groups
#'   should be centralized.
#' @param pow positive number indicating how deep the valleys
#'   are in the periodic structure for agents and locations. Overrides
#'   `pow.locations` and `pow.agents` arguments.
#' @param pow.locations positive number indicating how deep the valleys
#'   are in the periodic structure for creating locations.
#' @param pow.agents positive number indicating how deep the valleys
#'   are in the periodic structure when assigning agents to locations.
#' @param coef.score Numeric weights for the standardized attributes in the
#'   final agent score. Supply any subset of the names `resource`, `status`,
#'   and `group`. Omitted coefficients retain defaults of 1, 0, and 0,
#'   respectively. Unnamed entries are ignored.
#' @details The mean and Gini are distribution parameters; realized sample
#'   values vary. Category sizes are rounded to integers by assigning
#'   remaining agents to categories with the largest fractional remainders.
#'   All random draws use the caller's RNG stream; use [set.seed()] to
#'   reproduce a draw. Of `2 * N` uniformly sampled candidate locations, the
#'   `N` with highest noisy location utility are retained. Their final
#'   scores rank utility without the selection noise.
#' @return A data frame with one row per retained location and its assigned
#'   agent: `point_id`, `x`, `y`, `score`, `resource`, `status` (an ordered
#'   factor), and `group` (a factor). Omitted category assignments contain zero.
#'   Score 1 pairs the highest-utility location with the lowest-ranked agent.
#' @examples
#' set.seed(42)
#' pattern <- seg_pattern(N = 100, cycles = 2, noise.agents = 0.2)
#' plot(pattern$x, pattern$y, col = as.integer(pattern$group), asp = 1)
#' @export
seg_pattern <- function(
  N = 1000L,
  side = 1L,
  cycles = 1,
  resource.gini = 0.5,
  resource.mean = 10,
  status = c(low = 1, middle = 1, upper = 1),
  group = c(B = 1, A = 1),
  noise.agents = 1,
  noise.location = 1,
  center = 1,
  center.locations = NULL,
  center.agents = NULL,
  pow = 2,
  pow.locations = NULL,
  pow.agents = NULL,
  coef.score = NULL,
  grid = c("square", "hex", "runif")
) {
  coefficients <- c(resource = 1, status = 0, group = 0)
  grid <- match.arg(grid)
  if (!is.null(coef.score)) {
    if (!is.numeric(coef.score) || !is.null(dim(coef.score))) {
      stop("'coef.score' must be a numeric vector", call. = FALSE)
    }
    labels <- names(coef.score)
    if (!is.null(labels)) {
      if (anyNA(labels)) {
        stop("'coef.score' names must not be NA", call. = FALSE)
      }
      named <- nzchar(labels)
      keys <- labels[named]
      unknown <- setdiff(keys, names(coefficients))
      if (length(unknown)) {
        stop(
          "Unknown 'coef.score' name(s): ",
          paste(unknown, collapse = ", "),
          ". Use resource, status, or group.",
          call. = FALSE
        )
      }
      if (anyDuplicated(keys)) {
        stop("'coef.score' names must be unique", call. = FALSE)
      }
      if (any(!is.finite(coef.score[named]))) {
        stop("Named 'coef.score' values must be finite", call. = FALSE)
      }
      coefficients[keys] <- coef.score[named]
    }
  }
  if (is.null(pow.locations)) {
    pow.locations <- pow
  }
  if (is.null(pow.agents)) {
    pow.agents <- pow
  }
  if (is.null(center.locations)) {
    center.locations <- center
  }
  if (is.null(center.agents)) {
    center.agents <- center
  }

  scalars <- list(
    N = N,
    side = side,
    cycles = cycles,
    center = center,
    resource.gini = resource.gini,
    resource.mean = resource.mean,
    noise.agents = noise.agents,
    noise.location = noise.location,
    pow.locations = pow.locations,
    pow.agents = pow.agents
  )
  for (name in names(scalars)) {
    value <- scalars[[name]]
    if (
      !is.numeric(value) ||
        length(value) != 1L ||
        !is.finite(value) ||
        value < 0
    ) {
      stop("'", name, "' must be a finite nonnegative number", call. = FALSE)
    }
  }
  if (N < 1 || N != floor(N) || N > .Machine$integer.max %/% 2L) {
    stop(
      "'N' must be a positive integer with room for 2 * N locations",
      call. = FALSE
    )
  }
  if (side <= 0 || resource.mean <= 0 || resource.gini >= 1) {
    stop(
      "'side' and 'resource.mean' must be positive; 'resource.gini' must be below 1",
      call. = FALSE
    )
  }

  category_pool <- function(weights, name, ordered = FALSE) {
    if (is.null(weights)) {
      return(rep(0L, N))
    }
    labels <- names(weights)
    if (
      !is.numeric(weights) ||
        is.null(labels) ||
        anyNA(labels) ||
        any(!nzchar(labels)) ||
        anyDuplicated(labels) ||
        any(!is.finite(weights) | weights < 0) ||
        !any(weights > 0)
    ) {
      stop(
        "'",
        name,
        "' must have unique names and nonnegative weights with positive total",
        call. = FALSE
      )
    }
    weights <- weights / max(weights)
    expected <- N * (weights / sum(weights))
    counts <- floor(expected)
    remainder <- N - sum(counts)
    extra <- order(expected - counts, decreasing = TRUE)[seq_len(remainder)]
    counts[extra] <- counts[extra] + 1L
    factor(rep(labels, counts), levels = labels, ordered = ordered)
  }
  status_pool <- category_pool(status, "status", ordered = TRUE)
  group_pool <- category_pool(group, "group")

  standardize <- function(x) {
    spread <- if (length(x) > 1L) stats::sd(x) else 0
    if (spread == 0) {
      return(x)
    }
    (x - mean(x)) / spread
  }
  fuzzy_rank <- function(x, noise) {
    rank(
      standardize(x) + noise * stats::rnorm(length(x)),
      ties.method = "first"
    )
  }
  location_utility <- function(x, y, center, pow) {
    x <- x / side
    y <- y / side
    periodic <- cos(2 * pi * (cycles * x)) * cos(2 * pi * (cycles * y))
    # -((x^2 + y^2)^center) * periodic^pow
    -center * (standardize(x^2 + y^2)) + standardize(periodic^pow)
  }

  log_sd <- sqrt(2) * stats::qnorm(resource.gini / 2 + 0.5)
  if (!is.finite(log_sd)) {
    stop(
      "'resource.gini' is too close to 1 for numerical precision",
      call. = FALSE
    )
  }
  resource <- stats::rlnorm(N, log(resource.mean) - log_sd^2 / 2, log_sd)
  status <- status_pool[fuzzy_rank(
    rank(resource, ties.method = "min"),
    noise.agents
  )]
  group <- group_pool[fuzzy_rank(
    standardize(resource) + standardize(as.integer(status)),
    noise.agents
  )]
  score <- fuzzy_rank(
    coefficients[["resource"]] *
      standardize(resource) +
      coefficients[["status"]] * standardize(as.integer(status)) +
      coefficients[["group"]] * standardize(as.integer(group)),
    noise.agents
  )
  if (grid == 'runif') {
    locations <- data.frame(
      x = stats::runif(2 * N, -side / 2, side / 2),
      y = stats::runif(2 * N, -side / 2, side / 2)
    )
  }
  if (grid %in% c('hex', 'square')) {
    sfc <- st_sfc(st_polygon(list(rbind(
      c(-side / 2, -side / 2),
      c(side / 2, -side / 2),
      c(side / 2, side / 2),
      c(-side / 2, -side / 2)
    ))))
    if (grid == 'hex') {
      locations <- st_coordinates(st_make_grid(
        sfc,
        what = "centers",
        n = ceiling(1.4 * rep(sqrt(N), 2)),
        square = FALSE
      ))
    } else if (grid == "square") {
      locations <- st_coordinates(st_make_grid(
        sfc,
        what = "centers",
        n = ceiling(1.5 * rep(sqrt(N), 2))
      ))
    }
    locations <- as.data.frame(locations)
    names(locations) <- c("x", "y")
  }
  locations$point_id <- seq_len(nrow(locations))

  utility <- location_utility(
    locations$x,
    locations$y,
    center.locations,
    pow.locations
  )
  locations <- locations[
    fuzzy_rank(utility, noise.location) > (nrow(locations) - N),
    ,
    drop = FALSE
  ]
  locations$score <- rank(
    -location_utility(locations$x, locations$y, center.agents, pow.agents),
    ties.method = "first"
  )
  rownames(locations) <- NULL
  agent <- match(locations$score, score)
  locations$resource <- resource[agent]
  locations$status <- status[agent]
  locations$group <- group[agent]
  locations
}
