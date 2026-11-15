library(tidyverse)


seg_pattern <- function(
  N = 1000L,
  side = 1L,
  concentrate = 1,
  cycles = 1,
  center = 1,
  resource.gini = 0.5,
  resource.mean = 10,
  status = c(low = 1, middle = 1, upper = 1),
  group = c(B = 1, A = 1),
  noise.agents = 1,
  noise.location = 1,
  noise.assign = 1,
  tol = 0.1
) {
  require(dplyr)
  require(tidyr)

  if (any(is.null(status)) || any(is.na(status))) {
    status <- 0L
  }
  if (any(is.null(group)) || any(is.na(group))) {
    group <- 0L
  }

  # Helpers
  sd <- sqrt(2) * qnorm((1 + resource.gini) / 2)
  stdz <- function(x) {
    if (sd(x, na.rm = TRUE) == 0) {
      return(x)
    }
    (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
  }
  fuzzyrank <- function(x, noise) {
    rank(stdz(x) + noise * rnorm(length(x)))
  }
  wgts_to_factor <- function(x, N) {
    if (length(x) == 1L) {
      return(rep(x, N))
    }
    y <- round(N * x[-length(x)] / sum(x))
    y <- c(y, N - sum(y))
    factor(rep(names(x), times = y), levels = names(x))
  }
  loc_function <- function(x, y, cycles, center, tol, log = FALSE) {
    stopifnot(
      is.numeric(tol),
      length(tol) == 1L,
      is.finite(tol),
      tol > 0
    )
    if (log) {
      stdz(log(
        1 +
          tol +
          cos(2 * pi * cycles * x / side) * cos(2 * pi * cycles * y / side)
      )) -
        center * stdz(log((x^2 + y^2)))
    } else {
      stdz(
        cos(2 * pi * cycles * x / side) *
          cos(2 * pi * cycles * y / side)
      ) -
        center * stdz((x^2 + y^2))
    }
  }

  # output structure
  out <- list(agents = NULL, locations = NULL)

  # generate the agents, their attributes, and location scores
  agents <- tibble(
    resource = rlnorm(
      N,
      log(resource.mean) - sd^2 / 2,
      sd
    )
  ) |>
    mutate(
      status = wgts_to_factor(status, N)[fuzzyrank(
        min_rank(resource),
        noise.agents
      )],
      group = wgts_to_factor(group, N)[fuzzyrank(
        stdz(resource) + stdz(as.integer(status)),
        noise.agents
      )],
      score = fuzzyrank(
        stdz(resource) + stdz(as.integer(status) + stdz(as.integer(group))),
        noise.agents
      )
    )
  out$agents <- agents

  # generate the locations and their agent assignment scores
  locations <- tibble(
    point_id = seq_len(2 * N),
    x = runif(2 * N, -side / 2, side / 2),
    y = runif(2 * N, -side / 2, side / 2)
  ) |>
    filter(
      (fuzzyrank(
        stdz(loc_function(x, y, cycles, center, tol)),
        noise.location
      ) >
        N)
    ) |>
    mutate(score = min_rank(-loc_function(x, y, cycles, center, tol)))

  out$locations <- locations

  out
}

result <- seg_pattern(
  N = 1e4,
  cycles = 3,
  center = 1,
  noise.location = 1,
  noise.assign = 0,
  noise.agents = 0
)
ineq::Gini(result$agents$resource)

result$agents |>
  summarize(mean(resource), sd(resource), .by = c(status, group)) |>
  arrange(`mean(resource)`)

result$agents |>
  ggplot(aes(
    x = score,
    group = interaction(group, status),
    fill = interaction(group, status)
  )) +
  geom_density(alpha = 0.5, color = NA) +
  scale_x_continuous()

result$locations |>
  ggplot(aes(x = x, y = y)) +
  geom_point(size = 0.1) +
  coord_equal()

left_join(result$locations, result$agents) |>
  ggplot(aes(x = x, y = y, color = min_rank(resource))) +
  geom_point(size = 0.2) +
  coord_equal() +
  scale_color_viridis_c()

left_join(result$locations, result$agents) |>
  ggplot(aes(x = x, y = y, color = group)) +
  geom_point(size = 0.2) +
  coord_equal()
