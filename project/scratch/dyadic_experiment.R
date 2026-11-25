set.seed(2026)
N <- 2000L
side <- 1000
resource_mean <- 10
resource_gini <- 0.5
resource_sd <- sqrt(2) * qnorm((1 + resource_gini) / 2)
resource_noise <- 1

resource_sample <- sort(rlnorm(
  N,
  log(resource_mean) - resource_sd^2 / 2,
  resource_sd
))

stdz <- function(x) {
  (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
}

locations <- tibble(
  point_id = seq_len(N),
  x = runif(N, 0, side),
  y = runif(N, 0, side)
)

agents <- tibble(
  agent_id = seq_len(N),
  group = sample(rep(c("A", "B"), each = N / 2))
) |>
  mutate(
    status = stdz(1L * (group == "A")) + resource_noise * rnorm(N),
    status = rank(status) / n(),
    status = case_when(
      status < (1 / 3) ~ 1L,
      status < (2 / 3) ~ 2L,
      .default = 3L
    )
  ) |>
  mutate(
    resources = stdz(status + 1L * (group == "A")) + resource_noise * rnorm(N),
    resources = resource_sample[rank(resources, ties = "first")],
    resources_rank = rank(resources),
    status = factor(status, levels = c(1, 2, 3), ordered = TRUE),
    group = factor(group)
  )

location_noise <- rnorm(N)
assign_city <- function(cycles = 1, noise = 0.5) {
  periodic <- with(
    locations,
    cos(2 * pi * cycles * x / side) *
      cos(2 * pi * cycles * y / side)
  )
  utility <- log((1 + 1e-6 + periodic) / 2)
  score <- as.numeric(scale(utility)) + noise * location_noise
  destination <- order(score, locations$point_id)
  source <- order(agents$resources, agents$agent_id)
  assigned <- integer(N)
  assigned[destination] <- source
  bind_cols(locations, slice(agents, assigned))
}

scenarios <- crossing(cycles = c(1, 3), noise = c(0, 1)) |>
  mutate(pattern = paste0(cycles, " cycle(s); noise = ", noise))

patterns <- map2(scenarios$cycles, scenarios$noise, assign_city) |>
  set_names(scenarios$pattern) |>
  imap(\(city, label) mutate(city, pattern = label))


p <- 3
dyads <- patterns[[p]] |>
  select(point_id, x, y, group, resources) |>
  rename_with(~ paste0(.x, "_ref")) |>
  mutate(
    point_id_tar = lapply(
      point_id_ref,
      \(x) point_id_ref[point_id_ref > x]
    )
  ) |>
  unnest(point_id_tar) |>
  left_join(
    patterns[[p]] |>
      select(point_id, x, y, group, resources) |>
      rename_with(~ paste0(.x, "_tar"))
  ) |>
  mutate(dist = sqrt((x_ref - x_tar)^2 + (y_ref - y_tar)^2)) |>
  mutate(
    group_comp = paste0(
      pmin(as.character(group_ref), as.character(group_tar)),
      pmax(as.character(group_ref), as.character(group_tar))
    ),
    resource_diff = min_rank(resources_tar) - min_rank(resources_ref)
  )

dyads |>
  ggplot(aes(x = min_rank(dist), fill = group_comp, group = group_comp)) +
  geom_density(alpha = 0.2, color = NA)

dyads |>
  ggplot(aes(x = min_rank(dist), y = resource_diff)) +
  geom_hex()
