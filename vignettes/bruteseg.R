local_shares <- function(city, bands, power = 3) {
  groups <- unique(city$group)

  crossing(
    city,
    city |> select(x_tar = x, y_tar = y, group_tar = group)
  ) |>
    mutate(d = sqrt((x - x_tar)^2 + (y - y_tar)^2)) |>
    group_by(x, y) |>
    reframe(
      bind_rows(lapply(bands, function(bw) {
        w <- pmax(1 - (d / bw)^power, 0)^power

        tibble(
          band = bw,
          group = groups,
          share = vapply(
            groups,
            function(g) {
              sum(w * (group_tar == g)) / sum(w)
            },
            numeric(1)
          )
        )
      }))
    )
}


bands <- c(10, 50, seq(100, 1000, by = 100))


system.time(shares <- local_shares(city, bands = bands))

system.time(
  result <- seg_measure(
    city,
    variables = "group",
    coords = c("x", "y"),
    bands = bands,
    neighbors = "radius",
    measures = "all",
    weighting = "biweight",
    scope = "both",
    output = "localenv"
  )
)
