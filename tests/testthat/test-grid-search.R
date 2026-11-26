test_that("grid offsets reproduce radius and population-kNN analyses", {
  xy <- expand.grid(x = seq(0, 6, 2), y = seq(0, 6, 3))[-c(3, 8), ]
  xy <- rbind(xy, xy[c(1, 5), ])
  n <- nrow(xy)
  d <- transform(xy, a = seq_len(n) %% 3, b = seq_len(n) %% 4,
                 c = seq_len(n) %% 2, v = seq_len(n) %% 5)
  d[1, c("a", "b", "c")] <- 0
  variables <- list(
    group = seg_categorical(c("a", "b", "c")),
    ordinal = seg_ordinal(c("a", "b", "c")),
    interval = seg_continuous(c("a", "b", "c")),
    exact_counts = seg_continuous(c("a", "b", "c"), scores = c(1, 5, 10)),
    exact_rows = seg_continuous("v")
  )
  for (mode in c("radius", "knn")) {
    for (kernel in c("unweighted", "biweight", "inverse", "exponential")) {
      for (normalize in c(FALSE, TRUE)) {
        run <- function(search) seg_measure(d, variables, coords = c("x", "y"),
          population = rep(1, n), neighbors = mode,
          bands = if (mode == "radius") c(6, 0, 2, 3, 6) else c(5, 1, 9, 5),
          weighting = kernel, normalize = normalize, search = search,
          output = "full", keep_scores = TRUE, keep_contributions = TRUE)
        grid <- run("grid")
        expect_equal(grid$analyses, run("boost")$analyses, tolerance = 1e-11)
        expect_equal(grid$analyses, run("brute")$analyses, tolerance = 1e-11)
        expect_equal(grid$diagnostics$focal_searches, n)
        expect_equal(grid$diagnostics$grid$cell_size, c(2, 3))
        expect_equal(grid$diagnostics$grid$dimensions, c(4L, 3L))
      }
    }
  }
})

test_that("grid layouts handle absent cells, coincident rows, and single axes", {
  cases <- list(
    data.frame(x = c(0, 0, 0), y = 0, g = c("a", "b", "b")),
    data.frame(x = c(0, 0.25, 0.5) + 1000, y = 2e6, g = c("a", "b", "b")),
    data.frame(x = 0, y = c(0, 1, 3), g = c("a", "b", "b"))
  )
  for (d in cases) for (mode in c("radius", "knn")) {
    run <- function(search) seg_measure(d, list(g = seg_categorical("g")),
      coords = c("x", "y"), bands = if (mode == "radius") c(0, 1, 3) else c(1, 2),
      neighbors = mode, search = search, output = "full")
    expect_equal(run("grid")$analyses, run("boost")$analyses, tolerance = 1e-12)
  }
  d <- data.frame(x = c(0, 1e6), y = 0, g = c("a", "b"))
  r <- seg_measure(d, "g", coords = c("x", "y"), bands = 1,
                   search = "grid", grid_size = 1, output = "full")
  expect_identical(r$diagnostics$grid$lookup, "sparse")
  expect_equal(unname(r$analyses$g$env[[1]]), diag(2))
  d <- expand.grid(x = seq(0, 10, by = 0.1), y = c(0, 0.2))
  d$g <- rep(c("a", "b"), length.out = nrow(d))
  run <- function(search) seg_measure(d, list(g = seg_categorical("g")),
    coords = c("x", "y"), bands = c(0, 0.251, 1.001),
    search = search, output = "full")
  inferred <- run("grid")
  expect_equal(inferred$diagnostics$grid$cell_size, c(0.1, 0.2))
  expect_equal(inferred$analyses, run("boost")$analyses, tolerance = 1e-11)
})

test_that("constructed surfaces provide grid spacing to the offset backend", {
  boxes <- sf::st_make_grid(sf::st_as_sfc(sf::st_bbox(
    c(xmin = 0, ymin = 0, xmax = 6, ymax = 4))), n = c(2, 1))
  d <- sf::st_sf(a = c(10, 2), b = c(0, 8), geometry = boxes)
  for (surface in c("grid", "pycno")) {
    run <- function(search) seg_measure(d, list(g = seg_categorical(c("a", "b"))),
      bands = c(0, 2, 4), surface = surface, nrow = 4, ncol = 3,
      search = search, output = "full")
    grid <- run("grid")
    expect_equal(grid$diagnostics$grid$cell_size, c(2, 1))
    expect_equal(grid$analyses, run("boost")$analyses, tolerance = 1e-11)
    expect_identical(run("auto")$neighbors$engine, "grid")
  }
})

test_that("grid search rejects misalignment and excessive caches explicitly", {
  d <- data.frame(x = c(0, 1, 2.2), y = 0, g = c("a", "b", "b"))
  run <- function(...) seg_measure(d, "g", coords = c("x", "y"), bands = 1, ...)
  expect_error(run(search = "grid"), "rectangular lattice")
  expect_error(run(search = "brute", grid_size = 1), "requires search")
  for (size in list(0, -1, NA_real_, Inf, "one", 1:3, numeric()))
    expect_error(run(search = "grid", grid_size = size), "positive finite")
  d$x <- c(0, 1, 1e9)
  expect_error(seg_measure(d, "g", coords = c("x", "y"), bands = 1e9,
    search = "grid", grid_size = 1), "cache budget")
  d$x <- c(0, .Machine$double.eps, 1)
  expect_error(run(search = "grid", grid_size = 1), "Distinct coordinates")
})

test_that("automatic search chooses the neighborhood-specific backend", {
  irregular <- data.frame(x = c(0, .2, 1.3, 1.7, 3.1),
    y = c(0, .3, 1.1, .4, 2.7), g = c("a", "a", "b", "b", "a"))
  regular <- expand.grid(x = 0:3, y = 0:3)
  regular$g <- ifelse(regular$x < 2, "a", "b")
  for (d in list(irregular, regular)) for (mode in c("radius", "knn")) {
    expected <- if (identical(d, regular)) "grid" else
      if (mode == "radius") "brute" else "boost"
    run <- function(...) seg_measure(d, list(g = seg_categorical("g")),
      coords = c("x", "y"), neighbors = mode,
      bands = if (mode == "radius") c(0, 1.5) else c(1, 3),
      output = "full", ...)
    automatic <- run()
    expect_identical(automatic$neighbors$engine, expected)
    expect_equal(automatic$analyses, run(search = expected)$analyses)
    for (search in c("boost", "brute"))
      expect_identical(run(search = search)$neighbors$engine, search)
  }
})

test_that("automatic grids handle sparse spacing and coincident rows deliberately", {
  d <- data.frame(x = c(0, 1, 1e6), y = 0, g = c("a", "b", "b"))
  run <- function(...) seg_measure(d, list(g = seg_categorical("g")),
    coords = c("x", "y"), bands = 1, ...)
  expect_identical(run()$neighbors$engine, "brute")
  expect_identical(run(grid_size = 1)$neighbors$engine, "grid")
  d$x <- c(0, .Machine$double.eps, 1)
  expect_identical(run()$neighbors$engine, "brute")
  expect_error(run(grid_size = 1), "Distinct coordinates")
})
