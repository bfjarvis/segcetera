test_that("count preparation is reused during population inference", {
  d <- data.frame(x = 1:3, y = 0, a = c(1, NA, 3), b = c(3, 2, 1),
    c = c(2, 1, 2), e = c(2, 3, 2))
  original <- segcetera:::.seg_counts
  calls <- 0L
  testthat::local_mocked_bindings(.seg_counts = function(spec, data) {
    calls <<- calls + 1L
    original(spec, data)
  }, .package = "segcetera")
  # Omitted count rows still determine each variable's eligible population.
  specs <- list(first = seg_categorical(vars = c("a", "b"), na.action = "omit"),
    second = seg_categorical(vars = c("c", "e")))
  expect_error(seg_measure(d, specs, coords = c("x", "y"), bands = 0),
    "populations differ")
  expect_identical(calls, 2L)
  calls <- 0L
  result <- seg_measure(d, specs, coords = c("x", "y"), population = rep(4, 3),
    bands = 0, output = "full")
  expect_identical(calls, 2L)
  expect_equal(result$analyses$first$definition$eligible_population, 8)
  expect_true(all(is.na(result$analyses$first$env[[1]][2, ])))
})

test_that("case preparation resolves each supplied weight vector once", {
  d <- data.frame(x = rep(0:2, 2), y = 0, year = rep(1:2, each = 3),
                  g = rep(c("A", "B", "A"), 2), w = 1:6)
  original <- segcetera:::.seg_mass
  calls <- 0L
  testthat::local_mocked_bindings(.seg_mass = function(...) {
    calls <<- calls + 1L
    original(...)
  }, .package = "segcetera")
  result <- seg_measure(d, list(g = seg_categorical("g", weights = "w")),
    coords = c("x", "y"), population = d$w, case = "year", bands = 1)
  expect_identical(calls, 2L)
  expect_equal(vapply(result$cases, function(z)
    z$analyses$g$definition$eligible_population, numeric(1)), c(`1` = 6, `2` = 15))
  for (invalid in list(c(1, NA, 1, 1, 1, 1), rep(-1, 6), rep(Inf, 6), 1:2)) {
    expect_error(seg_measure(d,
      list(g = seg_categorical("g", weights = invalid)),
      coords = c("x", "y"), case = "year", bands = 1), "Population weights")
  }
})

test_that("effective measures retain aliases and variable overrides", {
  d <- data.frame(x = 0:2, y = 0, g = c("A", "B", "A"))
  run <- function(selected, global = "all") seg_measure(d,
    list(g = seg_categorical("g", measures = selected)),
    coords = c("x", "y"), bands = 1, measures = global)
  expect_equal(run("relative_diversity")$analyses, run("diversity")$analyses)
  expect_equal(run("variance_ratio")$analyses, run("diversity")$analyses)
  expect_equal(run("information", "unused")$analyses, run("information")$analyses)
  for (invalid in list(NULL, character(), NA_character_, "unknown")) {
    expect_error(run(invalid, "unknown"), "Unknown or empty measure")
  }
})

test_that("surface setup validates the current option vocabulary", {
  geom <- sf::st_make_grid(sf::st_as_sfc(sf::st_bbox(
    c(xmin = 0, ymin = 0, xmax = 2, ymax = 1))), n = c(2, 1))
  x <- sf::st_sf(a = c(4, 1), b = c(1, 4), geometry = geom)
  run <- function(...) seg_measure(x, list(race = seg_categorical(vars = c("a", "b"))),
    bands = 0, surface = "pycno", ...)
  expect_error(run(keep_geometry = TRUE), "surface options")
  expect_error(run(surface_geometry = "polygon"), "surface_geometry")
  expect_error(run(nrow = 2), "both 'nrow' and 'ncol'")
  expect_error(run(celldim = Inf), "positive finite")
  expect_error(run(celldim = NA_real_), "positive finite")
  expect_error(run(nrow = 2, ncol = 4, max_iter = 0), "positive integer")
  expect_error(run(nrow = 2, ncol = 4, converge = NA_real_), "finite numeric")
})
