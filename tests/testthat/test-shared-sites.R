test_that("shared sites match explicit aggregation for all variable types", {
  xy <- expand.grid(x = c(50, 150, 250), y = c(50, 150))
  rows <- rep(seq_len(nrow(xy)), each = 6)
  d <- xy[rows, ]
  d$g <- factor(rep(c("A", "B", "C", "A", "B", "A"), 6))
  d$o <- ordered(rep(c(1, 3, 2, 1, 2, 3), 6))
  d$v <- rep(c(10, 10, 20, 30, 50, 80), 6)
  d$w <- rep(c(.5, 1, 2, 1.5, 0, 1), 6)
  d$g[c(2, 9)] <- NA
  d$o[5] <- NA
  d$v[14] <- NA
  counts <- xy
  for (column in c("g", "o", "v")) {
    values <- sort(unique(d[[column]][!is.na(d[[column]])]))
    for (value in as.character(values)) {
      mass <- ifelse(!is.na(d[[column]]) & as.character(d[[column]]) == value, d$w, 0)
      counts[[value]] <- as.numeric(rowsum(mass, rows, reorder = FALSE))
    }
  }
  micro <- list(g = seg_categorical("g", weights = "w", na.action = "omit"),
    o = seg_ordinal("o", weights = "w", na.action = "omit"),
    v = seg_continuous("v", weights = "w", na.action = "omit",
                       summary = list(probs = c(.25, .5, .75))))
  aggregate <- list(g = seg_categorical(c("A", "B", "C")),
    o = seg_ordinal(as.character(1:3)),
    v = seg_continuous(as.character(c(10, 20, 30, 50, 80)),
                       representation = "exact", scores = c(10, 20, 30, 50, 80),
                       summary = list(probs = c(.25, .5, .75))))
  for (mode in c("radius", "knn")) for (search in c("brute", "boost", "grid")) {
    for (kernel in c("unweighted", "biweight", "inverse", "exponential")) {
      run <- function(x, variables, population) seg_measure(x, variables,
        coords = c("x", "y"), population = population, neighbors = mode,
        search = search, weighting = kernel,
        bands = if (mode == "radius") c(0, 100, 250) else c(1, 8, 15),
        output = "full", keep_scores = TRUE, keep_contributions = TRUE)
      actual <- run(d, micro, d$w)
      expected <- run(counts, aggregate, as.numeric(rowsum(d$w, rows)))
      for (what in c("indices", "pairwise", "exposure", "env", "scores",
                     "contributions", "diagnostics")) {
        expect_equal(as.data.frame(actual, what = what),
                     as.data.frame(expected, what = what), tolerance = 1e-11)
      }
      expect_equal(actual$diagnostics$focal_searches, nrow(xy))
      expect_equal(actual$site, rows)
      expect_equal(actual$data, d)
      expect_equal(actual$coords, as.matrix(xy), ignore_attr = TRUE)
    }
  }
})

test_that("site results are invariant to record order and retain exact quantiles", {
  d <- data.frame(x = c(0, 0, 1, 1, 2, 2), y = 0,
                  v = c(1, 100, 2, 80, 3, 60))
  variables <- list(v = seg_continuous("v",
    summary = list(type = "quantiles", probs = c(.25, .5, .75))))
  run <- function(x) seg_measure(x, variables, coords = c("x", "y"),
    neighbors = "knn", bands = c(1, 3), weighting = "unweighted", output = "full")
  r <- run(d)
  shuffled <- run(d[c(2, 1, 4, 3, 6, 5), ])
  expect_equal(r$analyses, shuffled$analyses)
  expect_equal(r$site, c(1, 1, 2, 2, 3, 3))
  expect_equal(nrow(sf::st_as_sf(r)), 3 * 2 * 3)
  expect_equal(r$analyses$v$knn$min_population, c(1, 3))
  expect_equal(r$analyses$v$knn$max_population, c(1, 3))
})

test_that("repeated geometries group by shape rather than representative point", {
  square <- function(size) sf::st_polygon(list(rbind(c(-size, -size),
    c(size, -size), c(size, size), c(-size, size), c(-size, -size))))
  d <- sf::st_sf(a = c(1, 2, 3), b = c(3, 2, 1),
                 geometry = sf::st_sfc(square(1), square(2), square(1)))
  r <- seg_measure(d, list(g = seg_categorical(c("a", "b"))),
    polygon_points = "centroid", bands = 0, output = "full")
  expect_equal(r$diagnostics$focal_searches, 2)
  expect_equal(r$site, c(1, 2, 1))
  expect_equal(length(r$geometry), 2)
  expect_equal(unname(r$analyses$g$env[[1]]), matrix(.5, 2, 2))
})

test_that("default bands use unique support with summed population weights", {
  xy <- data.frame(x = 1:600, y = 0)
  d <- xy[rep(seq_len(600), each = 2), ]
  d$g <- rep(c("A", "B"), 600)
  expected <- seg_default_bands(xy, coords = c("x", "y"), population = rep(2, 600))
  expect_equal(seg_default_bands(d, coords = c("x", "y")), expected)
  r <- seg_measure(d, list(g = seg_categorical("g")), coords = c("x", "y"))
  expect_equal(r$bands, expected)
  expect_null(r$site)
  expect_null(r$data)
  expect_equal(r$diagnostics$focal_searches, 600)
})

test_that("cases share sites despite different numbers of colocated records", {
  first <- data.frame(x = c(0, 0, 1, 2), y = 0, year = 1,
    g = factor(c("A", "B", "A", "B")), v = c(1, 5, 2, 9))
  second <- data.frame(x = c(0, 1, 1, 1, 2), y = 0, year = 2,
    g = factor(c("B", "A", "B", "B", "A")), v = c(8, 1, 3, 6, 4))
  specs <- list(g = seg_categorical("g"), v = seg_continuous("v"))
  for (mode in c("radius", "knn")) {
    run <- function(x, ...) seg_measure(x, specs, coords = c("x", "y"),
      neighbors = mode, bands = c(1, 2), output = "full", ...)
    result <- run(rbind(first, second), case = "year")
    expect_equal(result$diagnostics$focal_searches, 3)
    expect_equal(result$cases[[1]]$analyses, run(first)$analyses)
    expect_equal(result$cases[[2]]$analyses, run(second)$analyses)
    expect_equal(result$cases[[1]]$site, c(1, 1, 2, 3))
    expect_equal(result$cases[[2]]$site, c(1, 2, 2, 2, 3))
    expect_setequal(unique(sf::st_as_sf(result)$case), 1:2)
  }
})

test_that("high-cardinality discrete sites use sparse membership correctly", {
  site <- c(1:6, 1, 2)
  d <- data.frame(x = site, y = 0,
    g = factor(c(letters[1:6], "b", "f")))
  counts <- data.frame(x = 1:6, y = 0)
  for (category in letters[1:6])
    counts[[category]] <- as.numeric(rowsum(as.numeric(d$g == category), site))
  for (type in list(seg_categorical, seg_ordinal)) {
    if (identical(type, seg_ordinal)) d$g <- ordered(d$g)
    run <- function(x, spec) seg_measure(x, list(g = spec), coords = c("x", "y"),
      bands = c(0, 3), output = "full", keep_contributions = TRUE)
    actual <- run(d, type("g"))
    expected <- run(counts, type(letters[1:6]))
    for (what in c("indices", "pairwise", "exposure", "env", "contributions"))
      expect_equal(as.data.frame(actual, what = what), as.data.frame(expected, what = what))
  }
})

test_that("repeated polygons conserve counts before surface construction", {
  boxes <- sf::st_make_grid(sf::st_as_sfc(sf::st_bbox(
    c(xmin = 0, ymin = 0, xmax = 2, ymax = 1))), n = c(2, 1))
  d <- sf::st_sf(a = c(2, 3, 5), b = c(4, 6, 1),
                 geometry = boxes[c(1, 2, 1)])
  combined <- sf::st_sf(a = c(7, 3), b = c(5, 6), geometry = boxes)
  variables <- list(g = seg_categorical(c("a", "b")))
  for (surface in c("grid", "pycno")) {
    run <- function(x) seg_measure(x, variables, bands = c(0, 1),
      surface = surface, nrow = 2, ncol = 4, output = "full", surface_audit = TRUE)
    result <- run(d)
    expect_equal(result$analyses, run(combined)$analyses)
    expect_equal(colSums(result$data), c(a = 10, b = 11))
    expect_equal(result$surface_info$source_records$site, c(1, 2, 1))
    expect_equal(result$surface_info$source_records$data, sf::st_drop_geometry(d))
    expect_equal(result$site, seq_len(8))
  }
})
