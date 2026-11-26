test_that("cases share searches while retaining independent distributions", {
  xy <- expand.grid(x = 0:2, y = 0:2)
  first <- transform(xy, year = 2020L, a = 1:9, b = 9:1,
                     c = rep(1:3, 3), v = rep(1:3, 3), w = 1)
  second <- transform(first, year = 2030L, a = b + 2, b = a %% 4,
                      c = c + 1, v = rev(v) + 1)
  d <- rbind(first, second)
  variables <- list(
    race = seg_categorical(c("a", "b", "c")),
    status = seg_ordinal(c("a", "b", "c")),
    intervals = seg_continuous(c("a", "b", "c")),
    income = seg_continuous("v", weights = "w", summary = list(probs = c(.25, .75)))
  )
  for (search in c("brute", "boost", "grid")) {
    for (mode in c("radius", "knn")) {
      for (kernel in c("unweighted", "biweight", "inverse", "exponential")) {
        run <- function(x, ...) suppressWarnings(seg_measure(x, variables,
          coords = c("x", "y"), search = search, neighbors = mode,
          bands = if (mode == "radius") c(2, 0, 1, 2) else c(1, 4),
          weighting = kernel, population = x$a + x$b + x$c,
          output = "full", keep_scores = TRUE, keep_contributions = TRUE, ...))
        batch <- run(d, case = "year")
        expect_equal(batch$case_ids, c(2020L, 2030L))
        expect_equal(batch$diagnostics$focal_searches, 9)
        expect_equal(batch$cases[[1]]$analyses, run(first)$analyses)
        expect_equal(batch$cases[[2]]$analyses, run(second)$analyses)
      }
    }
  }
})

test_that("case tables identify every kind of retained result", {
  xy <- expand.grid(x = 0:2, y = 0:2)
  d <- transform(xy, year = as.Date("2020-01-01"), a = 1:9, b = 9:1)
  d <- rbind(d, transform(d, year = as.Date("2021-01-01"), a = b))
  variables <- list(group = seg_categorical(c("a", "b")),
    income = seg_continuous(c("a", "b"), summary = list(probs = .5)))
  r <- suppressWarnings(seg_measure(d, variables, coords = c("x", "y"),
    case = "year", bands = c(0, 2), population = d$a + d$b, output = "full",
    keep_scores = TRUE, keep_contributions = TRUE))
  for (what in c("indices", "pairwise", "exposure", "env", "diagnostics",
                 "interpolation", "scores", "contributions")) {
    table <- as.data.frame(r, what = what)
    expect_s3_class(table$case, "Date")
    expect_setequal(unique(table$case), unique(d$year))
  }
  spatial <- sf::st_as_sf(r, variables = "group")
  expect_setequal(unique(spatial$case), unique(d$year))
  expect_equal(nrow(spatial), 2 * 9 * 2 * 2)
  expect_output(print(r), "2 cases")
  single <- seg_measure(d[1:9, ], variables[1], coords = c("x", "y"), bands = 1)
  expect_false("case" %in% names(as.data.frame(single)))
  compact <- seg_measure(d, variables[1], coords = c("x", "y"), case = "year", bands = 1)
  expect_null(compact$geometry)
  expect_true(all(vapply(compact$cases, function(z) is.null(z$data), logical(1))))
})

test_that("case alignment handles row order, missing support, and repeated points", {
  d <- data.frame(x = c(0, 0, 1, 1, 0), y = 0,
    year = c(1, 1, 1, 2, 2), g = factor(c("A", "B", "B", "A", "B")))
  variables <- list(g = seg_categorical("g", weights = c(1, 2, 3, 4, 5)))
  for (mode in c("radius", "knn")) {
    bands <- if (mode == "radius") c(0, 1) else c(1, 4)
    r <- seg_measure(d, variables, coords = c("x", "y"), case = "year",
                      neighbors = mode, bands = bands, output = "full")
    reference <- data.frame(x = c(0, 0, 1), y = 0,
                            g = factor(c("B", "A", "A"), levels = levels(d$g)))
    expected <- seg_measure(reference,
      list(g = seg_categorical("g", weights = c(5, 0, 4))),
      coords = c("x", "y"), neighbors = mode, bands = bands, output = "full")
    expect_equal(r$cases[[2]]$analyses, expected$analyses)
    expect_equal(r$diagnostics$focal_searches, 3)
  }
})

test_that("network cases reuse graph searches", {
  net <- seg_network(data.frame(from = c("a", "b", "c"), to = c("b", "c", "a")))
  first <- data.frame(node = c("a", "b", "c"), year = 1, a = c(2, 0, 3), b = c(0, 4, 1))
  second <- transform(first, year = 2, a = b + 1, b = a + 2)
  variables <- list(g = seg_categorical(c("a", "b")))
  for (mode in c("radius", "knn")) {
    run <- function(x, ...) seg_measure(x, variables, network = net,
      node_id = "node", bands = c(1, 2), neighbors = mode, output = "full", ...)
    r <- run(rbind(first, second), case = "year")
    expect_equal(r$diagnostics$focal_searches, 3)
    expect_equal(r$cases[[1]]$analyses, run(first)$analyses)
    expect_equal(r$cases[[2]]$analyses, run(second)$analyses)
    expect_setequal(as.data.frame(r, what = "env")$node, first$node)
  }
})

test_that("grid and pycno cases preserve their prepared populations", {
  boxes <- sf::st_make_grid(sf::st_as_sfc(sf::st_bbox(
    c(xmin = 0, ymin = 0, xmax = 4, ymax = 4))), n = c(2, 1))
  first <- sf::st_sf(year = 1, a = c(10, 3), b = c(2, 8), geometry = boxes)
  second <- transform(first, year = 2, a = b, b = a)
  second <- second[2:1, ]
  variables <- list(g = seg_categorical(c("a", "b")))
  for (surface in c("grid", "pycno")) {
    run <- function(x, ...) seg_measure(x, variables, bands = c(0, 2),
      surface = surface, nrow = 4, ncol = 4, output = "full", ...)
    r <- run(rbind(first, second), case = "year")
    expect_equal(r$cases[[1]]$analyses, run(first)$analyses)
    expect_equal(r$cases[[2]]$analyses, run(second[2:1, ])$analyses)
    expect_equal(r$diagnostics$focal_searches, 16)
    expect_s3_class(sf::st_as_sf(r), "sf")
    expect_null(r$cases[[1]]$surface_info$locations)
    audit <- run(rbind(first, second), case = "year", surface_audit = TRUE)
    expect_equal(audit$cases[[1]]$surface_info$locations, 1:16)
    expect_equal(audit$cases[[2]]$surface_info$source_rows, c(4L, 3L))
  }
})

test_that("case output handles empty pairwise tables and variable overrides", {
  d <- data.frame(x = c(0, 1, 0, 1), y = 0, year = c(1, 1, 2, 2),
                  g = c("A", "A", "B", "B"), v = c(1, 2, 2, 3))
  r <- seg_measure(d, list(g = seg_categorical("g")),
    coords = c("x", "y"), case = "year", bands = 1)
  empty <- as.data.frame(r, what = "pairwise")
  expect_equal(nrow(empty), 0)
  expect_true("case" %in% names(empty))
  r <- seg_measure(d, list(g = seg_categorical("g", output = "localenv"),
    v = seg_continuous("v", output = "indices")), coords = c("x", "y"),
    case = "year", bands = 1)
  expect_equal(unique(as.data.frame(r)$variable), "v")
  expect_equal(unique(as.data.frame(r, what = "env")$variable), "g")
  expect_error(as.data.frame(r, what = "env", variables = "v"), "not stored")
})

test_that("case selectors validate IDs and share default bands", {
  d <- data.frame(x = rep(0:3, 2), y = 0, year = rep(c(1, 2), each = 4),
                  g = rep(c("A", "B"), 4))
  variables <- list(g = seg_categorical("g"))
  run <- function(...) seg_measure(d, variables, coords = c("x", "y"), ...)
  expect_error(run(case = c("year", "g")), "one 'case'")
  expect_error(run(case = "unknown"), "Unknown column")
  d$year[1] <- NA
  expect_error(run(case = "year"), "Case IDs")
  d$year <- rep(c(1, 2), each = 4)
  r <- run(case = "year")
  expect_equal(r$bands, seg_default_bands(cbind(x = 0:3, y = 0), population = rep(2, 4)))
  expect_equal(r$cases[[1]]$bands, r$cases[[2]]$bands)
  one <- seg_measure(d[1:4, ], variables, coords = c("x", "y"),
                     case = "year", bands = 1)
  expect_equal(length(one$cases), 1)
  expect_equal(unique(as.data.frame(one)$case), 1)
})

test_that("absent polygon rows are zero-filled before surface construction", {
  boxes <- sf::st_make_grid(sf::st_as_sfc(sf::st_bbox(
    c(xmin = 0, ymin = 0, xmax = 4, ymax = 4))), n = c(2, 1))
  first <- sf::st_sf(year = 1, a = c(10, 3), b = c(2, 8), geometry = boxes)
  second <- first[1, ]
  second$year <- 2
  padded <- first
  padded$year <- 2
  padded$a[2] <- padded$b[2] <- 0
  for (surface in c("grid", "pycno")) for (mode in c("radius", "knn")) {
    run <- function(x, ...) seg_measure(x, list(g = seg_categorical(c("a", "b"))),
      bands = c(1, 2), neighbors = mode, surface = surface,
      nrow = 4, ncol = 4, output = "full", ...)
    r <- run(rbind(first, second), case = "year")
    expect_equal(r$cases[[2]]$analyses, run(padded)$analyses)
  }
})
