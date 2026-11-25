test_that("synthetic patterns have reproducible one-to-one scores", {
  args <- list(N = 11, side = 2, status = c(low = 1, high = 2),
               group = c(B = 1, A = 0))
  pattern <- withr::with_seed(42, do.call(seg_pattern, args))
  expect_identical(pattern, withr::with_seed(42, do.call(seg_pattern, args)))
  expect_s3_class(pattern, "data.frame")
  expect_named(pattern, c("point_id", "x", "y", "score", "resource", "status", "group"))
  expect_true(is.ordered(pattern$status))
  expect_identical(as.integer(table(pattern$status)), c(4L, 7L))
  expect_identical(as.integer(table(pattern$group)), c(11L, 0L))
  expect_identical(sort(pattern$score), seq_len(11))
  expect_equal(anyDuplicated(pattern$point_id), 0)
  expect_true(all(abs(pattern$x) <= 1 & abs(pattern$y) <= 1))
})

test_that("constant utilities and omitted categories have well-defined matches", {
  for (n in c(1L, 10L)) {
    pattern <- withr::with_seed(42, seg_pattern(
      N = n, cycles = 0, center = 0, resource.gini = 0,
      status = NULL, group = NULL, noise.agents = 0, noise.location = 0
    ))
    expect_identical(pattern$status, rep(0L, n))
    expect_identical(pattern$group, rep(0L, n))
    expect_equal(pattern$resource, rep(10, n))
    expect_identical(pattern$score, seq_len(n))
  }
  expect_error(seg_pattern(N = 0), "'N'")
  expect_error(seg_pattern(group = c(A = 0, B = 0)), "positive total")
  expect_error(seg_pattern(status = c(1, 2)), "unique names")
})

test_that("equal scores pair low-resource agents with high-utility locations", {
  pattern <- withr::with_seed(42, seg_pattern(
    N = 10, cycles = 0, center = 1, status = NULL, group = NULL,
    noise.agents = 0, noise.location = 0
  ))
  expect_identical(order(pattern$resource), order(pattern$score))
  expect_identical(order(pattern$x^2 + pattern$y^2), order(pattern$score))
})
test_that("score coefficients override named defaults only", {
  run <- function(...) withr::with_seed(42, seg_pattern(N = 40, noise.agents = 0, ...))
  expect_identical(run(), run(coef.score = c(group = 0, resource = 1, status = 0)))
  expect_identical(run(), run(coef.score = c(9, 8)))
  expect_identical(run(coef.score = c(status = 2, 9)),
                   run(coef.score = c(resource = 1, status = 2, group = 0)))
  pattern <- run(coef.score = c(resource = -1))
  expect_identical(order(pattern$resource, decreasing = TRUE), order(pattern$score))
  expect_error(run(coef.score = c(resources = 1)), "Unknown.*resources")
  expect_error(run(coef.score = c(status = 1, status = 0)), "unique")
  expect_error(run(coef.score = c(group = Inf)), "finite")
  expect_error(run(coef.score = c(group = "1")), "numeric vector")
  expect_error(run(coef.score = setNames(1, NA_character_)), "must not be NA")
})

test_that("location patterns preserve their documented output and validate centers", {
  columns <- c("point_id", "x", "y", "score", "resource", "status", "group")
  for (grid in c("square", "hex", "runif")) {
    pattern <- withr::with_seed(42, seg_pattern(N = 20, grid = grid))
    expect_named(pattern, columns)
    expect_equal(nrow(pattern), 20)
    expect_true(all(abs(pattern$x) <= 0.5 & abs(pattern$y) <= 0.5))
  }
  expect_error(seg_pattern(center.locations = -1), "center.locations")
  expect_error(seg_pattern(center.agents = NA_real_), "center.agents")
})
