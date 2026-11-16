test_that("radius searches include closed boundaries and coincident points", {
  d <- data.frame(x=c(0,9,0,10),y=c(0,4,0,0),a=c(2,0,0,0),b=c(0,2,2,2))
  vars <- list(g=seg_categorical(vars=c("a","b")))
  for (band in c(0,sqrt(97))) {
    indexed <- seg_measure(d,vars,coords=c("x","y"),bands=band,
                     weighting="unweighted",output="full")
    brute <- seg_measure(d,vars,coords=c("x","y"),bands=band,
                   weighting="unweighted",output="full",search="brute")
    expect_equal(indexed$analyses,brute$analyses,tolerance=1e-12)
    expect_equal(unname(indexed$analyses$g$env[[1]][1,]),
                 if(band==0)c(.5,.5) else c(1/3,2/3))
  }
})

test_that("KNN includes tied boundary populations proportionately", {
  d <- data.frame(x=c(0,-1,1,0),y=c(0,0,0,1),a=c(10,20,0,0),b=c(0,0,10,30))
  vars <- list(g=seg_categorical(vars=c("a","b")))
  run <- function(d,search="boost") seg_measure(d,vars,coords=c("x","y"),
    neighbors="knn",bands=25,weighting="unweighted",output="full",search=search)
  indexed <- run(d)
  # Ten focal people, then 15/60 of each population at distance one.
  expect_equal(unname(indexed$analyses$g$env[[1]][1,]),c(.6,.4))
  expect_equal(indexed$analyses,run(d,"brute")$analyses,tolerance=1e-12)
  permutation <- c(1,4,2,3)
  reordered <- run(d[permutation,])
  expect_equal(reordered$analyses$g$env[[1]],
               indexed$analyses$g$env[[1]][permutation,],tolerance=1e-12)
  expect_equal(reordered$analyses$g$indices,indexed$analyses$g$indices,tolerance=1e-12)
})

test_that("coincident observations share a fractional boundary after the focal", {
  d <- data.frame(x=0,y=0,g=c("a","b","b","b"))
  r <- seg_measure(d,list(g=seg_categorical(vars="g")),coords=c("x","y"),
             neighbors="knn",bands=2,weighting="unweighted",output="full")
  expect_equal(unname(r$analyses$g$env[[1]][1,]),c(.5,.5))
  expect_equal(unname(r$analyses$g$env[[1]][2,]),c(1/6,5/6))
  expect_equal(r$analyses$g$knn$min_population,2)
  expect_equal(r$analyses$g$knn$max_population,2)
})
test_that("radius band reuse preserves independent results and requested order", {
  d <- data.frame(x = c(0, 0, 0.4, 1, 2), y = 0,
                  a = c(2, 0, 2, 0, 1), b = c(0, 2, 0, 3, 1))
  vars <- list(g = seg_categorical(c("a", "b")))
  bands <- c(2, 0, 0.4, 1, 0.4)
  for (weighting in c("unweighted", "biweight", "inverse", "exponential")) {
    for (normalize in c(FALSE, TRUE)) {
      run <- function(bands) seg_measure(
        d, vars, coords = c("x", "y"), bands = bands,
        weighting = weighting, normalize = normalize, output = "full"
      )
      joint <- run(bands)
      expect_identical(joint$bands, bands)
      expect_equal(joint$diagnostics$focal_searches, nrow(d))
      for (b in seq_along(bands)) {
        separate <- run(bands[b])$analyses$g
        expect_identical(joint$analyses$g$env[[b]], separate$env[[1]])
        expect_equal(joint$analyses$g$indices$overall$h[b],
                     separate$indices$overall$h)
      }
    }
  }
})
test_that("running counts and sparse exact ranks coexist across radius bands", {
  d <- data.frame(x = c(0, 0, .25, .5, 1), y = 0, value = c(1, 2, 2, 4, 8),
                  a = c(1, 0, 1, 3, 2), b = c(0, 1, 2, 1, 3))
  vars <- list(
    ordinal = seg_ordinal(c("a", "b")),
    interval = seg_continuous(c("a", "b")),
    exact_counts = seg_continuous(c("a", "b"), scores = c(1, 5)),
    exact_rows = seg_continuous("value")
  )
  run <- function(bands) seg_measure(
    d, vars, coords = c("x", "y"), population = d$a + d$b,
    bands = bands, weighting = "unweighted", output = "full",
    keep_scores = TRUE, keep_contributions = TRUE
  )
  bands <- c(1, 0, .25)
  joint <- run(bands)
  for (band in bands) {
    separate <- run(band)
    for (what in c("indices", "pairwise", "exposure", "env", "scores",
                   "contributions")) {
      expect_equal(as.data.frame(joint, what = what, bands = band),
                   as.data.frame(separate, what = what))
    }
  }
})
