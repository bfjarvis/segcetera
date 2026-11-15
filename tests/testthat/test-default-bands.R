test_that("default bands share coordinate and population semantics with segcetera", {
  points <- data.frame(x=seq_len(600), y=rep(0:5,100),
                        residents=rep(c(0,1,2,3),150))
  probs <- c(.1,.5,.7)
  set.seed(901)
  before <- .Random.seed
  actual <- seg_default_bands(points,coords=c("x","y"),population="residents",
                           n=25,probs=probs)
  expect_identical(.Random.seed,before)
  xy <- as.matrix(points[c("x","y")])
  expect_equal(seg_default_bands(xy,population=points$residents,n=25,probs=probs),actual)
  spatial <- sf::st_as_sf(points,coords=c("x","y"),sf_column_name="position",crs=3857)
  expect_equal(seg_default_bands(spatial,population="residents",n=25,probs=probs),actual)
  set.seed(1)
  chosen <- sample.int(nrow(points),25,prob=points$residents)
  expected <- as.numeric(stats::quantile(as.numeric(stats::dist(xy[chosen,])),
                                         probs=probs,type=8,names=FALSE))
  expect_equal(actual,expected)
  assign(".Random.seed",before,envir=.GlobalEnv)
})

test_that("band sampling does not create RNG state when none exists", {
  (function() {
    before <- if(exists(".Random.seed",.GlobalEnv,inherits=FALSE)) .Random.seed else NULL
    on.exit(if(is.null(before)) {
      if(exists(".Random.seed",.GlobalEnv,inherits=FALSE)) rm(".Random.seed",envir=.GlobalEnv)
    } else assign(".Random.seed",before,envir=.GlobalEnv))
    if(!is.null(before)) rm(".Random.seed",envir=.GlobalEnv)
    bands <- seg_default_bands(cbind(x=1:600,y=0),n=20)
    expect_false(exists(".Random.seed",.GlobalEnv,inherits=FALSE))
    expect_length(bands,9)
  })()
})

test_that("default bands validate geometry, sampling, and population inputs", {
  xy <- cbind(x=1:4,y=0)
  expect_error(seg_default_bands(xy,n=1),"at least 2")
  expect_error(seg_default_bands(xy,probs=c(-.1,.5)),"probabilities")
  expect_error(seg_default_bands(xy,population=rep(0,4)),"positive")
  expect_error(seg_default_bands(xy,population=c(1,NA,1,1)),"finite")
  expect_error(seg_default_bands(cbind(x=rep(0,4),y=0)),"failed to calculate")
  expect_error(seg_default_bands(cbind(x=c(0, 1e200), y=0)),
               "finite pairwise distances")
})

test_that("automatic bands match the public calculator without advancing RNG", {
  points <- data.frame(x = seq_len(600), y = 0, residents = rep(1:3, 200),
                       group = rep(c("A", "B"), 300))
  expected <- seg_default_bands(points, coords = c("x", "y"),
                                  population = "residents")
  withr::with_seed(42, {
    before <- .Random.seed
    result <- seg_measure(points, list(group = seg_categorical("group")),
                    coords = c("x", "y"), population = "residents",
                    measures = "information", scope = "overall")
    expect_identical(result$bands, expected)
    expect_identical(.Random.seed, before)
  })
})
