test_that("standalone grid and pycno match the surfaces used by segcetera", {
  bounds <- sf::st_bbox(c(xmin = 0, ymin = 0, xmax = 2, ymax = 1),
                        crs = sf::st_crs(3857))
  zones <- sf::st_sf(A = c(10, 20), B = c(20, 10),
    geometry = sf::st_make_grid(sf::st_as_sfc(bounds), n = c(2, 1)))
  for (surface in c("grid", "pycno")) {
    prepared <- seg_surface(zones, c("B", "A"), surface,
      nrow = 2, ncol = 4, surface_geometry = "polygons", surface_audit = TRUE)
    result <- seg_measure(zones, list(group = seg_categorical(c("B", "A"))),
      surface = surface, bands = 0, output = "full", nrow = 2, ncol = 4,
      surface_geometry = "polygons", surface_audit = TRUE)
    expect_identical(colnames(prepared$data), c("B", "A"))
    expect_equal(colSums(prepared$data), c(B = 30, A = 30))
    expect_equal(prepared$data, as.matrix(result$data))
    expect_equal(prepared$coords, result$coords)
    expect_equal(prepared$geometry, result$geometry)
    expect_equal(prepared$crs, sf::st_crs(zones))
    expect_equal(prepared$surface_info$audit, result$surface_info$audit)
  }

  one_column <- seg_surface(zones, 1L, nrow = 2, ncol = 4)
  expect_equal(colSums(one_column$data), c(A = 30))
  expect_null(one_column$geometry)
  expect_null(one_column$surface_info$audit)
  expect_error(seg_surface(zones), "at least one count column")
  expect_error(seg_surface(zones, "unknown"), "Unknown column")
  expect_error(seg_surface(sf::st_centroid(zones$geometry), "A"),
               "non-empty sf polygons")
  zones$A[1] <- NA_real_
  expect_error(seg_surface(zones, "A"), "Missing counts")
})
