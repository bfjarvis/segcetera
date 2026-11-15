# Resolve optional values without evaluating the fallback unless needed.
.seg_or_default <- function(x, default) if (is.null(x)) default else x

# Shared by seg_measure() and seg_default_bands() to align attributes with one
# coordinate pair per observation. Also handles optional network coordinates.
.seg_geometry <- function(x, data, coords, polygon_points = "point_on_surface") {
  if (inherits(x, "sf")) {
    if (!is.null(coords)) stop("For sf input, use the active geometry instead of 'coords'",
                              call. = FALSE)
    geometry <- sf::st_geometry(x)
    valid_geometry <- sf::st_make_valid(geometry)
    points <- suppressWarnings(sf::st_centroid(valid_geometry))
    polygonal <- sf::st_geometry_type(valid_geometry) %in% c("POLYGON", "MULTIPOLYGON")
    if (polygon_points == "point_on_surface" && any(polygonal)) {
      points[polygonal] <- suppressWarnings(
        sf::st_point_on_surface(valid_geometry[polygonal])
      )
    }
    xy <- sf::st_coordinates(points)
    source <- if (is.null(data)) sf::st_drop_geometry(x) else as.data.frame(data)
    crs <- sf::st_crs(x)
  } else {
    geometry <- NULL
    crs <- sf::st_crs(NA_character_)
    source <- if (is.null(data)) as.data.frame(x) else as.data.frame(data)
    if (is.null(coords)) {
      xy <- as.matrix(x)
      if (ncol(xy) != 2L)
        stop("Supply two coordinate columns, or select them using 'coords'",
             call. = FALSE)
    } else xy <- as.matrix(.seg_select(as.data.frame(x), coords))
  }
  if (!is.numeric(xy) || ncol(xy) != 2L || nrow(xy) < 1L ||
      nrow(xy) != nrow(source) || any(!is.finite(xy)))
    stop("Supply finite two-dimensional coordinates with one location per data row",
         call. = FALSE)
  storage.mode(xy) <- "double"
  list(coords = xy, data = source, geometry = geometry, crs = crs)
}

# Resolves the common population for seg_measure() radius analyses. When no weights
# are supplied, infer them from variable counts or weights that agree by row.
.seg_population <- function(population, variables, data, count_data) {
  n <- nrow(data)
  if (!is.null(population)) return(.seg_mass(population, data, n))
  candidates <- Map(function(spec, counts) {
    if (!is.null(counts)) rowSums(counts) else
      if (!is.null(spec$weights)) .seg_mass(spec$weights, data, n) else NULL
  }, variables, count_data)
  candidates <- Filter(Negate(is.null), candidates)
  if (!length(candidates)) return(rep(1, n))
  first <- candidates[[1L]]
  if (!all(vapply(candidates, function(z) isTRUE(all.equal(z, first,
                    check.attributes = FALSE)), logical(1))))
    stop("Variable populations differ; supply the common neighborhood 'population' explicitly",
         call. = FALSE)
  first
}

# Shared by seg_measure(), seg_default_bands(), and seg_network() to resolve
# column names or positions consistently while preserving table dimensions.
.seg_select <- function(data, columns) {
  if (is.character(columns)) {
    if (anyNA(columns) || !all(columns %in% names(data)))
      stop("Unknown column: ", paste(setdiff(columns, names(data)), collapse = ", "),
           call. = FALSE)
    index <- match(columns, names(data))
  } else if (is.numeric(columns) && !anyNA(columns) &&
             all(is.finite(columns)) && all(columns >= 1 & columns <= ncol(data)) &&
             all(columns == as.integer(columns))) {
    index <- as.integer(columns)
  } else stop("Column selectors must be names or valid integer positions", call. = FALSE)
  if (anyDuplicated(index)) stop("Columns must not be repeated", call. = FALSE)
  data[, index, drop = FALSE]
}

# Resolves row-aligned weights from a vector, a column, or a trusted default.
# Used by seg_measure() for population mass and by seg_default_bands() for sampling.
.seg_mass <- function(value, data, n, default = rep(1, n)) {
  if (is.null(value)) return(default)
  if (is.character(value) && length(value) == 1L)
    value <- .seg_select(data, value)[[1L]]
  if (!is.numeric(value) || length(value) != n ||
      any(!is.finite(value) | value < 0)) {
    stop("Population weights must be finite, nonnegative, and have one value per row",
         call. = FALSE)
  }
  as.numeric(value)
}

# Keeps retention checks consistent between seg_measure() and the three variable
# constructors. Only a per-variable setting can use NULL to inherit.
.seg_validate_retention <- function(keep_scores, keep_contributions,
                                      allow_null = FALSE) {
  options <- list(keep_scores = keep_scores, keep_contributions = keep_contributions)
  for (name in names(options)) {
    value <- options[[name]]
    if (allow_null && is.null(value)) next
    if (!is.logical(value) || length(value) != 1L || is.na(value)) {
      expected <- if (allow_null) "NULL, TRUE, or FALSE" else "TRUE or FALSE"
      stop("'", name, "' must be ", expected, call. = FALSE)
    }
  }
}
