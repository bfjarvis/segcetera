# Resolve optional values without evaluating the fallback unless needed.
.seg_or_default <- function(x, default) if (is.null(x)) default else x

# Validates one rectangular lattice for seg_measure() and supplies integer cell
# addresses to the native offset backend. Missing cells retain zero population.
# Non-strict calls return NULL for other layouts during automatic selection.
.seg_grid_layout <- function(coords, size = NULL, strict = TRUE) {
  if (!is.null(size) && (!is.numeric(size) || !length(size) ||
      length(size) > 2L || any(!is.finite(size) | size <= 0))) {
    stop("'grid_size' must contain one or two positive finite cell dimensions",
         call. = FALSE)
  }
  if (is.null(size)) {
    size <- vapply(seq_len(2), function(j) {
      values <- sort(unique(coords[, j]))
      gaps <- diff(values)
      if (!length(gaps)) return(NA_real_)
      # Refine the smallest gap over the full span to avoid cumulative roundoff.
      span <- values[length(values)] - values[1]
      span / round(span / min(gaps))
    }, numeric(1))
    size[is.na(size)] <- if (all(is.na(size))) 1 else size[!is.na(size)][1]
  } else size <- rep_len(as.numeric(size), 2L)
  origin <- apply(coords, 2, min)
  indices <- matrix(0L, nrow(coords), 2L)
  for (j in seq_len(2)) {
    position <- (coords[, j] - origin[j]) / size[j]
    index <- round(position)
    tolerance <- 8 * .Machine$double.eps *
      pmax(abs(coords[, j]), abs(origin[j]), size[j])
    if (any(!is.finite(index)) || any(index > .Machine$integer.max - 1) ||
        any(abs(coords[, j] - (origin[j] + index * size[j])) > tolerance)) {
      if (!strict) return(NULL)
      stop("search = 'grid' requires an axis-aligned rectangular lattice; ",
           "supply its cell dimensions with 'grid_size', or use 'boost'/'brute'",
           call. = FALSE)
    }
    if (length(unique(index)) != length(unique(coords[, j]))) {
      if (!strict) return(NULL)
      stop("Distinct coordinates map to the same grid cell; use 'boost' or 'brute'",
           call. = FALSE)
    }
    indices[, j] <- as.integer(index)
  }
  list(col = indices[, 1], row = indices[, 2], size = size,
       dimensions = as.integer(apply(indices, 2, max) + 1), origin = origin)
}

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

# Prepares one population's statistical inputs for seg_measure(). Case calls
# align spatial support afterward, preserving case-specific regional references.
.seg_compile_input <- function(checked, variables, population, measures,
                               neighbors, relational) {
  counts <- lapply(variables, function(spec) {
    if (!is.null(spec$counts)) .seg_counts(spec, checked$data)
  })
  population <- if (neighbors == "knn" || relational) {
    .seg_mass(population, checked$data, nrow(checked$data))
  } else .seg_population(population, variables, checked$data, counts)
  total <- if (neighbors == "radius") sum(population)
  if (neighbors == "radius" && (!is.finite(total) || total <= 0))
    stop("The common population total must be positive and finite", call. = FALSE)
  compiled <- Map(function(spec, name, count) {
    .seg_compile(spec, checked$data, population, measures, variable = name,
                 counts = count)
  }, variables, names(variables), counts)
  list(checked = checked, variables = variables, population = population,
       compiled = compiled,
       surface_info = list(type = "raw", fallback = list(count = 0L), audit = NULL))
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
