# ------------------------------------------------------------------------------
# Population surface construction for segcetera
# ------------------------------------------------------------------------------

#' Construct a gridded population surface
#'
#' Redistribute polygon counts over a regular grid, with optional
#' pycnophylactic smoothing. This uses the same surface construction as
#' [seg_measure()] and returns the grid without computing neighborhoods or indices.
#'
#' @param x An sf object with polygon or multipolygon geometry.
#' @param vars Names or integer positions of one or more numeric count columns.
#'   Counts must be finite and nonnegative, with one row per polygon.
#' @param surface Use `"grid"` for redistribution or `"pycno"` to smooth
#'   the resulting grid while preserving effective zone totals.
#' @param surface_audit Retain source counts, geometry, cell assignments,
#'   and intermediate grid counts.
#' @param verbose Report construction time.
#' @param ... Grid options: `cellsize` (or `celldim`) in coordinate units,
#'   or both `nrow` and `ncol`. Without either, the target cell size is
#'   one hundredth of the shorter bounding-box dimension.
#'   Set `surface_geometry` to `"points"`, `"polygons"`, or `"none"`
#'   (default). For pycno, `max_iter` and `converge` control iteration and
#'   convergence. See [seg_measure()] for surface and fallback details.
#' @return A list containing grid cell centers (`coords`), a matrix of
#'   redistributed counts (`data`), optional cell `geometry`, source-zone
#'   row numbers (`id`), the source `crs`, and `surface_info` with grid,
#'   fallback, and optional audit metadata. Rows follow the same order.
#' @examples
#' bounds <- sf::st_bbox(c(xmin = 0, ymin = 0, xmax = 2, ymax = 1),
#'                       crs = sf::st_crs(3857))
#' zones <- sf::st_sf(A = c(10, 20), B = c(20, 10),
#'   geometry = sf::st_make_grid(sf::st_as_sfc(bounds), n = c(2, 1)))
#' grid <- seg_surface(zones, c("A", "B"), nrow = 2, ncol = 4,
#'                       surface_geometry = "polygons")
#' sf::st_sf(as.data.frame(grid$data), geometry = grid$geometry)
#' @seealso [seg_measure()]
#' @export
seg_surface <- function(x, vars, surface = c("grid", "pycno"),
                          ..., surface_audit = FALSE, verbose = FALSE) {
  surface <- match.arg(surface)
  if (!inherits(x, "sf") || !nrow(x) ||
      !all(sf::st_geometry_type(x) %in% c("POLYGON", "MULTIPOLYGON")) ||
      any(sf::st_is_empty(x))) {
    stop("'x' must contain non-empty sf polygons or multipolygons", call. = FALSE)
  }
  if (missing(vars) || !length(vars))
    stop("Select at least one count column with 'vars'", call. = FALSE)
  if (!is.logical(surface_audit) || length(surface_audit) != 1L ||
      is.na(surface_audit))
    stop("'surface_audit' must be TRUE or FALSE", call. = FALSE)
  counts <- .seg_counts(list(counts = vars, na.action = "fail"),
                         sf::st_drop_geometry(x))
  result <- .seg_surface(x, coords = NULL, data = counts, surface = surface,
                          args = list(...), audit = surface_audit, verbose = verbose)
  result$crs <- sf::st_crs(x)
  result
}

# Prepares all variables in one gridding pass for seg_measure(), then recompiles
# each variable on the shared grid while retaining its reporting groups.
.seg_surface_variables <- function(x, checked, variables, compiled, population,
                                    measures, surface, dots, audit, verbose) {
  sizes <- vapply(compiled, function(z) ncol(z$counts), integer(1))
  ends <- cumsum(sizes)
  columns <- Map(seq.int, ends - sizes + 1L, ends)
  combined <- do.call(cbind, lapply(compiled, function(z) z$counts))
  colnames(combined) <- make.unique(unlist(lapply(compiled, function(z) z$metadata$labels)))
  reference <- which(vapply(compiled, function(z) identical(z$mass, population), logical(1)))
  if (!length(reference)) {
    reference_column <- ncol(combined) + 1L
    combined <- cbind(combined, population)
    colnames(combined)[reference_column] <- make.unique(c(colnames(combined)[-reference_column],
                                                         ".population"))[reference_column]
  }
  prepared <- .seg_surface(x, checked$coords, combined, surface, dots, audit, verbose)
  population <- if (length(reference)) {
    rowSums(prepared$data[, columns[[reference[1L]]], drop = FALSE])
  } else prepared$data[, reference_column]
  recompiled <- Map(function(spec, config, selected) {
    counts <- prepared$data[, selected, drop = FALSE]
    colnames(counts) <- config$metadata$labels
    if (config$type == 3L) spec$scores <- config$scores
    .seg_compile(spec, counts, population, measures,
                   resolution = config$summary_config, counts = counts)
  }, variables, compiled, columns)
  checked$coords <- as.matrix(prepared$coords)
  checked$data <- as.data.frame(prepared$data)
  checked$geometry <- prepared$geometry
  prepared$surface_info$variables <- setNames(columns, names(variables))
  list(checked = checked, population = population, compiled = recompiled,
       surface_info = prepared$surface_info)
}

#' Construct a population surface for spatial segregation calculations
#'
#' `.seg_surface()` prepares the coordinate and count matrix that will be passed
#' to the segregation engine. The `"raw"` surface returns the input coordinates
#' and data unchanged. The `"grid"` and `"pycno"` surfaces require polygonal
#' `sf` input and redistribute polygon counts to a regular square grid before
#' optionally applying pycnophylactic smoothing.
#' Shared by `seg_measure()` and `seg_surface()` after count preparation.
#'
#' @param x Input object supplied to `seg_measure()`. Grid-based surfaces require an
#'   `sf` object with polygon geometry.
#' @param coords Numeric two-column coordinate matrix extracted from `x`.
#' @param data Numeric matrix of group counts.
#' @param surface Character scalar. One of `"raw"`, `"grid"`, or `"pycno"`.
#' @param args Named list of surface options passed through `...` from `seg_measure()`.
#' @param audit Logical. If `TRUE`, retain surface-construction intermediates.
#' @param verbose Logical. If `TRUE`, report timing messages.
#'
#' @return A list with at least `coords`, `data`, and `geometry`. Grid-based
#'   surfaces may also include source-polygon identifiers and fallback metadata.
#'
#' @noRd
.seg_surface <- function(
  x,
  coords,
  data,
  surface,
  args,
  audit = FALSE,
  verbose = FALSE
) {
  if (surface == "raw") {
    return(list(
      coords = coords,
      data = data,
      geometry = NULL,
      surface_info = list(
        type = "raw",
        fallback = list(count = 0L, ids = integer(), cells = integer()),
        audit = NULL
      )
    ))
  }
  if (!inherits(x, "sf")) {
    stop("surface = '", surface, "' requires 'x' to be an sf object", call. = FALSE)
  }
  allowed <- c("cellsize", "celldim", "nrow", "ncol", "max_iter", "converge",
               "surface_geometry")
  if (length(args) && (is.null(names(args)) || anyNA(names(args)) ||
      any(!nzchar(names(args))) || anyDuplicated(names(args)) ||
      any(!names(args) %in% allowed))) {
    stop("Supply uniquely named surface options: ", paste(allowed, collapse = ", "),
         call. = FALSE)
  }
  geometry_type <- .seg_geometry_type(args)
  if (surface == "pycno") {
    max_iter <- .seg_or_default(args$max_iter, 1000L)
    converge <- .seg_or_default(args$converge, 3)
    if (!is.numeric(max_iter) || length(max_iter) != 1L ||
        !is.finite(max_iter) || max_iter < 1 || max_iter != floor(max_iter)) {
      stop("'max_iter' must be a positive integer", call. = FALSE)
    }
    if (!is.numeric(converge) || length(converge) != 1L || !is.finite(converge)) {
      stop("'converge' must be a finite numeric scalar", call. = FALSE)
    }
  }
  if (verbose) {
    started <- Sys.time()
    message("Constructing ", surface, " population surface ...")
  }
  grid <- .seg_grid(x, data, args, audit)
  .seg_warn_grid_assignment(grid$fallback, surface)
  if (surface == "pycno") {
    group_names <- colnames(grid$values)
    grid$values <- .seg_pycno_cpp(grid$zone_ids, grid$values, grid$ncol,
                                  grid$nrow, max_iter, converge)
    colnames(grid$values) <- group_names
  }
  out <- .seg_compact_grid(grid, st_crs(x), geometry_type, surface)
  if (verbose) {
    elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
    message("Population surface done [", elapsed, " seconds]")
  }
  out
}

#' Build the full raster-backed grid representation
#'
#' Creates a terra raster template, rasterizes source polygon identifiers,
#' distributes polygon counts over assigned grid cells, and adds the counts of
#' non-empty polygons without a centroid-based grid cell to the cell containing
#' a representative interior point.
#'
#' @param x Polygonal `sf` object.
#' @param data Numeric matrix of group counts, one row per source polygon.
#' @param args Named list of surface options.
#' @param audit Logical. If `TRUE`, retain source and initial-grid information.
#'
#' @return A list describing the full grid, including all cell coordinates,
#'   count values, source-zone identifiers, fallback assignments, grid dimensions,
#'   and cell size.
#'
#' @noRd
.seg_grid <- function(x, data, args, audit = FALSE) {
  template <- .seg_grid_template(x, args)
  zone_ids <- .seg_rasterize_zones(x, template)
  fallback <- .seg_fallback_counts(zone_ids, x, data, template)
  values <- .seg_distribute_counts(zone_ids, data)

  if (length(fallback$cells) > 0) {
    fallback_values <- rowsum(fallback$data, fallback$cells, reorder = FALSE)
    fallback_cells <- as.integer(rownames(fallback_values))
    values[fallback_cells, ] <- values[fallback_cells, , drop = FALSE] +
      fallback_values

    unassigned_cells <- fallback_cells[is.na(zone_ids[fallback_cells])]
    for (cell in unassigned_cells) {
      zone_ids[cell] <- fallback$ids[match(cell, fallback$cells)]
    }
  }

  out <- list(
    coords = terra::xyFromCell(template, seq_len(terra::ncell(template))),
    values = values,
    zone_ids = zone_ids,
    fallback = fallback,
    nrow = terra::nrow(template),
    ncol = terra::ncol(template),
    cell_width = terra::xres(template),
    cell_height = terra::yres(template)
  )
  if (audit) {
    out$initial_values <- values
    out$source_data <- data
    out$source_geometry <- st_geometry(x)
  }
  out
}

#' Construct a raster template from grid dimensions or cell size
#'
#' @param x Polygonal sf input.
#' @param args Named surface options: nrow and ncol, or celldim/cellsize.
#' @return A terra raster with the extent and CRS of x.
#' @noRd
.seg_grid_template <- function(x, args) {
  bbox <- st_bbox(x)
  if (!is.null(args$nrow) || !is.null(args$ncol)) {
    dimensions <- list(args$nrow, args$ncol)
    valid <- vapply(dimensions, function(z) {
      is.numeric(z) && length(z) == 1L && is.finite(z) &&
        z >= 1 && z == floor(z)
    }, logical(1))
    if (!all(valid)) {
      stop("Supply both 'nrow' and 'ncol' as positive integers", call. = FALSE)
    }
    rows <- args$nrow
    columns <- args$ncol
  } else {
    celldim <- .seg_or_default(args$celldim, args$cellsize)
    if (is.null(celldim))
      celldim <- min(bbox$xmax - bbox$xmin, bbox$ymax - bbox$ymin) / 100
    if (!is.numeric(celldim) || length(celldim) != 1L ||
        !is.finite(celldim) || celldim <= 0) {
      stop("'celldim' must be a positive finite numeric scalar", call. = FALSE)
    }
    rows <- max(1, ceiling((bbox$ymax - bbox$ymin) / celldim))
    columns <- max(1, ceiling((bbox$xmax - bbox$xmin) / celldim))
  }
  terra::rast(nrows = rows, ncols = columns,
    xmin = bbox$xmin, xmax = bbox$xmax, ymin = bbox$ymin, ymax = bbox$ymax,
    crs = st_crs(x)$wkt)
}

#' Rasterize source polygon identifiers
#'
#' Assigns each grid cell to the polygon containing its centroid. The returned
#' vector uses `NA` for cells outside the source polygons.
#'
#' @param x Polygonal `sf` object.
#' @param template Terra raster template.
#'
#' @return Integer vector of source polygon row numbers, one value per raster
#'   cell.
#'
#' @noRd
.seg_rasterize_zones <- function(x, template) {
  x$.seg_zone_id <- seq_len(nrow(x))
  zone_raster <- terra::rasterize(
    terra::vect(x),
    template,
    field = ".seg_zone_id",
    background = NA,
    touches = FALSE
  )
  as.integer(terra::values(zone_raster, mat = FALSE))
}

#' Distribute polygon counts over assigned grid cells
#'
#' Divides each source polygon's group counts evenly across all grid cells
#' assigned to that polygon.
#'
#' @param zone_ids Integer vector of source polygon identifiers, one per grid
#'   cell.
#' @param data Numeric matrix of group counts, one row per source polygon.
#'
#' @return Numeric matrix of grid-cell group counts.
#'
#' @noRd
.seg_distribute_counts <- function(zone_ids, data) {
  values <- matrix(0, nrow = length(zone_ids), ncol = ncol(data))
  assigned <- !is.na(zone_ids)
  if (!any(assigned)) {
    colnames(values) <- colnames(data)
    return(values)
  }

  ids <- zone_ids[assigned]
  cells_per_polygon <- table(ids)
  for (i in seq_len(ncol(data))) {
    values[assigned, i] <- data[ids, i] / cells_per_polygon[as.character(ids)]
  }
  colnames(values) <- colnames(data)
  values
}

#' Find fallback cells for non-empty polygons without grid-cell centroids
#'
#' Finds source polygons with positive population but no rasterized grid cell,
#' then identifies the grid cell containing a representative interior point.
#' Their counts are subsequently added directly to that grid cell. When several
#' missing polygons select the same cell, their counts are summed.
#' Called by `.seg_grid()` for both `seg_measure()` gridded surface modes.
#'
#' @param zone_ids Integer vector of source polygon identifiers.
#' @param x Polygonal `sf` object.
#' @param data Numeric matrix of group counts.
#' @param template Terra raster template.
#'
#' @return Assigned fallback cells, source polygon ids, and count rows,
#'   plus the ids of polygons for which fallback assignment failed.
#'
#' @noRd
.seg_fallback_counts <- function(zone_ids, x, data, template) {
  positive <- which(rowSums(data) > 0)
  represented <- unique(zone_ids[!is.na(zone_ids)])
  missing <- setdiff(positive, represented)
  if (!length(missing)) {
    return(list(cells = integer(), ids = integer(),
                data = data[integer(), , drop = FALSE], failed_ids = integer()))
  }

  fallback_points <- suppressWarnings(
    st_point_on_surface(st_geometry(x[missing, , drop = FALSE]))
  )
  fallback_xy <- st_coordinates(fallback_points)[, 1:2, drop = FALSE]
  cells <- terra::cellFromXY(template, fallback_xy)
  valid <- !is.na(cells)

  list(cells = as.integer(cells[valid]), ids = missing[valid],
       data = data[missing[valid], , drop = FALSE],
       failed_ids = missing[!valid])
}

#' Drop unused grid cells and attach optional geometry
#'
#' Keeps cells assigned to source polygons and returns the compact
#' representation expected by `seg_measure()`.
#'
#' @param grid Full grid list from `.seg_grid()`.
#' @param crs Coordinate reference system for returned geometry.
#' @param geometry_type Character scalar: `"none"`, `"points"`, or `"polygons"`.
#' @param surface_type Character scalar identifying the constructed surface.
#'
#' @return A compact surface list with coordinates, counts, geometry, source ids,
#'   and fallback metadata.
#'
#' @noRd
.seg_compact_grid <- function(grid, crs, geometry_type, surface_type) {
  keep <- !is.na(grid$zone_ids)
  cell_ids <- which(keep)

  coords <- grid$coords[keep, , drop = FALSE]
  values <- grid$values[keep, , drop = FALSE]
  if (any(!is.finite(values) | values < 0))
    stop("Surface construction must produce finite nonnegative counts", call. = FALSE)
  zone_ids <- grid$zone_ids[keep]
  geometry <- .seg_make_geometry(
    coords,
    crs,
    geometry_type,
    grid$cell_width,
    grid$cell_height
  )

  colnames(coords) <- c("x", "y")
  colnames(values) <- colnames(grid$values)
  surface_info <- list(
    type = surface_type,
    grid = list(
      nrow = grid$nrow,
      ncol = grid$ncol,
      cell_width = grid$cell_width,
      cell_height = grid$cell_height,
      retained_cells = length(cell_ids)
    ),
    fallback = list(
      count = length(grid$fallback$ids),
      ids = grid$fallback$ids,
      cells = grid$fallback$cells,
      failed_ids = grid$fallback$failed_ids
    ),
    audit = NULL
  )

  if (!is.null(grid$initial_values)) {
    fallback_audit <- data.frame(
      source_id = grid$fallback$ids,
      cell_id = grid$fallback$cells,
      zone_id = grid$zone_ids[grid$fallback$cells]
    )
    if (nrow(fallback_audit) > 0) {
      fallback_audit <- cbind(
        fallback_audit,
        as.data.frame(grid$fallback$data)
      )
    }

    surface_info$audit <- list(
      source = list(
        data = grid$source_data,
        geometry = grid$source_geometry
      ),
      cells = data.frame(
        row = seq_along(cell_ids),
        cell_id = cell_ids,
        zone_id = zone_ids,
        x = coords[, 1],
        y = coords[, 2]
      ),
      fallback = fallback_audit,
      initial_counts = grid$initial_values[keep, , drop = FALSE],
      smoothed_counts = if (identical(surface_type, "pycno")) {
        values
      } else {
        NULL
      }
    )
  }

  list(
    coords = coords,
    data = values,
    geometry = geometry,
    id = zone_ids,
    surface_info = surface_info
  )
}

#' Report the outcome of grid fallback assignment
#'
#' Called by `.seg_surface()` so `seg_measure()` can report successes and
#' failures together and explain the smoothing implications for pycno.
#'
#' @param fallback Fallback assignment list from `.seg_fallback_counts()`.
#' @param surface Character scalar identifying the surface type.
#' @return Invisibly returns `NULL`; called for its warning side effect.
#' @noRd
.seg_warn_grid_assignment <- function(fallback, surface) {
  assigned <- length(fallback$ids)
  failed <- length(fallback$failed_ids)
  if (!assigned && !failed) return(invisible(NULL))

  warning(
    if (assigned) paste0(assigned, " populated source polygon(s) received no ",
      "grid-cell centroid and were assigned to fallback grid cells. "),
    if (assigned && surface == "pycno")
      paste0("These polygons are not preserved as separate zonal constraints ",
             "during pycnophylactic smoothing. "),
    if (failed) paste0(failed, " populated source polygon(s) could not be ",
      "assigned to fallback grid cells; their counts are absent from the grid. "),
    "Consider using a finer grid and inspect 'result$surface_info$fallback'.",
    call. = FALSE
  )
  invisible(NULL)
}

#' Normalize the requested surface geometry type
#'
#' Accepts logical and character inputs for `surface_geometry`; logical values
#' select polygon geometry or no stored geometry.
#'
#' @param args Named list of surface options.
#'
#' @return One of `"none"`, `"points"`, or `"polygons"`.
#'
#' @noRd
.seg_geometry_type <- function(args) {
  geometry <- .seg_or_default(args$surface_geometry, FALSE)
  if (is.logical(geometry)) {
    if (length(geometry) != 1 || is.na(geometry)) {
      stop(
        "'surface_geometry' must be TRUE, FALSE, or a geometry type",
        call. = FALSE
      )
    }
    return(if (geometry) "polygons" else "none")
  }

  geometry <- as.character(geometry)
  if (length(geometry) != 1L || is.na(geometry) ||
      !geometry %in% c("none", "points", "polygons")) {
    stop("'surface_geometry' must be 'none', 'points', or 'polygons'",
         call. = FALSE)
  }
  geometry
}

#' Build optional grid-cell geometry
#'
#' @param coords Numeric matrix of grid-cell center coordinates.
#' @param crs Coordinate reference system for the returned geometry.
#' @param type Character scalar: `"none"`, `"points"`, or `"polygons"`.
#' @param width Numeric grid-cell width.
#' @param height Numeric grid-cell height.
#'
#' @return An `sfc` geometry vector, or `NULL` when `type = "none"`.
#'
#' @noRd
.seg_make_geometry <- function(coords, crs, type, width, height = width) {
  if (identical(type, "none")) {
    return(NULL)
  }
  if (nrow(coords) == 0) {
    return(st_sfc(crs = crs))
  }

  if (identical(type, "points")) {
    points <- as.data.frame(coords)
    names(points) <- c("x", "y")
    return(st_geometry(st_as_sf(points, coords = c("x", "y"), crs = crs)))
  }

  cells <- lapply(seq_len(nrow(coords)), function(i) {
    x <- coords[i, 1]
    y <- coords[i, 2]
    st_polygon(list(rbind(
      c(x - width / 2, y - height / 2),
      c(x + width / 2, y - height / 2),
      c(x + width / 2, y + height / 2),
      c(x - width / 2, y + height / 2),
      c(x - width / 2, y - height / 2)
    )))
  })
  st_sfc(cells, crs = crs)
}
