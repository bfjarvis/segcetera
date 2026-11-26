#' Spatial and Network Segregation Measures
#'
#' Calculates spatial segregation measures from point coordinates or [sf::sf]
#' geometries, or a supplied network, using local population environments.
#'
#' `seg_measure()` is the primary interface for the package. Local environments
#' and segregation indices are produced by a unified C++ engine. The default
#' calculation uses normalized biweight local environments with `power = 3`
#' and radius neighborhoods.
#'
#' By default, `seg_measure()` returns indices over default distance bands. These
#' default bands are distance quantiles from a reproducible population-weighted
#' sample of support points; see [seg_default_bands()] for the quantile defaults.
#' If `output = "indices"`, coordinates, input data, and local environments
#' are not stored in the result.
#'
#' Set `bands = 0` to compute indices without spatial weighting. In that
#' case, local environments contain the observed composition at coincident
#' support points, or just the focal unit when coordinates are distinct.
#'
#' When `neighbors = "knn"`, `bands` is interpreted as population-count
#' thresholds before distance weighting. Each variable
#' selects the requested eligible population, truncating the focal unit when
#' needed and including equidistant boundary units proportionately. Neighboring
#' units are added by representative-point distance. An unattainable threshold raises an error.
#'
#' Each measurement stores named variable results in `analyses`. With `case`,
#' measurements are stored in `result$cases`; an ordinary call returns the
#' measurement directly, with variables in `result$analyses`. Its overall
#' measures are stored in `indices$overall` and requested pairwise
#' measures in `indices$pairwise`.
#' Pairwise scalar measures are returned as one symmetric matrix per
#' bandwidth with one row and column for each population group. Units with
#' no population in either group for a given pair do not contribute to that
#' pair's accumulators.
#'
#' The `surface` argument controls population-surface construction before
#' local environments are calculated. `surface = "grid"` redistributes
#' areal counts over a square grid. `surface = "pycno"` applies internal
#' groupwise pycnophylactic interpolation: each group is smoothed separately
#' by averaging each grid cell with neighboring cells, then rescaled within
#' each source polygon after every iteration so known group totals are
#' preserved.
#'
#' When computing the spatial information theory index, zero proportions are
#' evaluated using the limiting value $x \\log(x) = 0$ as $x$ approaches
#' zero. This avoids artificial perturbation of the data and prevents
#' `0 * log(0)` from producing `NaN`.
#'
#' @param x A numeric matrix or data frame with two coordinate columns, or
#'   an object of class [sf::sf]. Polygon inputs use representative points
#'   selected by `polygon_points` as focal coordinates. With `variables`, a plain
#'   data frame may contain both attributes and the columns selected by `coords`.
#'   With `network`, supply an observation table; coordinates are optional.
#' @param data Optional attribute data with one row per coordinate or feature.
#'   This may contain category memberships, observed
#'   numerical values, and counts selected by the specifications. If omitted,
#'   attributes come from `x` (excluding sf geometry).
#' @param case Optional single column name or position identifying repeated
#'   measurements on common geography. Each case has its own population,
#'   distributions, and normalizers; searches are shared across cases. Support
#'   rows are matched by exact geometry, coordinates, or network node.
#'   Network input uses the node together with any supplied geometry/coordinates.
#'   Repeated rows at a location are matched by occurrence within each case;
#'   absent rows have zero population. Omitted radius bands use the pooled
#'   case populations to select one shared set of distance quantiles. Results
#'   retain individual `seg_result` objects in `cases`; extracted tables identify
#'   their case in a `case` column. With `case = NULL`, the usual result structure
#'   and tables are unchanged. Grid/pycnophylactic surfaces are prepared per case
#'   on the common source geometry, then aligned for shared searches.
#' @param measures A character vector naming measures to compute: "all",
#'   "exposure", "information", "diversity", "dissimilarity", or "square_root".
#'   `"all"` includes all five families. Aliases `"relative_diversity"` and
#'   `"variance_ratio"` select diversity.
#' @param bands Neighborhood thresholds. With `neighbors = "radius"`, these
#'   are distance bandwidths in the units of `x`. With `neighbors = "knn"`,
#'   these are population-count thresholds. If omitted for radius
#'   neighborhoods, [seg_default_bands()] is used.
#' @param weighting Spatial weighting scheme used within each neighborhood.
#'   One of "biweight", "unweighted", "inverse", or "exponential".
#' @param power Power parameter for the biweight and inverse-distance
#'   weighting schemes.
#' @param normalize Logical. If `TRUE`, distances are divided by the
#'   bandwidth before weights are calculated. The bandwidth still defines
#'   the search radius when `neighbors = "radius"`.
#' @param neighbors Neighborhood strategy. "radius" uses distance
#'   bandwidths. "knn" uses count-based adaptive neighborhoods.
#' @param search Neighbor-search backend. The default, "auto", selects "grid"
#'   for regular axis-aligned lattices, "brute" for other Euclidean radius
#'   searches, and "boost" for other Euclidean kNN searches and all network
#'   searches. Automatic lattice detection requires at least one observation per
#'   16 bounding-grid cells on average; explicit `grid_size` also identifies
#'   known sparse grids. Grid/pycnophylactic surfaces use "grid" automatically.
#'   "boost" can be faster for narrow radius neighborhoods; choose an explicit
#'   backend to compare timings for your application. "boost" uses the exact
#'   Boost.Geometry R-tree or Boost.Graph traversal. "brute" computes Euclidean
#'   distances directly. "grid" reuses distance-ordered offsets on an
#'   axis-aligned square or rectangular lattice. Missing cells contribute zero;
#'   repeated observations at a cell remain separate. The backend precomputes
#'   radius weights across locations; kNN retains population-based selection,
#'   proportional distance ties, and locally adaptive weighting. Boundaries are
#'   finite, with offsets outside the input grid skipped. Results are returned
#'   only for supplied focal rows. The offset and weight cache is limited to
#'   512 MiB; use smaller/fewer bands or another backend if
#'   the cache would exceed this budget. Sparse grids can favor "boost".
#'   The reusable-offset approach follows EquiPop (Osth, 2014).
#' @param grid_size Optional positive cell dimensions for `search = "auto"` or
#'   `search = "grid"`:
#'   one value for square spacing, or `c(x_spacing, y_spacing)` for rectangular
#'   spacing. Raw inputs default to the smallest positive coordinate gap on each
#'   axis; specify the known spacing for grids with gaps or large coordinate
#'   offsets. Single-row/column grids use the available axis spacing, and one
#'   unique location defaults to unit spacing. Constructed surfaces supply their
#'   own cell dimensions. Alignment permits floating-point roundoff; arbitrary
#'   points are never silently gridded. Distances use the validated cell spacing.
#' @param surface Population surface used before local environments are
#'   calculated. "raw" uses the input support directly. "grid"
#'   redistributes areal counts over a square grid. "pycno" applies
#'   groupwise pycnophylactic smoothing on that grid while preserving source
#'   polygon group totals.
#' @param surface_audit Logical. If `TRUE`, retain auditable surface
#'   intermediates in `result$surface_info$audit`, including source counts and
#'   geometry, source-zone assignments for retained grid cells, fallback
#'   contributions, the initial gridded counts, and pycnophylactically smoothed
#'   counts when applicable. The default retains only lightweight grid and
#'   fallback metadata.
#' @param scope Whether to calculate "overall" measures, "pairwise"
#'   two-group measures, or "both" for categorical/ordinal categories and resolved
#'   continuous global groups. Ordinal pairwise indices compare individual
#'   categories using the same pair-conditioned formulas as categorical input.
#'   Overall indices use the full categorical distribution, ordinal cumulative
#'   thresholds, or continuous rank distribution. Continuous pairwise indices
#'   use binary group shares within the original full-variable neighborhoods,
#'   independently of summary mode. Exposure is selected through `measures`
#'   independently of scope. This global default can be overridden per variable.
#' @param output Return mode. "indices" stores only segregation indices.
#'   "full" stores both indices and local environments. "localenv"
#'   stores local environments without calculating indices.
#'   Categorical and ordinal environments default to category shares; continuous
#'   variables default to shares of resolved global quartile groups. Use
#'   `summary = list(type = "quantiles")` for synchronized local/grouped rank
#'   summaries, with `mean_rank = TRUE` in that list for a separate mean.
#'   Customize continuous summaries with `summary`, or ordinal summaries with
#'   `localenv`, independently of full-resolution scalar indices.
#'   This global default can be overridden per variable. Shared input geometry,
#'   coordinates, data, and population are retained once if any variable stores
#'   environments.
#' @param keep_scores Logical. Retain per-location H/R/S local variation and
#'   D discrepancies, for the requested measures and scopes. Extract with
#'   `as.data.frame(result, what = "scores")`.
#' @param keep_contributions Logical. Retain additive per-location index
#'   contributions, extracted with `what = "contributions"`. They sum to the
#'   corresponding regional index, using original focal-population weights and
#'   pair-specific weights for pairwise results. Both retention options default
#'   to `FALSE` and work independently of `output`; only requested arrays are
#'   allocated. Exposure and rank profiles have separate output.
#'   Retention memory grows with locations, bands, requested scalar measures,
#'   and group pairs; each enabled option stores its own arrays.
#'   Both retention options can be overridden per variable, including by `FALSE`.
#' @param verbose Logical. If `TRUE`, report progress messages.
#' @param variables Column names for a single categorical variable, or a uniquely
#'   named list of [seg_variables] specifications. A single column name is
#'   interpreted as observed category memberships; two or more column
#'   names are interpreted as category counts (aggregated data). A message reports
#'   this assumption. The result is named after the single column, or `categorical`
#'   for multiple count columns. Use an explicit specification for other variable
#'   types, a custom result name, or population weights for long-form counts.
#'   Single-column count input is not supported. All variables share spatial
#'   discovery; kNN uses each variable's eligible population and adaptive
#'   distance weights. Each specification can override `measures`.
#'   Grid and pycnophylactic surfaces require count specifications for all
#'   variables and polygonal sf input; raw support also accepts observed values.
#'   See `vignette("segcetera-measure-definitions")` for formulas and sources,
#'   and `vignette("segcetera-multiple-variables")` for workflows.
#' @param coords With plain data frames,
#'   two coordinate column names or positions. For sf input, omit this argument
#'   to use the active geometry. Distances are Euclidean in coordinate units;
#'   use an appropriate projected CRS for metric bandwidths.
#' @param polygon_points Representative locations for polygonal `sf` inputs.
#'   Defaults to `"point_on_surface"`, using [sf::st_point_on_surface()] to
#'   select interior points. Use `"centroid"` for [sf::st_centroid()]. This
#'   controls raw Euclidean polygon neighborhoods and their default bands.
#'   Grid surfaces use grid-cell centers; their polygon-allocation and interior
#'   fallback rules are independent of this option. Point inputs and network
#'   shortest-path distances are unaffected. Stored polygon geometry is preserved.
#' @param network Optional [seg_network()] specification. If omitted,
#'   neighborhoods use Euclidean distances. Supplying a network uses shortest
#'   paths exclusively, regardless of any geometry. Network calculations require
#'   explicit `bands` and `surface = "raw"`; `search = "auto"` selects boost,
#'   which can also be chosen explicitly. Radius bands
#'   are in hops or cost units; kNN bands remain variable-specific populations.
#'   Unreachable population does not count toward kNN; an unattainable target
#'   raises an error identifying the variable and focal observation. The global
#'   reference distribution still includes all supplied observations.
#' @param node_id Single attribute column name or position mapping each
#'   observation to a network node. Required with `network`; otherwise omit.
#' @param population Default population-weight column or numeric vector for
#'   observed-value specifications.
#'   Descriptor `weights` override it; count specifications use their row totals.
#'   kNN defaults to one per observed-value row and selects each variable's
#'   own eligible population after omissions. For radius neighborhoods,
#'   consistent count/weight totals can supply the default; specify this
#'   reference population explicitly when those totals differ. It also controls
#'   population-based default distance bands for radius neighborhoods.
#'   Network observed-value specifications default to one per row unless weights
#'   are supplied; count specifications always use their own row totals.
#' @param ... Additional population-surface options, such as `cellsize`,
#'   `celldim`, `nrow`, `ncol`, `max_iter`, `converge`,
#'   `surface_geometry`. Grid and pycnophylactic surfaces use
#'   [terra::rasterize()] to assign source polygons to generated grid-cell
#'   centroids. Non-empty polygons that receive no centroid-assigned cell
#'   contribute their counts to the cell containing an interior point of the
#'   polygon. Set `surface_geometry = "points"` to store centroid point
#'   geometries, or `surface_geometry = "polygons"`/`TRUE` to construct and
#'   store square grid-cell geometries for mapping. By default, generated
#'   surfaces store only centroid coordinates and population counts.
#'
#' @return A [seg_result] with named variable results in `analyses`, plus
#'   shared `bands`, `neighbors`, `diagnostics`, and `surface_info`. Input data,
#'   coordinates, geometry, and environments are retained for `output = "full"`
#'   or `"localenv"`. Surface audits are retained when requested.
#'
#' @note
#' In shares mode the exposure/isolation measure, P, is presented in matrix form. The spatial
#' exposure of group 'm' to group 'n' is located in row 'm' and column 'n'.
#' The matrix is rarely symmetric in practice.
#' In continuous quantiles mode, `measures = "exposure"` requests focal-group
#' averages of local rank statistics, stored as `rank_summary` rather than P.
#' These descriptive profiles are neither P matrices nor segregation indices.
#'
#' @references
#' Osth, J. (2014) Introducing the EquiPop software: an application for the
#'   calculation of k-nearest neighbour contexts/neighbourhoods. Uppsala University.
#'   \url{https://www.uu.se/en/department/human-geography/research/equipop}
#'
#' Osth, J., Clark, W. A. V., and Malmberg, B. (2015) Measuring the scale of
#'   segregation using k-nearest neighbor aggregates. *Geographical Analysis*,
#'   **47**, 34-49. \doi{10.1111/gean.12053}
#'
#' Reardon, S. F. and O'Sullivan, D. (2004) Measures of spatial segregation.
#'   *Sociological Methodology*, **34**, 121-162.
#'
#' Reardon, S. F., Farrell, C. R., Matthews, S. A., O'Sullivan, D., Bischoff,
#'   K., and Firebaugh, G. (2009) Race and space in the 1990s: Changes in the
#'   geographic scale of racial residential segregation, 1990-2000.
#'   *Social Science Research*, **38**, 55-70.
#'
#' @seealso [seg_default_bands()], [seg_surface()], [as.data.frame.seg_result()]
#'
#' @examples
#' data(seg_data)
#' xy <- as.matrix(expand.grid(x = seq(0.5, 9.5), y = seq(0.5, 9.5)))
#' variables <- list(ethnicity = seg_categorical(vars = c("A1", "A2")))
#' result <- seg_measure(xy, variables, data = seg_data, bands = c(0, 2), output = "full")
#' as.data.frame(result)
#' head(as.data.frame(result, what = "env"))
#' result$analyses$ethnicity$env[[1]]
#'
#' # Shorthand for one categorical variable stored as category counts:
#' seg_measure(xy, c("A1", "A2"), data = seg_data, bands = c(0, 2))
#' # For microdata, use variables = "race" to select observed memberships.
#'
#' geom <- sf::st_polygon(list(rbind(c(0,0), c(10,0), c(10,10),
#'   c(0,10), c(0,0)))) |> sf::st_make_grid(cellsize = 1)
#' blocks <- sf::st_sf(seg_data, geometry = geom)
#' gridded <- seg_measure(blocks, variables, surface = "grid", nrow = 10, ncol = 10,
#'                  bands = 2, output = "localenv")
#' sf::st_as_sf(gridded)
#'
#' # Reuse offsets and radius weights on a regular coordinate grid:
#' seg_measure(xy, variables, data = seg_data, bands = c(0, 2),
#'             search = "grid", grid_size = 1)
#'
#' @export
seg_measure <- function(
  x,
  variables,
  coords = NULL,
  data = NULL,
  case = NULL,
  population = NULL,
  measures = "all",
  bands = NULL,
  weighting = c("biweight", "unweighted", "inverse", "exponential"),
  power = 3,
  normalize = TRUE,
  neighbors = c("radius", "knn"),
  search = c("auto", "brute", "boost", "grid"),
  grid_size = NULL,
  surface = c("raw", "grid", "pycno"),
  surface_audit = FALSE,
  scope = c("both", "overall", "pairwise"),
  output = c("indices", "full", "localenv"),
  keep_scores = FALSE,
  keep_contributions = FALSE,
  verbose = FALSE,
  network = NULL,
  node_id = NULL,
  polygon_points = c("point_on_surface", "centroid"),
  ...
) {
  if (missing(variables)) {
    stop("Supply column names or a named 'variables' list of seg_variable specifications",
         call. = FALSE)
  }
  if (is.character(variables) && is.null(dim(variables))) {
    if (!length(variables) || anyNA(variables) || any(!nzchar(variables)) ||
        anyDuplicated(variables)) {
      stop("'variables' column names must be non-empty, non-missing, and unique",
           call. = FALSE)
    }
    columns <- unname(variables)
    if (length(columns) == 1L) {
      variables <- setNames(list(seg_categorical(vars = columns)), columns)
      message("Assuming categorical observed memberships: '", columns,
              "' contains observed category memberships; result named '",
              columns, "'.")
    } else {
      variables <- list(categorical = seg_categorical(vars = columns))
      message("Assuming aggregated categorical data: selected columns contain ",
              "category counts for one variable; result named 'categorical'.")
    }
  }
  call <- match.call()
  dots <- list(...)
  if (!is.list(variables) || !length(variables) || is.null(names(variables)) ||
      anyNA(names(variables)) || any(!nzchar(names(variables))) ||
      anyDuplicated(names(variables)) ||
      !all(vapply(variables, inherits, logical(1), "seg_variable")))
    stop("'variables' must be a uniquely named list of seg_variable specifications",
         call. = FALSE)
  neighbors <- match.arg(neighbors, c("radius", "knn"))
  search <- match.arg(search, c("auto", "brute", "boost", "grid"))
  if (!is.null(grid_size) && !search %in% c("auto", "grid"))
    stop("'grid_size' requires search = 'auto' or 'grid'", call. = FALSE)
  weighting <- match.arg(weighting, c("biweight", "unweighted", "inverse", "exponential"))
  surface <- match.arg(surface, c("raw", "grid", "pycno"))
  polygon_points <- match.arg(polygon_points)
  relational <- !is.null(network)
  if (relational && surface != "raw")
    stop("Network input requires surface = 'raw'", call. = FALSE)
  scope <- match.arg(scope, c("both", "overall", "pairwise"))
  output <- match.arg(output, c("indices", "full", "localenv"))
  .seg_validate_retention(keep_scores, keep_contributions)
  if (!is.logical(surface_audit) || length(surface_audit) != 1L || is.na(surface_audit))
    stop("'surface_audit' must be TRUE or FALSE", call. = FALSE)
  if (surface == "raw" && length(dots))
    stop("Surface options require surface = 'grid' or 'pycno'", call. = FALSE)
  if (surface != "raw" && any(vapply(variables, function(z) is.null(z$counts), logical(1))))
    stop("Grid and pycnophylactic surfaces require category counts in multiple 'vars' columns for every variable",
         call. = FALSE)
  if (length(power) != 1L || !is.numeric(power) || !is.finite(power) || power <= 0)
    stop("'power' must be a positive finite number", call. = FALSE)
  if (!is.logical(normalize) || length(normalize) != 1L || is.na(normalize))
    stop("'normalize' must be TRUE or FALSE", call. = FALSE)
  prepared_network <- NULL
  if (relational) {
    if (search == "auto" && is.null(grid_size)) search <- "boost"
    if (search != "boost") stop("Network input requires search = 'boost'", call. = FALSE)
    if (is.null(bands)) stop("Supply explicit 'bands' for network neighborhoods", call. = FALSE)
    prepared_network <- .seg_network_input(x, data, coords, network, node_id)
    checked <- prepared_network$checked
  } else {
    if (!is.null(node_id)) stop("'node_id' requires a network", call. = FALSE)
    checked <- .seg_geometry(x, data, coords, polygon_points)
  }
  supplied_population <- !is.null(population)
  # Normalize row-aligned inputs and data-independent options once. Case slices
  # and surface recompilation reuse these validated values.
  n <- nrow(checked$data)
  if (supplied_population) population <- .seg_mass(population, checked$data, n)
  valid_measures <- c("exposure", "information", "diversity", "dissimilarity", "square_root")
  variables <- lapply(variables, function(spec) {
    if (!is.null(spec$weights)) spec$weights <- .seg_mass(spec$weights, checked$data, n)
    selected <- .seg_or_default(spec$measures, measures)
    selected[selected %in% c("relative_diversity", "variance_ratio")] <- "diversity"
    if (!is.character(selected) || !length(selected) || anyNA(selected) ||
        any(!selected %in% c(valid_measures, "all")))
      stop("Unknown or empty measure selection", call. = FALSE)
    spec$measures <- if ("all" %in% selected) valid_measures else unique(selected)
    spec
  })
  case_info <- if (!is.null(case)) .seg_case_inputs(checked, case, variables,
    population, if (relational) prepared_network$graph)
  inputs <- if (is.null(case_info)) {
    list(list(checked = checked, variables = variables, population = population))
  } else case_info$inputs
  layout <- if (!is.null(case_info)) .seg_case_layout(inputs)
  prepared_cases <- lapply(seq_along(inputs), function(k) {
    input <- inputs[[k]]
    prepared <- tryCatch(.seg_prepare_population(input$checked, input$variables,
      input$population, neighbors, relational), error = function(e) {
        if (is.null(case_info)) stop(e)
        stop("Case '", as.character(case_info$ids[k]), "': ", conditionMessage(e),
             call. = FALSE)
      })
    if (!is.null(layout))
      prepared <- .seg_align_case(prepared, layout$positions[[k]], layout)
    if (surface != "raw") {
      if (neighbors == "knn" && !supplied_population)
        prepared$population <- prepared$compiled[[1L]]$mass
      source <- if (is.null(layout)) x else
        sf::st_sf(geometry = layout$geometry)
      smoothed <- .seg_surface_variables(source, prepared$checked,
        prepared$variables, prepared$compiled, prepared$population,
        surface, dots, surface_audit, verbose)
      prepared[names(smoothed)] <- smoothed
      if (!is.null(layout) && surface_audit)
        prepared$surface_info$source_rows <- input$rows[
          match(seq_len(layout$n), layout$positions[[k]])]
    }
    prepared
  })
  if (!is.null(case_info) && surface != "raw") {
    # Fallback cells can differ by case; their union supplies shared grid support.
    layout <- .seg_case_layout(prepared_cases)
    prepared_cases <- Map(.seg_align_case, prepared_cases,
                          layout$positions, MoreArgs = list(layout = layout))
  }
  if (!is.null(case_info)) {
    checked <- prepared_cases[[1L]]$checked
    if (relational) {
      prepared_network$graph$node <- layout$nodes
      prepared_network$ids <- prepared_network$graph$labels[layout$nodes + 1L]
    }
  } else checked <- prepared_cases[[1L]]$checked
  population <- Reduce(`+`, lapply(prepared_cases, `[[`, "population"))
  compiled <- do.call(c, lapply(prepared_cases, `[[`, "compiled"))
  flat_variables <- do.call(c, lapply(prepared_cases, `[[`, "variables"))
  if (!is.null(case_info)) {
    names(compiled) <- make.unique(paste(
      rep(as.character(case_info$ids), each = length(variables)),
      names(compiled), sep = ":"))
  }
  surface_info <- prepared_cases[[1L]]$surface_info
  if (is.null(bands)) {
    if (neighbors == "knn") stop("Supply population-count 'bands' for kNN", call. = FALSE)
    bands <- .seg_distance_bands(checked$coords, population)
  } else if (!is.numeric(bands) || !length(bands) || any(!is.finite(bands)) ||
             any(bands < 0) || (neighbors == "knn" && any(bands <= 0))) {
    stop("Supply finite nonnegative radius bands or positive population-count kNN bands",
         call. = FALSE)
  }
  if (neighbors == "knn") {
    for (name in names(compiled)) {
      available <- sum(compiled[[name]]$mass)
      # Allow accumulated rounding of fractional weights in R's sum; the
      # native selector separately verifies the population actually selected.
      tolerance <- (length(compiled[[name]]$mass) + 16) *
        .Machine$double.eps * max(bands)
      if (max(bands) > available + tolerance)
        stop("Variable '", name, "' has eligible population ", format(available),
             ", below the requested kNN threshold ", format(max(bands)), call. = FALSE)
    }
  }
  if (verbose) message("Calculating ", length(compiled), " variables with shared neighborhood searches")
  grid <- NULL
  if (search %in% c("auto", "grid")) {
    if (surface != "raw") {
      surface_size <- c(surface_info$grid$cell_width, surface_info$grid$cell_height)
      if (!is.null(grid_size) && !isTRUE(all.equal(rep_len(grid_size, 2L), surface_size)))
        stop("'grid_size' must match the constructed surface's cell dimensions", call. = FALSE)
      grid_size <- surface_size
    }
    # Unique-axis counts provide a cheap lower bound on bounding-grid cells.
    # Most arbitrary point sets can bypass lattice validation entirely.
    candidate_grid <- search == "grid" || !is.null(grid_size) ||
      as.double(length(unique(checked$coords[, 1]))) *
        length(unique(checked$coords[, 2])) <= 16 * nrow(checked$coords)
    if (candidate_grid) {
      grid <- .seg_grid_layout(checked$coords, grid_size,
                              strict = search == "grid" || !is.null(grid_size))
    }
    if (search == "auto") {
      # Match the native dense-lookup bound; explicit spacing identifies sparse
      # lattices whose large bounding rectangle would otherwise favor a tree.
      if (!is.null(grid) && is.null(grid_size) &&
          prod(as.double(grid$dimensions)) > 16 * nrow(checked$coords)) grid <- NULL
      search <- if (!is.null(grid)) "grid" else
        if (neighbors == "radius") "brute" else "boost"
    }
  }
  if (verbose) message("Using ", search, " neighborhood search")
  compiled <- Map(function(config, spec) {
    config$settings <- list(scope = .seg_or_default(spec$scope, scope),
      output = .seg_or_default(spec$output, output),
      keep_scores = .seg_or_default(spec$keep_scores, keep_scores),
      keep_contributions = .seg_or_default(spec$keep_contributions, keep_contributions))
    config
  }, compiled, flat_variables)
  engine <- .seg_engine_cpp(.seg_or_default(checked$coords, matrix(numeric(), 0L, 2L)),
                         compiled, as.numeric(bands),
                         power,
                         match(weighting, c("unweighted", "biweight",
                                            "inverse", "exponential")) - 1L,
                         normalize,
                         match(neighbors, c("radius", "knn")) - 1L,
                         match(search, c("boost", "brute", "grid")) - 1L,
                         if (relational) prepared_network$graph else list(), grid)
  .seg_prepare_result(engine, prepared_cases, compiled,
    settings = list(bands = as.numeric(bands), weighting = weighting, power = power,
      normalize = normalize, surface = surface, output = output,
      keep_scores = keep_scores, keep_contributions = keep_contributions,
      polygon_points = polygon_points,
      neighbors = list(type = neighbors, values = bands,
        space = if (relational) "network" else "euclidean",
        units = if (neighbors == "knn") "population" else
          if (relational) prepared_network$summary$units else "distance",
        distance_units = if (relational) prepared_network$summary$units else "coordinate",
        population = if (neighbors == "knn") "variable" else NULL,
        engine = search, scope = scope)),
    case_info = case_info, network = prepared_network, call = call)
}

# Prepares one population for seg_measure(). Weights and measure selections
# arrive validated; case-specific totals and inferred populations are checked here.
.seg_prepare_population <- function(checked, variables, population, neighbors, relational) {
  counts <- lapply(variables, function(spec) {
    if (!is.null(spec$counts)) .seg_counts(spec, checked$data)
  })
  if (is.null(population)) {
    candidates <- if (neighbors == "radius" && !relational) {
      Map(function(spec, count) {
        if (!is.null(count)) rowSums(count) else spec$weights
      }, variables, counts)
    } else list()
    candidates <- Filter(Negate(is.null), candidates)
    population <- if (length(candidates)) candidates[[1L]] else rep(1, nrow(checked$data))
    if (!all(vapply(candidates, function(z) isTRUE(all.equal(z, population,
                         check.attributes = FALSE)), logical(1)))) {
      stop("Variable populations differ; supply the common neighborhood 'population' explicitly",
           call. = FALSE)
    }
  }
  total <- if (neighbors == "radius") sum(population)
  if (neighbors == "radius" && (!is.finite(total) || total <= 0))
    stop("The common population total must be positive and finite", call. = FALSE)
  compiled <- Map(function(spec, name, count) {
    .seg_compile(spec, checked$data, population, variable = name, counts = count)
  }, variables, names(variables), counts)
  list(checked = checked, variables = variables, population = population,
       compiled = compiled,
       surface_info = list(type = "raw", fallback = list(count = 0L), audit = NULL))
}

# Packs measurement descriptors into one surface pass for seg_measure(), then
# recompiles their statistics with the original reporting groups. Surface
# construction itself is owned by seg_surface.R.
.seg_surface_variables <- function(x, checked, variables, compiled, population,
                                    surface, dots, audit, verbose) {
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
  prepared <- .seg_surface(x, combined, surface, dots, audit, verbose)
  population <- if (length(reference)) {
    rowSums(prepared$data[, columns[[reference[1L]]], drop = FALSE])
  } else prepared$data[, reference_column]
  recompiled <- Map(function(spec, config, selected) {
    counts <- prepared$data[, selected, drop = FALSE]
    colnames(counts) <- config$metadata$labels
    if (config$type == 3L) spec$scores <- config$scores
    .seg_compile(spec, counts, population,
                   resolution = config$summary_config, counts = counts)
  }, variables, compiled, columns)
  checked$coords <- as.matrix(prepared$coords)
  checked$data <- as.data.frame(prepared$data)
  checked$geometry <- prepared$geometry
  prepared$surface_info$variables <- setNames(columns, names(variables))
  list(checked = checked, population = population, compiled = recompiled,
       surface_info = prepared$surface_info)
}
