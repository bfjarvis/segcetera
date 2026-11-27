# Builds complete single- or multi-case results for seg_measure(). Native
# analysis labels, retained inputs, and shared metadata are assembled here once.
.seg_prepare_result <- function(engine, prepared_cases, compiled, settings,
                                case_info, network, call) {
  bands <- settings$bands
  analyses <- engine$analyses
  for (i in seq_along(analyses)) {
    analysis <- analyses[[i]]
    config <- compiled[[i]]
    keep_env <- config$settings$output != "indices"
    analysis$settings <- config$settings
    if (!is.null(analysis$scores) || !is.null(analysis$contributions)) {
      labels <- config$pair_labels
      pairs <- if (length(labels) > 1L) utils::combn(labels, 2) else
        matrix(character(), nrow = 2, ncol = 0)
      analysis$score_pairs <- data.frame(group_a = pairs[1, ], group_b = pairs[2, ])
      for (what in c("scores", "contributions")) {
        if (is.null(analysis[[what]])) next
        for (scope in c("overall", "pairwise")) {
          analysis[[what]][[scope]] <- lapply(analysis[[what]][[scope]], function(z) {
            names(z) <- as.character(bands)
            z
          })
        }
      }
    }
    for (group in c("overall", "pairwise")) {
      analysis$indices[[group]] <- lapply(analysis$indices[[group]], function(z) {
        names(z) <- as.character(bands)
        z
      })
    }
    for (key in names(analysis$indices$pairwise)) {
      analysis$indices$pairwise[[key]] <- lapply(analysis$indices$pairwise[[key]], function(z) {
        dimnames(z) <- list(config$pair_labels, config$pair_labels)
        z
      })
    }
    if (length(analysis$rank_summary)) {
      names(analysis$rank_summary) <- as.character(bands)
      analysis$rank_summary <- lapply(analysis$rank_summary, function(z) {
        dimnames(z) <- list(config$exposure_labels, config$metadata$localenv$labels)
        z
      })
    }
    if (keep_env) {
      names(analysis$env) <- as.character(bands)
      analysis$env <- lapply(analysis$env, function(z) {
        colnames(z) <- config$metadata$localenv$labels
        z
      })
      if (length(analysis$interpolation)) {
        names(analysis$interpolation) <- as.character(bands)
        analysis$interpolation <- lapply(analysis$interpolation, function(z) {
          colnames(z) <- c("bottom_category_share", "top_category_share",
                           "interpolated_share")
          z
        })
      }
    }
    metadata <- config$metadata
    # Millions of exact values are calculation inputs, not lightweight metadata.
    if (config$type == 3L && !keep_env) {
      metadata$distinct_values <- sum(config$totals > 0)
      metadata$labels <- metadata$scores <- metadata$shares <- NULL
    }
    analysis$definition <- metadata
    if (!is.null(analysis$knn)) analysis$knn <- data.frame(band = bands, analysis$knn)
    names(analysis$status) <- as.character(bands)
    analyses[[i]] <- analysis
  }
  n_variables <- length(prepared_cases[[1L]]$variables)
  store_input <- any(vapply(compiled, function(z)
    z$settings$output != "indices", logical(1)))
  store_local_keys <- store_input || any(vapply(compiled, function(z)
    z$settings$keep_scores || z$settings$keep_contributions, logical(1)))
  results <- lapply(seq_along(prepared_cases), function(k) {
    prepared <- prepared_cases[[k]]
    indices <- seq.int((k - 1L) * n_variables + 1L, k * n_variables)
    case_analyses <- analyses[indices]
    names(case_analyses) <- names(prepared$variables)
    checked <- prepared$checked
    result <- list(analyses = case_analyses, bands = bands,
      coords = if (store_input) checked$coords else NULL,
      data = if (store_input) .seg_or_default(prepared$input$data, checked$data) else NULL,
      site = if (store_input) .seg_or_default(prepared$input$site,
                                             seq_len(nrow(checked$data))) else NULL,
      geometry = if (store_input) checked$geometry else NULL,
      population = if (store_input) prepared$population else NULL,
      crs = checked$crs, weighting = settings$weighting, power = settings$power,
      normalize = settings$normalize, surface = settings$surface,
      output = settings$output, keep_scores = settings$keep_scores,
      keep_contributions = settings$keep_contributions,
      polygon_points = settings$polygon_points,
      surface_info = prepared$surface_info,
      network = if (!is.null(network)) network$summary else NULL,
      node_id = if (!is.null(network) && store_local_keys) network$ids else NULL,
      neighbors = settings$neighbors,
      diagnostics = engine$diagnostics, call = call)
    class(result) <- "seg_result"
    if (!is.null(case_info)) {
      result$case <- case_info$ids[k]
      if (!is.null(result$data) && case_info$column %in% names(result$data) &&
          anyNA(result$data[[case_info$column]]))
        result$data[[case_info$column]] <- rep(result$case, nrow(result$data))
      result$diagnostics$variables <- n_variables
      result$diagnostics$local_environments <- n_variables * length(bands) *
        nrow(prepared$checked$data)
    }
    result
  })
  if (is.null(case_info)) return(results[[1L]])
  names(results) <- make.unique(as.character(case_info$ids))
  result <- results[[1L]]
  result$case <- result$analyses <- result$data <- result$population <- result$site <- NULL
  result$cases <- results
  result$case_ids <- case_info$ids
  result$surface_info <- list(type = settings$surface,
                             cases = lapply(results, `[[`, "surface_info"))
  result$diagnostics <- engine$diagnostics
  result$diagnostics$analyses <- length(compiled)
  result$diagnostics$variables <- n_variables
  result$diagnostics$cases <- length(results)
  result
}

#' Spatial segregation results
#'
#' [seg_measure()] returns a `seg_result`. Each element of `analyses` contains `indices`, optional
#' `env`, `definition`, `normalizers`, and a status for each band. Stored
#' environments are compact matrices: categorical and ordinal shares, or local
#' shares of global groups for continuous variables by default. Continuous
#' `summary` selects synchronized shares or rank summaries; ordinal `localenv`
#' selects quantiles or category shares. Settings are recorded
#' in `definition$localenv`. The regional reference distribution is retained
#' once when environments are stored. Each analysis records effective `scope`,
#' `output`, `keep_scores`, and `keep_contributions` in `settings`. Global result
#' settings record defaults; variable-specific settings take precedence. Shared
#' inputs are retained once when any analysis stores environments.
#' With a `case` specification, `cases` contains one `seg_result` per case and
#' `case_ids` preserves the original ID values. Each case uses the same support
#' row numbering; absent observations have zero population. Definitions and
#' normalizers remain case-specific. Extracted tables include `case` for these
#' results, including when a single case was supplied. Ordinary results omit it.
#' Local results contain one row per unique spatial site. When input data are
#' retained, `data` preserves the supplied records for repeated locations and
#' `site` maps each stored data row to a row of `coords`, `geometry`, `population`,
#' and the local-result matrices. Scores and additive contributions are site-level;
#' their focal weights are the summed eligible populations at each site.
#' Probability-valued local summaries, exposure/isolation, rank summaries, and
#' interpolation shares are bounded to [0, 1] to correct floating-point endpoint
#' overshoot. Interior values retain full precision; undefined results remain
#' missing. Endpoint corrections use an accumulation-aware roundoff allowance.
#' Larger violations or non-finite computed
#' probabilities raise an error.
#' Each definition records its `eligible_population`. For kNN, each analysis
#' also contains a `knn` table: requested population, minimum/maximum selected
#' population before distance weighting, and minimum/maximum neighborhood
#' radius, by band across focal rows with positive eligible population.
#'
#' @param x A `seg_result`.
#' @param variables Optional names of analyses to extract.
#' @param what Table to extract: `"indices"` for overall scalar indices,
#'   `"pairwise"` for two-group indices, `"exposure"` for focal-group summaries,
#'   `"env"` for stored local summaries, `"diagnostics"` for per-band status and
#'   kNN diagnostics, or `"interpolation"` for stored neighborhood interpolation
#'   diagnostics. `"scores"` extracts retained local variation/deviation scores;
#'   `"contributions"` extracts their additive index contributions. Both identify
#'   `scope`, `location`, `measure`, and pair members where applicable, and require
#'   their corresponding `keep_scores` or `keep_contributions` option in [seg_measure()].
#'   Categorical and ordinal pairwise indices compare original categories;
#'   continuous pairwise indices compare resolved global groups.
#'   Probability-defined continuous groups are labeled by population percentile
#'   ranges, such as `P0-P25`. Explicit value cuts retain value-based labels.
#'   Exposure tables identify the focal population as `reference` and the
#'   neighboring population category as `target`. Matching categories describe
#'   isolation. Rank profiles identify `reference` and a rank statistic instead.
#'   With `variables = NULL`, mixed-output results include analyses with the
#'   requested stored results. Explicitly selecting an analysis without those
#'   results raises an error; so does requesting results unavailable everywhere.
#' @param bands Optional subset of neighborhood bands.
#' @param row.names,optional Standard data-frame generic arguments.
#' @param ... Reserved for generic methods. Extraction rejects unused arguments.
#' @return Printing returns `x` invisibly. Extraction returns a long data frame
#'   with identifiers and one numeric `value` column for results:
#'   \itemize{
#'     \item Indices: `variable`, `band`, `measure`, `value`.
#'     \item Pairwise indices additionally identify `group_a` and `group_b`.
#'     \item Share exposure: `variable`, `band`, `reference`, `target`, `value`.
#'     \item Quantile exposure: `variable`, `band`, `reference`, `statistic`,
#'       `prob`, `value`.
#'     \item Local shares: `variable`, `band`, `location`, `category`, `value`.
#'     \item Local ranks: `variable`, `band`, `location`, `statistic`, `prob`, `value`.
#'     \item Scores/contributions: `variable`, `band`, `scope`, `location`,
#'       `measure`, `value`, plus pair members when pairwise results are requested.
#'   }
#'   Columns are chosen from the selected analyses' summary specifications, not
#'   observed missingness. Mixed share/rank tables contain the union of their
#'   required columns; `statistic = "share"` identifies share rows in those tables.
#'   Quantile rows use `statistic = "quantile_rank"`; optional mean ranks use
#'   `"mean_rank"` and have no quantile probability. Undefined results retain
#'   their rows and missing `value`; empty tables retain their typed schema.
#'
#'   With `what = "exposure"`, shares are exposure/isolation probabilities, but
#'   quantiles are focal-population-weighted averages of local global-rank
#'   statistics, not probabilities of encountering a target group. Mean ranks
#'   summarize the full local rank distribution. These rank profiles are stored
#'   separately in `analyses[[name]]$rank_summary`, not in the P matrix. The
#'   summary specification remains in `analyses[[name]]$definition`.
#'
#'   Type, representation, method, normalizers, and calculation settings remain
#'   in the result object rather than being repeated in result tables.
#'   Retained scores and contributions are stored separately in each analysis,
#'   grouped by scope, measure, and band. `st_as_sf()` supports `what = "env"`,
#'   `"scores"`, and `"contributions"` when `output = "full"` or `"localenv"`
#'   has retained geometry or coordinates.
#'   `what = "diagnostics"` returns one row per variable/band with `status` and
#'   any kNN population/radius diagnostics. `what = "interpolation"` returns
#'   one row per stored interval neighborhood/band with `location`,
#'   `bottom_category_share`, `top_category_share`, and `interpolated_share`.
#'   These are weighted local population shares in the first category, last
#'   category, and categories split across requested groups. They are not
#'   uncertainty estimates; allocation fractions remain in
#'   `definition$summary$groups$allocation`. Non-interval analyses contribute
#'   no interpolation rows. Interpolation diagnostics require retained environments.
#'
#'   Network local and interpolation tables also identify `node`. Spatial
#'   extraction returns the compact local table with corresponding geometry.
#'   Network-only results require stored geometry or coordinates for mapping.
#' @name seg_result
NULL

#' @rdname seg_result
#' @method print seg_result
#' @export
print.seg_result <- function(x, ...) {
  if (!is.null(x$cases)) {
    cat("\nSpatial segregation across", length(x$cases), "cases and",
        length(x$cases[[1L]]$analyses), "variables\n")
    cat("Neighborhoods:", x$neighbors$type, "| search:", x$neighbors$engine, "\n")
    cat("Bands:", paste(x$bands, collapse = ", "), "\n")
    cat("Cases: ", paste(utils::head(as.character(x$case_ids), 6), collapse = ", "),
        if (length(x$cases) > 6) ", ..." else "", "\n", sep = "")
    return(invisible(x))
  }
  cat("\nSpatial segregation across", length(x$analyses), "variables\n")
  cat("Neighborhoods:", x$neighbors$type, "| bands:", paste(x$bands, collapse = ", "), "\n")
  cat("Space:", x$neighbors$space, "| search:", x$neighbors$engine, "\n")
  for (name in names(x$analyses)) {
    a <- x$analyses[[name]]
    cat(" ", name, ": ", a$definition$type, " (", a$definition$method, ")\n", sep = "")
    cat("    Scope:", a$settings$scope, "| output:", a$settings$output, "\n")
    if (a$definition$type == "continuous") {
      interpretation <- if (a$definition$localenv$type == "quantiles") {
        "focal-group averages of local global-rank statistics, not exposure probabilities"
      } else "exposure/isolation probabilities"
      cat("    Summary:", a$definition$localenv$type, "-", interpretation, "\n")
    }
    bad <- unique(a$status[a$status != "ok"])
    if (length(bad)) cat("    Status:", paste(bad, collapse = ", "), "\n")
  }
  invisible(x)
}

#' @rdname seg_result
#' @method as.data.frame seg_result
#' @export
as.data.frame.seg_result <- function(x, row.names = NULL, optional = FALSE,
    ..., what = c("indices", "pairwise", "exposure", "env", "diagnostics",
                 "interpolation", "scores", "contributions"), variables = NULL, bands = NULL) {
  .seg_unused_arguments(list(...))
  what <- match.arg(what)
  if (!is.null(x$cases)) {
    blocks <- lapply(x$cases, function(one) as.data.frame(one, what = what,
      variables = variables, bands = bands, optional = optional))
    out <- .seg_bind_result_rows(blocks)
    rownames(out) <- row.names
    return(out)
  }
  extract <- switch(what, indices = .seg_indices_rows,
    pairwise = .seg_pairwise_rows, exposure = .seg_exposure_rows,
    env = .seg_env_rows, diagnostics = .seg_diagnostic_rows,
    interpolation = .seg_interpolation_rows,
    scores = function(...) .seg_score_rows(..., what = "scores"),
    contributions = function(...) .seg_score_rows(..., what = "contributions"))
  selected_variables <- if (is.null(variables)) names(x$analyses) else variables
  if (!is.character(selected_variables) || !length(selected_variables) ||
      anyNA(selected_variables) || !all(selected_variables %in% names(x$analyses)))
    stop("'variables' must name stored analyses", call. = FALSE)
  selected_bands <- if (is.null(bands)) seq_along(x$bands) else match(bands, x$bands)
  if (!length(selected_bands) || anyNA(selected_bands))
    stop("'bands' must match stored bands", call. = FALSE)
  if (is.null(variables)) {
    available <- vapply(x$analyses[selected_variables], function(a) switch(what,
      indices = , pairwise = , exposure = a$settings$output != "localenv",
      env = !is.null(a$env),
      scores = !is.null(a$scores), contributions = !is.null(a$contributions),
      interpolation = !identical(a$definition$representation, "intervals") ||
        !is.null(a$env), TRUE), logical(1))
    if (any(available)) selected_variables <- selected_variables[available]
  }
  analyses <- x$analyses[selected_variables]
  schema <- .seg_result_schema(analyses, what)
  rows <- list()
  keys <- c(d = "dissimilarity", h = "information", r = "diversity", s = "square_root")
  for (name in selected_variables) {
    analysis <- x$analyses[[name]]
    if (what %in% c("scores", "contributions") && is.null(analysis[[what]])) {
      stop("Local ", what, " were not retained; set keep_", what, " = TRUE",
           call. = FALSE)
    }
    if (what %in% c("indices", "pairwise", "exposure") &&
        analysis$settings$output == "localenv") {
      stop("Indices were not calculated; use output = 'indices' or 'full'", call. = FALSE)
    }
    if (what == "env" && is.null(analysis$env)) {
      stop("Local environments were not stored; use output = 'full' or 'localenv'",
           call. = FALSE)
    }
    if (what == "interpolation" &&
        identical(analysis$definition$representation, "intervals") &&
        is.null(analysis$env)) {
      stop("Interpolation diagnostics were not stored; use output = 'full' or 'localenv'",
           call. = FALSE)
    }
    for (b in selected_bands) {
      base <- list(variable = name, band = x$bands[b])
      rows <- c(rows, extract(analysis, base, b, keys))
    }
  }
  out <- .seg_bind_result_rows(rows, schema)
  if (what %in% c("env", "interpolation", "scores", "contributions") && !is.null(x$node_id)) {
    out$node <- x$node_id[out$location]
  }
  if (!is.null(x$case)) {
    out <- data.frame(case = rep(x$case, nrow(out)), out, check.names = FALSE)
  }
  rownames(out) <- row.names
  out
}

# Binds result blocks for as.data.frame() at both case and variable/band levels.
# Typed schemas preserve empty tables, ID classes, and structural missing fields.
.seg_bind_result_rows <- function(rows, schema = NULL) {
  if (is.null(schema)) {
    columns <- unique(unlist(lapply(rows, names), use.names = FALSE))
    schema <- as.data.frame(setNames(lapply(columns, function(column) {
      first <- which(vapply(rows, function(z) column %in% names(z), logical(1)))[1L]
      rows[[first]][[column]][FALSE]
    }), columns), stringsAsFactors = FALSE)
  }
  if (!length(rows)) return(schema)
  rows <- lapply(rows, function(block) {
    for (column in setdiff(names(schema), names(block)))
      block[[column]] <- schema[[column]][rep(NA_integer_, nrow(block))]
    block[names(schema)]
  })
  do.call(rbind, rows)
}

# Defines as.data.frame() columns from the selected analyses so empty tables
# and mixed share/rank outputs use the same types and column order.
.seg_result_schema <- function(analyses, what) {
  base <- list(variable = character(), band = numeric())
  schema <- switch(what,
    indices = c(base, list(measure = character(), value = numeric())),
    pairwise = c(base, list(measure = character(), group_a = character(),
      group_b = character(), value = numeric())),
    exposure = {
      ranks <- vapply(analyses, function(a) {
        a$definition$type == "continuous" && a$definition$localenv$type == "quantiles"
      }, logical(1))
      c(base, list(reference = character()),
        if (any(!ranks)) list(target = character()),
        if (any(ranks)) list(statistic = character(), prob = numeric()),
        list(value = numeric()))
    },
    env = {
      ranks <- vapply(analyses, function(a) a$definition$localenv$type == "quantiles",
                      logical(1))
      c(base, list(location = integer()),
        if (any(!ranks)) list(category = character()),
        if (any(ranks)) list(statistic = character(), prob = numeric()),
        list(value = numeric()))
    },
    diagnostics = {
      knn <- unique(unlist(lapply(analyses, function(a) setdiff(names(a$knn), "band"))))
      c(base, list(status = character()), setNames(rep(list(numeric()), length(knn)), knn))
    },
    interpolation = c(base, list(location = integer(),
      bottom_category_share = numeric(), top_category_share = numeric(),
      interpolated_share = numeric())),
    scores = , contributions = {
      has_pairs <- any(vapply(analyses, function(a)
        length(a[[what]]$pairwise) > 0, logical(1)))
      c(base, list(scope = character(), location = integer(), measure = character()),
        if (has_pairs) list(group_a = character(), group_b = character()),
        list(value = numeric()))
    }
  )
  as.data.frame(schema, stringsAsFactors = FALSE)
}

# Flattens the selected overall scalar measures for one band in as.data.frame().
.seg_indices_rows <- function(analysis, base, b, keys) {
  lapply(intersect(names(keys), names(analysis$indices$overall)), function(key) {
    data.frame(base, measure = unname(keys[key]),
      value = unname(analysis$indices$overall[[key]][b]), row.names = NULL)
  })
}

# Serves both local scores and contributions in as.data.frame(), expanding
# location and pair identifiers only where applicable.
.seg_score_rows <- function(analysis, base, b, keys, what) {
  rows <- list()
  for (scope in c("overall", "pairwise")) {
    store <- analysis[[what]][[scope]]
    for (key in intersect(names(keys), names(store))) {
      values <- store[[key]][[b]]
      if (scope == "overall") {
        rows[[length(rows) + 1L]] <- data.frame(base, scope = scope,
          location = seq_along(values), measure = unname(keys[key]), value = values)
      } else if (ncol(values)) {
        pair <- rep(seq_len(ncol(values)), each = nrow(values))
        rows[[length(rows) + 1L]] <- data.frame(base, scope = scope,
          location = rep(seq_len(nrow(values)), times = ncol(values)),
          measure = unname(keys[key]), group_a = analysis$score_pairs$group_a[pair],
          group_b = analysis$score_pairs$group_b[pair], value = as.vector(values))
      }
    }
  }
  rows
}

# Extracts one triangle of each symmetric pairwise matrix for as.data.frame(),
# so each category pair appears once.
.seg_pairwise_rows <- function(analysis, base, b, keys) {
  rows <- list()
  for (key in intersect(names(keys), names(analysis$indices$pairwise))) {
    mat <- analysis$indices$pairwise[[key]][[b]]
    index <- which(upper.tri(mat), arr.ind = TRUE)
    if (nrow(index)) {
      rows[[length(rows) + 1L]] <- data.frame(base,
        measure = unname(keys[key]), group_a = rownames(mat)[index[, 1]],
        group_b = colnames(mat)[index[, 2]], value = mat[index], row.names = NULL)
    }
  }
  rows
}

# Converts exposure matrices or rank profiles to as.data.frame() rows while
# preserving their different identifiers and meanings.
.seg_exposure_rows <- function(analysis, base, b, keys) {
  if (length(analysis$rank_summary)) {
    mat <- analysis$rank_summary[[b]]
    if (!length(mat)) return(list())
    stats <- analysis$definition$localenv$statistics
    column <- rep(seq_len(ncol(mat)), each = nrow(mat))
    return(list(data.frame(base,
      reference = rep(rownames(mat), times = ncol(mat)),
      statistic = stats$statistic[column], prob = stats$prob[column],
      value = as.vector(mat))))
  }
  mat <- analysis$indices$overall$p[[b]]
  if (is.null(mat) || !length(mat)) return(list())
  list(data.frame(base,
    expand.grid(reference = rownames(mat), target = colnames(mat),
                stringsAsFactors = FALSE), statistic = "share", value = as.vector(mat)))
}

# Expands a band matrix into location-level as.data.frame() rows. Retains
# the location index so st_as_sf() can attach the matching geometry.
.seg_env_rows <- function(analysis, base, b, keys) {
  env <- analysis$env[[b]]
  if (!length(env)) return(list())
  definition <- analysis$definition$localenv
  location <- rep(seq_len(nrow(env)), times = ncol(env))
  column <- rep(seq_len(ncol(env)), each = nrow(env))
  if (definition$type == "quantiles") {
    stats <- definition$statistics
    return(list(data.frame(base, location = location,
      statistic = stats$statistic[column], prob = stats$prob[column],
      value = as.vector(env))))
  }
  list(data.frame(base, location = location, category = colnames(env)[column],
    statistic = "share", value = as.vector(env)))
}

# Builds one as.data.frame() row per analysis and band, adding kNN diagnostics
# when stored.
.seg_diagnostic_rows <- function(analysis, base, b, keys) {
  block <- data.frame(base, status = unname(analysis$status[b]))
  if (!is.null(analysis$knn)) {
    block <- cbind(block, analysis$knn[b, setdiff(names(analysis$knn), "band"), drop = FALSE])
  }
  list(block)
}

# Extracts stored interval interpolation diagnostics as location-level rows
# for as.data.frame(); other analyses contribute no rows.
.seg_interpolation_rows <- function(analysis, base, b, keys) {
  if (!length(analysis$interpolation)) return(list())
  mat <- analysis$interpolation[[b]]
  if (!nrow(mat)) return(list())
  list(data.frame(base, location = seq_len(nrow(mat)), mat))
}

#' @rdname seg_result
#' @method st_as_sf seg_result
#' @export
st_as_sf.seg_result <- function(x, ..., variables = NULL, bands = NULL,
                                what = c("env", "scores", "contributions")) {
  .seg_unused_arguments(list(...))
  if (is.null(x$coords) && is.null(x$geometry))
    stop("No stored geometry or coordinates; network-only results cannot be mapped",
         call. = FALSE)
  what <- match.arg(what)
  out <- as.data.frame(x, what = what, variables = variables, bands = bands)
  if (!nrow(out)) stop("No stored local result rows to extract", call. = FALSE)
  if (!is.null(x$geometry)) sf::st_sf(out, geometry = x$geometry[out$location]) else {
    out$x <- x$coords[out$location, 1]
    out$y <- x$coords[out$location, 2]
    sf::st_as_sf(out, coords = c("x", "y"), crs = x$crs)
  }
}

#' @rdname seg_result
#' @method as.list seg_result
#' @export
as.list.seg_result <- function(x, ...) unclass(x)

# Shares the strict dots check between as.data.frame() and st_as_sf() so
# misspelled extraction options produce an error.
.seg_unused_arguments <- function(dots) {
  if (length(dots)) stop("Unused extraction arguments: ",
                         paste(names(dots), collapse = ", "), call. = FALSE)
}
