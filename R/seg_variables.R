#' Describe variables for a shared spatial segregation calculation
#'
#' These constructors describe data columns; they do not search neighborhoods.
#' Select columns with `vars`: one column supplies observed memberships or values;
#' multiple columns supply counts in category/value order. Column selectors are
#' names or integer positions. A single count column is not supported.
#' Continuous variables accept exact observations or interval membership/counts
#' through `representation`. Interval CDFs are interpolated directly; individual
#' values are not imputed.
#' Ordinal overall indices combine cumulative category thresholds. With
#' `scope = "both"` or `"pairwise"`, ordinal categories also receive pairwise
#' indices using the categorical two-group formulas within the same neighborhoods.
#' Ordinal exposure/isolation uses the original categories. Pairwise indices and
#' exposures use full category shares even when retained environments are quantiles.
#'
#' @param vars One column containing category membership or observed values, or
#'   two or more columns containing nonnegative population counts. Count columns
#'   must be supplied in category/value order.
#' @param levels Category order for observed ordinal or interval membership.
#'   Ordered factors supply their own order. For categorical values, factor
#'   levels or sorted distinct values are used when this is omitted.
#' @param measures Optional measure selection overriding [seg_measure()].
#' @param scope Optional `"overall"`, `"pairwise"`, or `"both"`, overriding
#'   the global scope in [seg_measure()]. `NULL` inherits the global setting.
#' @param output Optional `"indices"`, `"full"`, or `"localenv"`, overriding
#'   the global return mode. `NULL` inherits the global setting.
#' @param keep_scores,keep_contributions Optional logical retention settings
#'   overriding [seg_measure()]. `NULL` inherits the global setting; `FALSE` explicitly
#'   disables retention for this variable. Effective settings are recorded in
#'   `result$analyses[[name]]$settings`.
#' @param weights Population-weight column or numeric vector when `vars` selects
#'   one column. This supports survey weights and long-form aggregated counts,
#'   where each row records a location/category combination. Rows at the same
#'   location share coordinates or geometry. With multiple `vars`, the count
#'   columns already supply population mass and `weights` must be omitted.
#'   Defaults to the `population` supplied to [seg_measure()], or its inferred default.
#' @param na.action Missing-value policy: `"fail"` or `"omit"`. Omission applies
#'   only to this variable, including its kNN stopping population. Other
#'   variables and radius bandwidths are unaffected.
#' @param representation Continuous input representation. The default `"auto"`
#'   selects `"exact"` for numeric observations and for count columns with
#'   `scores`. Count columns without scores, ordered-factor observations, and
#'   observations with explicit `levels` select `"intervals"`. Unordered
#'   memberships require an explicit order through `levels` or an ordered factor.
#'   Select `"exact"` or `"intervals"` to specify the interpretation directly.
#'   Resolved representation is recorded in the analysis definition. Interval
#'   membership requires an ordered factor or explicit `levels`; count columns
#'   supply interval order.
#' @param interpolation For `representation = "intervals"`, CDF interpolation.
#'   `NULL` selects `"linear"` in overall percentile space, not in the original
#'   measurement units. Exact input does not accept this option.
#' @param scores Exact numerical values corresponding to continuous count
#'   columns with `representation = "exact"`. These are observed values, not
#'   interval midpoints. Supplying scores selects exact-value counts in automatic
#'   mode. Interval input does not accept scores.
#' @param summary Continuous summary specification, defaulting to shares of the
#'   supplied categories for interval input and aggregate counts, or global
#'   quartile groups for numeric microdata. Explicit `probs` requests quantile
#'   groups for either input. Shares use the same
#'   global groups for retained local shares and group-averaged exposure.
#'   `type = "quantiles"` returns local global-rank quantiles and their
#'   focal-group averages instead, evaluated at quartiles by default. With
#'   omitted `probs`, these averages and pairwise indices retain the supplied
#'   categories for interval input and aggregate counts. `mean_rank = TRUE`
#'   additionally reports the mean of the full local rank distribution. `probs` is one strictly increasing
#'   vector in [0, 1]; optional 0/1 delimiters define tails but are not evaluated
#'   as local extrema. Exact shares may use finite absolute `cuts` instead of
#'   `probs`, with implicit infinite tails. Cuts require exact values and shares
#'   mode. Tied values remain together, so achieved group sizes may differ.
#'   Probability-defined groups use percentile labels such as `"P0-P25"` and
#'   `"P75-P100"`. For exact data these reflect achieved population shares after
#'   preserving ties. Explicit `cuts` use value-based labels. Resolved value
#'   cutoffs and labels remain in `definition$summary$groups$breaks` and
#'   `value_labels` for exact data; `boundaries` records the percentile boundaries
#'   on the [0, 1] scale for both representations.
#'   Default category groups retain their input labels. When `probs` is supplied,
#'   interval input retains the requested grid. Categories crossing requested
#'   boundaries are proportionally apportioned in global-percentile space, using
#'   the same fractions for original focal populations and local target shares.
#'   This assumes the same within-category rank distribution at every location;
#'   it does not require monetary tail bounds. Allocation is resolved once before
#'   surface construction and recorded in `definition$summary$groups$allocation`.
#'   One warning identifies when interpolation is required. Extraction with
#'   `what = "interpolation"` reports `bottom_category_share`, `top_category_share`, and `interpolated_share`
#'   rather than printing numerical details in the warning. The last is the local
#'   population share in input categories split across requested groups, not an
#'   uncertainty estimate. Choose probabilities on actual cumulative interval
#'   boundaries to avoid apportionment. Summary choices do not change continuous
#'   scalar indices.
#' @param localenv Ordinal retained-environment summary: `"shares"` (default)
#'   or `"quantiles"`, optionally supplied as a list with `type` and `probs`.
#'   Quantiles default to `seq(.1, .9, .1)` and report weighted empirical
#'   quantiles of global population midranks. Probabilities increase strictly
#'   within [0, 1]; endpoints select local support limits. Shares default to
#'   supplied categories; optional probabilities must run from 0 to 1 and group
#'   categories by global population quantiles, keeping ties together.
#'   These options affect retained environments only, not indices or exposure.
#'   Continuous variables use `summary` instead.
#' @return A `seg_variable` specification for `seg_measure(variables = ...)`.
#' @seealso [seg_measure()]
#' @examples
#' seg_categorical(vars = "ethnicity")
#' seg_ordinal(vars = c("primary", "secondary", "tertiary"))
#' seg_continuous(vars = c("low", "middle", "high"))
#' seg_continuous(vars = "income", summary = list(probs = c(.2, .8)))
#' seg_continuous(vars = "income",
#'                  summary = list(type = "quantiles", mean_rank = TRUE))
#' seg_ordinal(vars = "education", localenv = "quantiles")
#' seg_categorical(vars = c("group_a", "group_b"),
#'                   scope = "both", output = "full")
#' # Long-form counts: one membership column and a population weight per row.
#' seg_categorical(vars = "race", weights = "population_count")
#' @name seg_variables
NULL

# Builds the common descriptor for seg_categorical(), seg_ordinal(), and
# seg_continuous(), so column selection and overrides follow the same rules.
.seg_variable <- function(type, vars, levels, measures, weights,
                          na.action, interpolation = NULL, scores = NULL,
                          localenv = "shares", representation = NULL,
                          scope = NULL, output = NULL, keep_scores = NULL,
                          keep_contributions = NULL) {
  if (!typeof(vars) %in% c("character", "integer", "double") ||
      !is.null(dim(vars)) ||
      !length(vars) || anyNA(vars) || anyDuplicated(vars) ||
      (is.character(vars) && any(!nzchar(vars))) ||
      (is.numeric(vars) && any(!is.finite(vars) | vars < 1 | vars != floor(vars)))) {
    stop("'vars' must select non-empty, non-missing, unique column names or positive integer positions",
         call. = FALSE)
  }
  values <- if (length(vars) == 1L) vars else NULL
  counts <- if (length(vars) > 1L) vars else NULL
  if (!is.null(counts) && !is.null(weights)) {
    stop("'weights' applies only when 'vars' selects one column; count columns already contain mass",
         call. = FALSE)
  }
  if (!is.null(scope)) scope <- match.arg(scope, c("both", "overall", "pairwise"))
  if (!is.null(output)) output <- match.arg(output, c("indices", "full", "localenv"))
  .seg_validate_retention(keep_scores, keep_contributions, allow_null = TRUE)
  structure(list(type = type, values = values, counts = counts,
                 levels = levels, measures = measures, weights = weights,
                 scope = scope, output = output, keep_scores = keep_scores,
                 keep_contributions = keep_contributions,
                 na.action = match.arg(na.action, c("fail", "omit")),
                 interpolation = interpolation, scores = scores,
                 representation = representation,
                 localenv = if (type == "continuous") NULL else
                   .seg_localenv_options(localenv)),
            class = "seg_variable")
}

#' @rdname seg_variables
#' @export
seg_categorical <- function(vars, levels = NULL,
                            measures = NULL, weights = NULL,
                            na.action = c("fail", "omit"), scope = NULL,
                            output = NULL, keep_scores = NULL,
                            keep_contributions = NULL) {
  .seg_variable("categorical", vars, levels, measures, weights,
                na.action, scope = scope, output = output,
                keep_scores = keep_scores, keep_contributions = keep_contributions)
}

#' @rdname seg_variables
#' @export
seg_ordinal <- function(vars, levels = NULL,
                       measures = NULL, weights = NULL,
                       na.action = c("fail", "omit"), localenv = "shares",
                       scope = NULL, output = NULL, keep_scores = NULL,
                       keep_contributions = NULL) {
  .seg_variable("ordinal", vars, levels, measures, weights,
                na.action, localenv = localenv, scope = scope,
                output = output, keep_scores = keep_scores,
                keep_contributions = keep_contributions)
}

#' @rdname seg_variables
#' @export
seg_continuous <- function(vars,
                          representation = c("auto", "exact", "intervals"),
                          levels = NULL, scores = NULL,
                          measures = NULL, weights = NULL,
                          na.action = c("fail", "omit"), interpolation = NULL,
                          summary = NULL, scope = NULL, output = NULL,
                          keep_scores = NULL, keep_contributions = NULL) {
  representation <- match.arg(representation)
  .seg_validate_representation(representation, levels, scores, interpolation)
  if (!is.null(interpolation) && !identical(interpolation, "linear")) {
    stop("Only 'linear' interpolation is supported", call. = FALSE)
  }
  if (representation == "intervals") interpolation <- .seg_or_default(interpolation, "linear")
  if (length(vars) == 1L && !is.null(scores)) {
    stop("'scores' applies only to continuous count columns", call. = FALSE)
  }
  out <- .seg_variable("continuous", vars, levels, measures, weights,
                na.action, scores = scores,
                representation = representation,
                interpolation = interpolation, scope = scope, output = output,
                keep_scores = keep_scores, keep_contributions = keep_contributions)
  out$summary <- .seg_summary_options(summary)
  out
}

# Shares representation-specific rules between seg_continuous() and
# seg_measure(), which can resolve 'auto' only once it has the data.
.seg_validate_representation <- function(representation, levels, scores, interpolation) {
  if (representation == "intervals" && !is.null(scores))
    stop("'scores' applies only to representation = 'exact'", call. = FALSE)
  if (representation == "exact") {
    if (!is.null(levels))
      stop("Continuous 'levels' applies only to representation = 'intervals'", call. = FALSE)
    if (!is.null(interpolation))
      stop("'interpolation' applies only to representation = 'intervals'", call. = FALSE)
  }
}

# Prepares count columns for seg_measure() once so population inference and
# compilation can share the matrix and missing-count policy.
.seg_counts <- function(spec, data) {
  out <- .seg_select(data, spec$counts)
  if (!all(vapply(out, is.numeric, logical(1))))
    stop("Count columns must be numeric", call. = FALSE)
  out <- as.matrix(out)
  bad <- !stats::complete.cases(out)
  if (any(bad)) {
    if (spec$na.action == "fail") stop("Missing counts; choose na.action = 'omit' explicitly",
                                       call. = FALSE)
    out[bad, ] <- 0
  }
  if (any(!is.finite(out) | out < 0))
    stop("Counts must be finite and nonnegative", call. = FALSE)
  storage.mode(out) <- "double"
  out
}

# Builds a native variable descriptor for seg_measure() from counts or observed
# values. Weights, measures, and count matrices arrive validated by its preparation
# stage; surface analyses reuse it with gridded counts and resolved groups.
.seg_compile <- function(spec, data, population, resolution = NULL,
                           variable = NULL, counts = NULL) {
  n <- nrow(data)
  if (spec$type == "continuous" && spec$representation == "auto") {
    if (!is.null(spec$counts)) {
      spec$representation <- if (is.null(spec$scores)) "intervals" else "exact"
    } else {
      values <- .seg_select(data, spec$values)[[1L]]
      spec$representation <- if (!is.null(spec$levels) || is.ordered(values)) {
        "intervals"
      } else if (is.numeric(values)) {
        "exact"
      } else {
        stop("Continuous memberships require an ordered factor or explicit 'levels'; ",
             "use numeric observations for exact values", call. = FALSE)
      }
    }
    .seg_validate_representation(spec$representation, spec$levels, spec$scores,
                                   spec$interpolation)
    if (spec$representation == "intervals")
      spec$interpolation <- .seg_or_default(spec$interpolation, "linear")
  }
  intervals <- identical(spec$representation, "intervals")
  exact <- spec$type == "continuous" && !intervals
  selected <- spec$measures
  code <- rep(-1L, n)
  scores <- numeric()
  if (!is.null(spec$counts)) {
    labels <- colnames(counts)
    mass <- rowSums(counts)
    if (exact) {
      scores <- spec$scores
      if (!is.numeric(scores) || length(scores) != ncol(counts) ||
          any(!is.finite(scores)))
        stop("Supply finite 'scores' matching the continuous count columns", call. = FALSE)
      if (anyDuplicated(scores)) {
        # Collapse duplicate exact values, retaining the population mass.
        keys <- sort(unique(scores))
        counts <- vapply(keys, function(z) rowSums(counts[, scores == z, drop = FALSE]),
                         numeric(n))
        dim(counts) <- c(n, length(keys))
        scores <- keys
      } else {
        order <- order(scores)
        scores <- scores[order]
        counts <- counts[, order, drop = FALSE]
      }
      labels <- as.character(scores)
    }
    totals <- colSums(counts)
  } else {
    counts <- matrix(numeric(), 0L, 0L)
    values <- .seg_select(data, spec$values)[[1L]]
    mass <- .seg_or_default(spec$weights, population)
    bad <- is.na(values)
    if (any(bad) && spec$na.action == "fail")
      stop("Missing observed values; choose na.action = 'omit' explicitly", call. = FALSE)
    mass[bad] <- 0
    if (exact) {
      if (!is.numeric(values) || any(!is.finite(values[!bad])))
        stop("Continuous values must be finite numeric observations", call. = FALSE)
      scores <- sort(unique(values[mass > 0]))
      labels <- as.character(scores)
      code <- match(values, scores) - 1L
    } else {
      levels <- spec$levels
      if (is.null(levels)) {
        if (spec$type != "categorical" && !is.ordered(values))
          stop("Supply ordered 'levels' or an ordered factor", call. = FALSE)
        levels <- if (is.factor(values)) levels(values) else sort(unique(values[!bad]))
      }
      if (anyNA(levels) || anyDuplicated(levels) || !length(levels))
        stop("Category levels must be nonempty, distinct, and nonmissing", call. = FALSE)
      labels <- as.character(levels)
      code <- match(values, levels) - 1L
      if (any(is.na(code) & !bad))
        stop("Some observed categories are absent from 'levels'", call. = FALSE)
    }
    code[is.na(code)] <- -1L
    totals <- numeric(length(labels))
    ok <- code >= 0L & mass > 0
    if (any(ok)) {
      sums <- rowsum(mass[ok], code[ok] + 1L, reorder = FALSE)
      totals[as.integer(rownames(sums))] <- sums[, 1L]
    }
  }
  total <- sum(mass)
  if (!is.finite(total) || total <= 0)
    stop("Each variable must have a positive finite eligible population", call. = FALSE)
  shares <- totals / total
  type <- switch(spec$type, categorical = 0L, ordinal = 1L,
                 continuous = if (intervals) 2L else 3L)
  method <- c("categorical", "ordinal_cutoffs", "linear_cdf", "empirical_midrank")[type + 1L]
  continuous <- spec$type == "continuous"
  resolved <- NULL
  exposure <- if (continuous) {
    summary_labels <- if (exact && !is.null(spec$counts)) {
      .seg_or_default(colnames(counts), labels)
    } else labels
    resolved <- .seg_or_default(resolution, .seg_resolve_summary(
      spec$summary, summary_labels, totals, scores, intervals, variable,
      source_groups = intervals || !is.null(spec$counts)))
    list(map = resolved$map, labels = resolved$labels,
         definition = resolved$definition$groups)
  } else list(map = seq_along(labels) - 1L, labels = labels, definition = NULL)
  localenv <- if (continuous) {
    definition <- resolved$definition
    list(type = if (definition$type == "shares") 0L else 1L,
      probs = .seg_or_default(definition$effective_probs, numeric()),
      map = resolved$map,
      definition = list(type = definition$type, probs = definition$effective_probs,
        labels = definition$statistics$label, statistics = definition$statistics,
        scale = if (definition$type == "quantiles") "global_rank" else "group_share"))
  } else .seg_localenv(spec, labels, totals)
  list(type = type, counts = counts, code = as.integer(code), mass = mass,
       summary_config = resolved,
       group_fractions = if (continuous) resolved$fractions else matrix(numeric(), 0, 0),
       split_intervals = if (continuous) {
         as.integer(resolved$definition$interpolation$split_intervals - 1L)
       } else integer(),
       env_mean = continuous && spec$summary$mean_rank,
       pair_labels = if (continuous) resolved$labels else labels,
       totals = totals, shares = shares, scores = scores,
       z = cumsum(shares) - shares / 2,
       exposure_map = as.integer(exposure$map), exposure_labels = exposure$labels,
       env_type = localenv$type, env_probs = localenv$probs,
       env_map = localenv$map, env_columns = length(localenv$definition$labels),
       flags = as.integer(c("exposure", "information", "diversity", "dissimilarity",
                            "square_root") %in% selected),
       metadata = list(type = spec$type, representation = spec$representation,
                       interpolation = spec$interpolation,
                       method = method, measures = selected,
                       eligible_population = total,
                       labels = labels, shares = shares, scores = scores,
                       ties = if (type == 3L) "weighted_midrank" else NULL,
                       na.action = spec$na.action,
                       summary = if (continuous) resolved$definition else NULL,
                       exposure = exposure$definition, localenv = localenv$definition))
}
