# Validates and normalizes seg_continuous() summary options before data
# are available. Records whether probabilities were explicitly supplied.
.seg_summary_options <- function(summary) {
  if (is.null(summary)) {
    summary <- list()
  }
  if (is.character(summary) && length(summary) == 1L) {
    summary <- list(type = summary)
  }
  if (!is.list(summary) ||
      (length(summary) && (is.null(names(summary)) || anyNA(names(summary)) ||
       any(!nzchar(names(summary))) || anyDuplicated(names(summary)) ||
       any(!names(summary) %in% c("type", "probs", "cuts", "mean_rank"))))) {
    stop("'summary' must contain type, probs or cuts, and optional mean_rank",
         call. = FALSE)
  }
  type <- .seg_or_default(summary$type, "shares")
  if (!is.character(type) || length(type) != 1L || is.na(type) ||
      !type %in% c("shares", "quantiles")) {
    stop("summary 'type' must be 'shares' or 'quantiles'", call. = FALSE)
  }
  if (!is.null(summary$probs) && !is.null(summary$cuts)) {
    stop("Supply only one of summary 'probs' and 'cuts'", call. = FALSE)
  }
  mean_rank <- .seg_or_default(summary$mean_rank, FALSE)
  if (!is.logical(mean_rank) || length(mean_rank) != 1L || is.na(mean_rank) ||
      (mean_rank && type != "quantiles")) {
    stop("summary 'mean_rank' must be TRUE or FALSE and requires quantiles mode",
         call. = FALSE)
  }
  if (!is.null(summary$cuts)) {
    cuts <- summary$cuts
    if (type != "shares") {
      stop("summary 'cuts' requires shares mode; use 'probs' for rank quantiles",
           call. = FALSE)
    }
    if (!is.numeric(cuts) || !length(cuts) || any(!is.finite(cuts)) ||
        is.unsorted(cuts, strictly = TRUE)) {
      stop("summary 'cuts' must be finite and strictly increasing", call. = FALSE)
    }
    return(list(type = type, probs = NULL, cuts = cuts, mean_rank = mean_rank))
  }
  probs <- .seg_or_default(summary$probs, c(.25, .5, .75))
  if (!is.numeric(probs) || !length(probs) || any(!is.finite(probs)) ||
      any(probs < 0 | probs > 1) || is.unsorted(probs, strictly = TRUE)) {
    stop("summary 'probs' must increase strictly within [0, 1]", call. = FALSE)
  }
  list(type = type, probs = as.numeric(probs[probs > 0 & probs < 1]),
       cuts = NULL, mean_rank = mean_rank,
       probs_explicit = !is.null(summary$probs))
}

# Resolves reporting groups for seg_measure() from the regional distribution.
# Keeps the rank cutoffs, interval allocations, and output labels aligned.
.seg_resolve_summary <- function(options, labels, totals, scores, intervals,
                                   variable = NULL, source_groups = intervals) {
  requested <- options$probs
  effective <- requested
  shares <- as.numeric(totals / sum(totals))
  cumulative <- cumsum(shares)
  reasons <- character()
  fractions <- matrix(numeric(), 0L, 0L)
  split_intervals <- integer()
  use_categories <- source_groups && !isTRUE(options$probs_explicit) &&
    is.null(options$cuts)
  if (use_categories) {
    requested <- NULL
    if (options$type == "shares") effective <- NULL
    group_labels <- labels
    map <- seq_along(labels) - 1L
    if (intervals) {
      fractions <- diag(length(labels))
      dimnames(fractions) <- list(labels, labels)
      map <- integer()
    }
    group_definition <- list(labels = labels, shares = shares,
      method = "input_categories")
    if (intervals) {
      group_definition$interval_ids <- as.list(seq_along(labels))
      group_definition$interval_labels <- labels
      group_definition$zero_mass_intervals <- which(totals == 0)
      group_definition$source_shares <- shares
      group_definition$allocation <- fractions
    }
  } else if (intervals) {
    if (!is.null(options$cuts)) {
      stop("Absolute summary 'cuts' require exact values, not interval labels",
           call. = FALSE)
    }
    lower <- c(0, utils::head(cumulative, -1L))
    edges <- c(0, requested, 1)
    groups <- seq_len(length(edges) - 1L)
    fractions <- matrix(0, length(labels), length(groups))
    # Allocate population mass, not income values; open monetary tails need no bounds.
    for (j in groups) {
      overlap <- pmax(0, pmin(cumulative, edges[j + 1L]) - pmax(lower, edges[j]))
      fractions[shares > 0, j] <- overlap[shares > 0] / shares[shares > 0]
    }
    positive <- which(shares > 0)
    allocated <- rowSums(fractions)
    representable <- positive[allocated[positive] > 0]
    fractions[representable, ] <- fractions[representable, , drop = FALSE] /
      allocated[representable]
    # Tiny positive tails can collapse to one rank at floating-point precision.
    collapsed <- positive[allocated[positive] == 0]
    if (length(collapsed)) {
      containing_group <- findInterval(cumulative[collapsed], edges, all.inside = TRUE)
      fractions[cbind(collapsed, containing_group)] <- 1
    }
    split_intervals <- which(rowSums(fractions > 8 * .Machine$double.eps) > 1L)
    if (length(split_intervals)) {
      reasons <- "linear_interpolation"
    }
    map <- integer()
    group_definition <- list(shares = as.numeric(colSums(fractions * shares)),
      interval_ids = lapply(groups, function(j) which(fractions[, j] > 0)),
      interval_labels = labels, zero_mass_intervals = which(totals == 0),
      source_shares = shares,
      method = "linear_global_rank")
  } else {
    exposure <- .seg_rank_groups(scores, totals,
      if (is.null(options$cuts)) list(probs = c(0, requested, 1)) else
        list(breaks = c(-Inf, options$cuts, Inf)))
    map <- exposure$map
    group_labels <- exposure$labels
    group_definition <- exposure$definition
  }
  group_definition$requested_groups <- if (use_categories) {
    length(labels)
  } else if (is.null(options$cuts)) {
    length(requested) + 1L
  } else length(options$cuts) + 1L
  group_definition$effective_groups <- length(group_definition$shares)
  group_definition$boundaries <- c(0,
    utils::head(cumsum(group_definition$shares), -1L), 1)
  group_definition$label_scale <- if (use_categories) "category" else
    if (is.null(options$cuts)) "percentile" else "value"
  if (!intervals) group_definition$value_labels <- group_labels
  if (!use_categories && is.null(options$cuts)) {
    percentiles <- vapply(100 * group_definition$boundaries, function(p)
      format(p, digits = 12, trim = TRUE, scientific = FALSE), character(1))
    group_labels <- paste0("P", utils::head(percentiles, -1L), "-P",
                           utils::tail(percentiles, -1L))
    group_definition$labels <- group_labels
    if (intervals) {
      dimnames(fractions) <- list(labels, group_labels)
      group_definition$allocation <- fractions
    }
  }
  statistics <- if (options$type == "shares") {
    data.frame(statistic = "share", prob = NA_real_, label = group_labels)
  } else {
    data.frame(statistic = c(rep("quantile_rank", length(effective)),
                            if (options$mean_rank) "mean_rank" else character()),
               prob = c(effective, if (options$mean_rank) NA_real_),
               label = c(if (length(effective)) paste0("q", effective) else character(),
                         if (options$mean_rank) "mean_rank"))
  }
  definition <- list(type = options$type, requested_probs = requested,
    effective_probs = effective, cuts = options$cuts, mean_rank = options$mean_rank,
    groups = group_definition, reasons = reasons,
    interpolation = list(needed = length(split_intervals) > 0L,
                         split_intervals = split_intervals),
    status = if (length(split_intervals)) "interpolated" else "ok",
    statistics = statistics)
  if (length(reasons)) {
    warning("Interval summary interpolation",
      if (!is.null(variable)) paste0(" for '", variable, "'"),
      ": requested groups proportionally split income categories in global-rank space. ",
      "The same fractions apply at every location; see what = 'interpolation' ",
      "and definition$summary$groups$allocation.", call. = FALSE)
  }
  list(map = as.integer(map), fractions = fractions, labels = group_labels,
       definition = definition)
}
# Discrete and continuous reporting share this module; their separate resolvers
# preserve whole-category ties and interval interpolation respectively.

# Normalizes ordinal summary shorthand for seg_ordinal() through the
# shared constructor, so compilation receives validated options.
.seg_localenv_options <- function(options) {
  if (is.character(options) && length(options) == 1L && !is.na(options))
    options <- list(type = options)
  if (!is.list(options) || is.null(names(options)) || anyDuplicated(names(options)) ||
      any(!names(options) %in% c("type", "probs")) ||
      !is.character(options$type) || length(options$type) != 1L ||
      is.na(options$type) || !options$type %in% c("quantiles", "shares"))
    stop("'localenv' must select 'quantiles' or 'shares', optionally with probs",
         call. = FALSE)
  if (options$type == "quantiles") {
    options$probs <- .seg_or_default(options$probs, seq(.1, .9, .1))
  }
  p <- options$probs
  if (!is.null(p) && (!is.numeric(p) || !length(p) || any(!is.finite(p)) ||
      any(p < 0 | p > 1) || is.unsorted(p, strictly = TRUE)))
    stop("Local-environment probabilities must increase strictly within [0, 1]",
         call. = FALSE)
  if (options$type == "shares" && !is.null(p) &&
      (length(p) < 2L || p[1L] != 0 || p[length(p)] != 1))
    stop("Local share probabilities must increase from 0 to 1", call. = FALSE)
  options
}

# Builds discrete local-environment settings for seg_measure(), keeping reporting
# groups separate from the categories used for indices and exposure.
.seg_localenv <- function(spec, labels, totals) {
  options <- spec$localenv
  if (options$type == "quantiles") {
    statistics <- data.frame(
      statistic = "quantile_rank", prob = options$probs,
      label = paste0("q", format(options$probs, trim = TRUE))
    )
    return(list(type = 1L, probs = as.numeric(options$probs), map = integer(),
      definition = list(type = "quantiles", probs = options$probs,
        labels = statistics$label, statistics = statistics,
        scale = "global_rank", method = "inverse_midrank_cdf")))
  }
  if (!is.null(options$probs)) {
    groups <- .seg_rank_groups(seq_along(labels), totals,
                                 list(probs = options$probs))
    definition <- c(list(type = "shares"), groups$definition)
    map <- groups$map
  } else {
    definition <- list(type = "shares", labels = labels, method = "supplied_categories")
    map <- seq_along(labels) - 1L
  }
  list(type = 0L, probs = numeric(), map = as.integer(map), definition = definition)
}

# Shared by seg_measure()'s ordinal and exact continuous reporting paths to group
# ordered values without splitting ties. Takes validated probabilities
# spanning [0, 1] or breaks padded to cover the observed values.
.seg_rank_groups <- function(scores, totals, options) {
  present <- totals > 0
  v <- scores[present]
  mass <- totals[present]
  if (!length(v)) return(list(map = rep(-1L, length(scores)), labels = character(),
                              definition = list(requested = options)))
  if (!is.null(options$probs)) {
    p <- options$probs
    cumulative <- cumsum(mass) / sum(mass)
    inner <- v[vapply(p[-c(1L, length(p))],
                       function(q) which(cumulative >= q - 8*.Machine$double.eps)[1L], integer(1))]
    breaks <- c(-Inf, unique(inner[inner < max(v)]), Inf)
  } else {
    breaks <- options$breaks
  }
  groups <- cut(scores, breaks, include.lowest = TRUE, right = TRUE)
  original <- as.integer(groups)
  shares <- vapply(seq_len(length(breaks) - 1L),
                   function(k) sum(totals[which(original == k)]), numeric(1))
  keep <- which(shares > 0)
  list(map = match(original, keep) - 1L, labels = levels(groups)[keep],
       definition = list(requested = options, breaks = breaks, right = TRUE,
                         labels = levels(groups)[keep], shares = shares[keep] / sum(shares)))
}
