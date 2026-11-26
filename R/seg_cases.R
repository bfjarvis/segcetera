# Splits seg_measure() input into measurement cases, keeping row-aligned weights
# with their observations. Labels retain their input type and first-seen order.
.seg_case_inputs <- function(checked, case, variables, population, graph = NULL) {
  if (length(case) != 1L) stop("Select one 'case' column", call. = FALSE)
  selected <- .seg_select(checked$data, case)
  ids <- selected[[1L]]
  if (!is.atomic(ids) || is.matrix(ids) || anyNA(ids) ||
      (is.numeric(ids) && any(!is.finite(ids))) ||
      any(!nzchar(as.character(ids)))) {
    stop("Case IDs must be nonmissing, nonempty finite values", call. = FALSE)
  }
  labels <- unique(ids)
  groups <- split(seq_along(ids), factor(match(ids, labels),
                                        levels = seq_along(labels)))
  inputs <- lapply(groups, function(rows) {
    part <- checked
    part$data <- checked$data[rows, , drop = FALSE]
    if (!is.null(part$coords)) part$coords <- part$coords[rows, , drop = FALSE]
    if (!is.null(part$geometry)) part$geometry <- part$geometry[rows]
    specs <- lapply(variables, function(spec) {
      if (!is.null(spec$weights)) spec$weights <- spec$weights[rows]
      spec
    })
    list(checked = part, variables = specs,
         population = if (!is.null(population)) population[rows],
         nodes = if (!is.null(graph)) graph$node[rows], rows = rows)
  })
  list(ids = labels, column = names(selected)[1L], inputs = inputs)
}

# Aligns support rows across cases for seg_measure(). Exact geometry/coordinate
# keys avoid numerical rounding; occurrence slots preserve coincident micro rows.
.seg_case_layout <- function(inputs) {
  sizes <- vapply(inputs, function(z) nrow(z$checked$data), integer(1))
  checks <- lapply(inputs, `[[`, "checked")
  reference <- checks[[1L]]
  same <- all(sizes == sizes[1L]) && all(vapply(inputs, function(z) {
    identical(z$nodes, inputs[[1L]]$nodes) &&
      identical(z$checked$geometry, reference$geometry) &&
      all(z$checked$coords == reference$coords)
  }, logical(1)))
  if (same) {
    return(list(positions = rep(list(seq_len(sizes[1L])), length(inputs)),
      n = sizes[1L], coords = reference$coords, geometry = reference$geometry,
      nodes = inputs[[1L]]$nodes, crs = reference$crs))
  }
  coords <- if (!is.null(checks[[1L]]$coords))
    do.call(rbind, lapply(checks, `[[`, "coords"))
  geometry <- if (!is.null(checks[[1L]]$geometry))
    do.call(c, lapply(checks, `[[`, "geometry"))
  nodes <- if (!is.null(inputs[[1L]]$nodes))
    unlist(lapply(inputs, `[[`, "nodes"), use.names = FALSE)
  if (!is.null(geometry)) {
    binary <- sf::st_as_binary(geometry)
    keys <- match(binary, unique(binary))
  } else if (!is.null(coords)) {
    # Integer axis IDs give exact, compact keys without formatting coordinates.
    keys <- paste(match(coords[, 1], unique(coords[, 1])),
                  match(coords[, 2], unique(coords[, 2])), sep = ":")
  } else keys <- nodes
  if (!is.null(nodes) && (!is.null(coords) || !is.null(geometry)))
    keys <- paste(nodes, keys, sep = ":")
  groups <- split(seq_along(keys), rep(seq_along(sizes), sizes))
  occurrence <- integer(length(keys))
  for (rows in groups) {
    sorted <- order(keys[rows])
    occurrence[rows[sorted]] <- sequence(rle(keys[rows[sorted]])$lengths)
  }
  slot_keys <- paste(keys, occurrence, sep = "/")
  slots <- match(slot_keys, unique(slot_keys))
  first <- match(seq_len(max(slots)), slots)
  positions <- lapply(groups, function(rows) slots[rows])
  list(positions = positions, n = length(first),
       coords = if (!is.null(coords)) coords[first, , drop = FALSE],
       geometry = if (!is.null(geometry)) geometry[first],
       nodes = if (!is.null(nodes)) nodes[first],
       crs = checks[[1L]]$crs)
}

# Pads already-compiled case inputs onto shared support. Compilation occurs on
# observed rows first, so absent rows add no categories, ranks, or population.
.seg_align_case <- function(prepared, positions, layout) {
  n <- layout$n
  if (!identical(positions, seq_len(n))) {
    prepared$compiled <- lapply(prepared$compiled, function(config) {
      mass <- numeric(n)
      mass[positions] <- config$mass
      code <- rep(-1L, n)
      code[positions] <- config$code
      if (nrow(config$counts)) {
        counts <- matrix(0, n, ncol(config$counts))
        counts[positions, ] <- config$counts
        colnames(counts) <- colnames(config$counts)
        config$counts <- counts
      }
      config$mass <- mass
      config$code <- code
      config
    })
    population <- numeric(n)
    population[positions] <- prepared$population
    prepared$population <- population
    original <- match(seq_len(n), positions)
    prepared$checked$data <- prepared$checked$data[original, , drop = FALSE]
    # Count columns represent zero population on absent support rows.
    count_columns <- if (prepared$surface_info$type == "raw") {
      unique(unlist(lapply(prepared$variables, function(spec) {
        if (!is.null(spec$counts)) names(.seg_select(prepared$checked$data, spec$counts))
      })))
    } else names(prepared$checked$data)
    for (column in count_columns) prepared$checked$data[is.na(original), column] <- 0
  }
  prepared$checked$coords <- layout$coords
  prepared$checked$geometry <- layout$geometry
  prepared$checked$crs <- layout$crs
  if (!is.null(prepared$surface_info$audit))
    prepared$surface_info$locations <- positions
  prepared
}

