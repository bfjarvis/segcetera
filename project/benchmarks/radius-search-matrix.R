# Run from the repository root with an optimized installed package.
library(segcetera)
library(bench)
library(ggplot2)

agent_counts <- c(1000L, 2500L, 5000L, 10000L, 20000L, 40000L, 80000L)
band_basis <- "distance_quantile" # Or "side_fraction" for the earlier scan.
probs <- eval(formals(seg_default_bands)$probs)
radius_fractions <- c(0.01, 0.025, 0.05, 0.10, 0.15, 0.20, 0.30, 0.40)
side <- 1000
seed <- 2026L
iterations <- 3L
sample_focals <- 256L
tolerance <- 1e-10
resume <- TRUE
output_dir <- file.path("project", "benchmarks", "results")
run_name <- "radius-default-band-matrix"
pattern_args <- list(
  side = side, grid = "square", cycles = 2,
  group = c(A = 1, B = 1, C = 1),
  coef.score = c(group = 1, resource = 0)
)
weighting <- "biweight"
power <- 3
normalize <- TRUE
scope <- "overall"
band_basis <- match.arg(band_basis, c("distance_quantile", "side_fraction"))
band_values <- if (band_basis == "distance_quantile") probs else radius_fractions
band_label <- if (band_basis == "distance_quantile")
  "Sampled pairwise-distance percentile (%)" else "Radius (% of square side)"

settings <- list(
  agent_counts = agent_counts, band_basis = band_basis, band_values = band_values,
  seed = seed, iterations = iterations, sample_focals = sample_focals,
  tolerance = tolerance, pattern = pattern_args, weighting = weighting,
  power = power, normalize = normalize, scope = scope, output = "indices"
)
build <- tools::md5sum(c(
  file.path(find.package("segcetera"), "R", "segcetera.rdb"),
  getLoadedDLLs()[["segcetera"]][["path"]]
))
metadata <- list(
  started = Sys.time(),
  commit = system2("git", c("rev-parse", "HEAD"), stdout = TRUE),
  working_tree = system2("git", c("status", "--porcelain"), stdout = TRUE),
  package_version = as.character(packageVersion("segcetera")),
  package_path = find.package("segcetera"),
  compiler = system2(file.path(R.home("bin"), "R"),
                     c("CMD", "config", "CXX17"), stdout = TRUE),
  compiler_flags = system2(file.path(R.home("bin"), "R"),
                           c("CMD", "config", "CXX17FLAGS"), stdout = TRUE),
  system = Sys.info(), session = sessionInfo()
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
prefix <- file.path(output_dir, run_name)
checkpoint <- paste0(prefix, ".rds")
rows <- list()
timings <- list()
if (resume && file.exists(checkpoint)) {
  previous <- readRDS(checkpoint)
  if (!identical(previous$settings, settings) || !identical(previous$build, build)) {
    stop("Checkpoint settings/build differ; choose a new 'run_name' or set resume = FALSE")
  }
  rows <- previous$rows
  timings <- previous$timings
  metadata <- previous$metadata
  message("Resuming ", length(rows), " completed comparisons.")
}

run_search <- function(search, radius) {
  seg_measure(
    pattern, variables, coords = c("x", "y"), bands = radius,
    neighbors = "radius", search = search, weighting = weighting,
    power = power, normalize = normalize, scope = scope, output = "indices"
  )
}
expressions <- list(
  brute = quote(run_search("brute", radius)),
  boost = quote(run_search("boost", radius))
)

for (N in agent_counts) {
  set.seed(seed)
  pattern <- do.call(seg_pattern, c(list(N = N), pattern_args))
  variables <- list(group = seg_categorical("group"))
  radii <- if (band_basis == "distance_quantile") {
    seg_default_bands(pattern, coords = c("x", "y"), probs = probs)
  } else side * radius_fractions
  message("N = ", N, ": ", length(radii), " single-band comparisons; radii = ",
          paste(format(round(radii, 2), trim = TRUE), collapse = ", "))
  focals <- sample.int(N, min(sample_focals, N))
  counts <- vapply(focals, function(i) {
    d <- sqrt((pattern$x - pattern$x[i])^2 + (pattern$y - pattern$y[i])^2)
    findInterval(radii, sort(d))
  }, integer(length(radii)))
  mean_neighbors <- rowMeans(matrix(counts, nrow = length(radii)))

  for (i in seq_along(radii)) {
    key <- paste(N, i, sep = "-")
    if (key %in% names(rows)) next
    radius <- radii[i]
    # Rotate backend order between cells. Check indices; measure elapsed time
    # including GC and setup. Native allocations are outside bench's memory view.
    timing <- bench::mark(
      exprs = if ((match(N, agent_counts) + i) %% 2L) expressions else rev(expressions),
      iterations = iterations, filter_gc = FALSE, memory = FALSE,
      check = function(x, y) {
        isTRUE(all.equal(x$analyses, y$analyses, tolerance = tolerance))
      }
    )
    elapsed <- setNames(as.numeric(timing$median), as.character(timing$expression))
    rows[[key]] <- data.frame(
      N = N, band_value = band_values[i], radius = radius,
      radius_pct_side = 100 * radius / side,
      mean_neighbors = mean_neighbors[i], neighbor_pct = 100 * mean_neighbors[i] / N,
      brute_ms = 1000 * elapsed[["brute"]],
      boost_ms = 1000 * elapsed[["boost"]],
      boost_time_saved_pct = 100 * (1 - elapsed[["boost"]] / elapsed[["brute"]]),
      boost_time_saved_ms = 1000 * (elapsed[["brute"]] - elapsed[["boost"]]),
      boost_speedup = elapsed[["brute"]] / elapsed[["boost"]]
    )
    timings[[key]] <- timing
    summary <- do.call(rbind, rows)
    rownames(summary) <- NULL
    saveRDS(list(settings = settings, build = build, metadata = metadata,
                 rows = rows, timings = timings), checkpoint)
    write.csv(summary, paste0(prefix, ".csv"), row.names = FALSE)
    print(rows[[key]], row.names = FALSE, digits = 4)
    flush.console()
  }
}

summary <- do.call(rbind, rows)
rownames(summary) <- NULL
axes <- ggplot(summary, aes(
  x = factor(100 * band_value, levels = 100 * band_values),
  y = factor(N, levels = agent_counts, labels = format(agent_counts, big.mark = ","))
)) +
  labs(x = band_label, y = "Agents",
       subtitle = paste(iterations, "timed repeats per backend; positive values favor boost")) +
  theme_minimal() +
  theme(panel.grid = element_blank())
relative_plot <- axes +
  geom_tile(aes(fill = boost_time_saved_pct), colour = "white") +
  geom_text(aes(label = sprintf("%+.0f%%", boost_time_saved_pct)), size = 3.5) +
  scale_fill_gradient2(low = "#cc6666", mid = "white", high = "#65b5b0",
                       midpoint = 0, name = "Boost time\nsaved (%)") +
  labs(title = "Radius search: boost versus brute")
absolute_plot <- axes +
  geom_tile(aes(fill = boost_time_saved_ms), colour = "white") +
  geom_text(aes(label = ifelse(abs(boost_time_saved_ms) < 1,
    sprintf("%+.2f", boost_time_saved_ms),
    sprintf("%+.0f", boost_time_saved_ms))), size = 3.5) +
  scale_fill_gradient2(low = "#cc6666", mid = "white", high = "#65b5b0",
                       midpoint = 0, name = "Boost time\nsaved (ms)") +
  labs(title = "Radius search: absolute timing difference")
if (interactive()) print(relative_plot)
ggsave(paste0(prefix, "-relative.png"), relative_plot, width = 10, height = 5.5, dpi = 160)
ggsave(paste0(prefix, "-absolute.png"), absolute_plot, width = 10, height = 5.5, dpi = 160)
writeLines(capture.output(print(metadata)), paste0(prefix, "-session.txt"))
message("Saved matrix results and heatmaps to ", prefix)
