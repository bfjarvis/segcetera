# Run from the repository root with an optimized installed package.
library(segcetera)
library(bench)
library(ggplot2)

N <- 10000L
side <- 1000
seed <- 2026L
radius_fractions <- c(0.005, 0.01, 0.02, 0.03, 0.05, 0.075,
                      0.10, 0.125, 0.15, 0.20, 0.30, 0.40)
iterations <- 5L
sample_focals <- 256L
tolerance <- 1e-10
output_dir <- file.path("project", "benchmarks", "results")
pattern_args <- list(
  N = N, side = side, grid = "square", cycles = 2,
  group = c(A = 1, B = 1, C = 1),
  coef.score = c(group = 1, resource = 0)
)
weighting <- "biweight"
power <- 3
normalize <- TRUE
scope <- "overall"

# Optional command-line arguments: N, iterations, then radius/side fractions.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1L) N <- as.integer(args[1])
if (length(args) >= 2L) iterations <- as.integer(args[2])
if (length(args) >= 3L) radius_fractions <- as.numeric(args[-c(1, 2)])
if (any(!is.finite(radius_fractions) | radius_fractions <= 0)) {
  stop("Radius fractions must be positive finite numbers")
}
pattern_args$N <- N

set.seed(seed)
pattern <- do.call(seg_pattern, pattern_args)
variables <- list(group = seg_categorical("group"))
radii <- side * radius_fractions

# Estimate actual neighborhood occupancy, including boundaries and clustering.
# One distance vector per sampled focal row keeps this outside the timed calls.
focals <- sample.int(N, min(sample_focals, N))
counts <- vapply(focals, function(i) {
  d <- sqrt((pattern$x - pattern$x[i])^2 + (pattern$y - pattern$y[i])^2)
  findInterval(radii, sort(d))
}, integer(length(radii)))
mean_neighbors <- rowMeans(matrix(counts, nrow = length(radii)))

run_search <- function(search, radius, output = "indices") {
  seg_measure(
    pattern, variables, coords = c("x", "y"), bands = radius,
    neighbors = "radius", search = search, weighting = weighting,
    power = power, normalize = normalize, scope = scope, output = output
  )
}
expressions <- list(
  brute = quote(run_search("brute", radius)),
  boost = quote(run_search("boost", radius))
)
timings <- vector("list", length(radii))
rows <- vector("list", length(radii))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
run_id <- format(Sys.time(), "%Y%m%d-%H%M%S")
prefix <- file.path(output_dir, paste0("radius-crossover-", run_id))
for (i in seq_along(radii)) {
  radius <- radii[i]
  reference <- run_search("brute", radius, output = "full")
  candidate <- run_search("boost", radius, output = "full")
  agreement <- all.equal(reference$analyses, candidate$analyses,
                         tolerance = tolerance)
  if (!isTRUE(agreement)) {
    stop("Searches disagree at radius ", radius, ": ",
         paste(agreement, collapse = "; "))
  }
  rm(reference, candidate)
  # Alternate backend order between radii; include all GC time.
  timing <- bench::mark(
    exprs = if (i %% 2L) expressions else rev(expressions),
    iterations = iterations, filter_gc = FALSE,
    check = function(x, y) {
      isTRUE(all.equal(x$analyses, y$analyses, tolerance = tolerance))
    }
  )
  timings[[i]] <- timing
  elapsed <- setNames(as.numeric(timing$median), as.character(timing$expression))
  rows[[i]] <- data.frame(
    radius = radius,
    radius_pct_side = 100 * radius_fractions[i],
    mean_neighbors = mean_neighbors[i],
    neighbor_pct = 100 * mean_neighbors[i] / N,
    brute_ms = 1000 * elapsed[["brute"]],
    boost_ms = 1000 * elapsed[["boost"]],
    boost_speedup = elapsed[["brute"]] / elapsed[["boost"]]
  )
  print(rows[[i]], row.names = FALSE, digits = 4)
  write.csv(do.call(rbind, rows[seq_len(i)]), paste0(prefix, ".csv"),
            row.names = FALSE)
  flush.console()
}
summary <- do.call(rbind, rows)
rownames(summary) <- NULL
print(summary, row.names = FALSE, digits = 4)

comparison_plot <- ggplot(summary, aes(radius_pct_side, boost_speedup)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50") +
  geom_line(colour = "#008080") +
  geom_point(colour = "#008080", size = 2) +
  scale_x_log10(breaks = c(0.5, 1, 2, 5, 10, 20, 40)) +
  labs(x = "Radius (% of square side)",
       y = "Boost speedup relative to brute",
       title = paste(format(N, big.mark = ","), "agents; one radius per call")) +
  theme_minimal()
if (interactive()) print(comparison_plot)

settings <- list(
  seed = seed, pattern = pattern_args, radius_fractions = radius_fractions,
  iterations = iterations, sampled_focals = focals, tolerance = tolerance,
  weighting = weighting, power = power, normalize = normalize,
  scope = scope, output = "indices"
)
metadata <- list(
  timestamp = Sys.time(),
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
saveRDS(list(settings = settings, metadata = metadata, summary = summary,
             timings = timings), paste0(prefix, ".rds"))
write.csv(summary, paste0(prefix, ".csv"), row.names = FALSE)
writeLines(capture.output(print(metadata)), paste0(prefix, "-session.txt"))
ggsave(paste0(prefix, ".png"), comparison_plot, width = 7, height = 4, dpi = 160)
message("Saved crossover results to ", prefix)
