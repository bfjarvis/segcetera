# Run from the repository root with an optimized installed package.
library(segcetera)
library(bench)

# Benchmark settings. Data generation and agreement checks precede timing.
N <- 10000L
side <- 1000
seed <- 2026L
bands <- c(50, 100, 200)
weighting <- "biweight"
power <- 3
normalize <- TRUE
neighbors <- "radius"
scope <- "overall"
iterations <- 10L
tolerance <- 1e-10
output_dir <- file.path("project", "benchmarks", "results")
pattern_args <- list(
  N = N, side = side, grid = "square", cycles = 2,
  group = c(A = 1, B = 1, C = 1),
  coef.score = c(group = 1, resource = 0)
)

set.seed(seed)
pattern <- do.call(seg_pattern, pattern_args)
variables <- list(group = seg_categorical("group"))

run_search <- function(search, output = "indices") {
  seg_measure(
    pattern, variables, coords = c("x", "y"), bands = bands,
    weighting = weighting, power = power, normalize = normalize,
    neighbors = neighbors, search = search, scope = scope, output = output
  )
}

# Compare environments as well as indices before timing indices-only runs.
reference <- run_search("brute", output = "full")
for (search in c("boost", "grid")) {
  candidate <- run_search(search, output = "full")
  agreement <- all.equal(reference$analyses, candidate$analyses,
                         tolerance = tolerance)
  if (!isTRUE(agreement)) {
    stop(search, " disagrees with brute: ", paste(agreement, collapse = "; "))
  }
}
grid <- candidate$diagnostics$grid
rm(reference, candidate)
message("Environments and indices agree across brute, boost, and grid.")

# Include setup and neighborhood construction on every call, and retain GC time.
# Search diagnostics differ, so compare the scientific results within each object.
timings <- bench::mark(
  brute = run_search("brute"),
  boost = run_search("boost"),
  grid = run_search("grid"),
  iterations = iterations,
  filter_gc = FALSE,
  check = function(x, y) {
    isTRUE(all.equal(x$analyses, y$analyses, tolerance = tolerance))
  }
)
print(timings[, c("expression", "min", "median", "itr/sec",
                 "mem_alloc", "n_itr", "n_gc")])

summary <- data.frame(
  search = as.character(timings$expression),
  median_seconds = as.numeric(timings$median),
  min_seconds = as.numeric(timings$min),
  r_allocated_bytes = as.numeric(timings$mem_alloc),
  iterations = timings$n_itr,
  garbage_collections = timings$n_gc
)
summary$speedup_vs_brute <- summary$median_seconds[summary$search == "brute"] /
  summary$median_seconds
print(summary, row.names = FALSE)

settings <- list(
  N = N, side = side, seed = seed, bands = bands,
  weighting = weighting, power = power, normalize = normalize,
  neighbors = neighbors, scope = scope, output = "indices",
  iterations = iterations, tolerance = tolerance,
  pattern = pattern_args,
  grid = grid
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
  system = Sys.info(),
  session = sessionInfo()
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
prefix <- file.path(output_dir, paste0("grid-search-",
                                     format(metadata$timestamp, "%Y%m%d-%H%M%S")))
saveRDS(list(settings = settings, metadata = metadata, summary = summary,
             timings = timings), paste0(prefix, ".rds"))
write.csv(summary, paste0(prefix, ".csv"), row.names = FALSE)
writeLines(capture.output(print(metadata)), paste0(prefix, "-session.txt"))
message("Saved benchmark results to ", prefix)
