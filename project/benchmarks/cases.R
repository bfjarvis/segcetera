# Run from the repository root with an optimized installed package.
library(segcetera)
library(bench)

N <- 2500L
case_count <- 4L
iterations <- 3L
seed <- 2026L
side <- 1000
radius_bands <- c(50, 100, 200)
knn_bands <- c(10, 100, 300)
output_dir <- file.path("project", "benchmarks", "results")
pattern_args <- list(N = N, side = side, grid = "square",
  group = c(A = 1, B = 1, C = 1), coef.score = c(group = 1, resource = 0))

set.seed(seed)
pattern <- do.call(seg_pattern, pattern_args)
observations <- do.call(rbind, lapply(seq_len(case_count), function(i) {
  d <- pattern
  d$case_id <- i
  shuffled <- sample.int(N, floor(N * (i - 1) / case_count))
  d$group[shuffled] <- sample(levels(d$group), length(shuffled), replace = TRUE)
  d
}))
variables <- list(group = seg_categorical("group"))
separate_inputs <- split(observations, observations$case_id)
rows <- list()
timings <- list()
for (mode in c("radius", "knn")) for (search in c("brute", "boost", "grid")) {
  bands <- if (mode == "radius") radius_bands else knn_bands
  measure <- function(x, ...) seg_measure(x, variables, coords = c("x", "y"),
    neighbors = mode, search = search, bands = bands, scope = "overall", ...)
  run_separate <- function() lapply(separate_inputs, function(d) measure(d)$analyses)
  run_shared <- function() lapply(measure(observations, case = "case_id")$cases,
                                  `[[`, "analyses")
  proof <- measure(observations, case = "case_id")
  stopifnot(proof$diagnostics$focal_searches == N)
  key <- paste(mode, search, sep = "-")
  timing <- bench::mark(separate = run_separate(), shared = run_shared(),
    iterations = iterations, filter_gc = FALSE, memory = FALSE, check = TRUE)
  elapsed <- setNames(as.numeric(timing$median), as.character(timing$expression))
  rows[[key]] <- data.frame(neighbors = mode, search = search, N = N,
    cases = case_count, separate_ms = 1000 * elapsed[["separate"]],
    shared_ms = 1000 * elapsed[["shared"]],
    speedup = elapsed[["separate"]] / elapsed[["shared"]],
    separate_searches = N * case_count,
    shared_searches = proof$diagnostics$focal_searches)
  timings[[key]] <- timing
  print(rows[[key]], row.names = FALSE)
}
summary <- do.call(rbind, rows)
rownames(summary) <- NULL
metadata <- list(timestamp = Sys.time(), session = sessionInfo(), system = Sys.info(),
  package_path = find.package("segcetera"),
  compiler = system2(file.path(R.home("bin"), "R"),
                     c("CMD", "config", "CXX17"), stdout = TRUE),
  compiler_flags = system2(file.path(R.home("bin"), "R"),
                           c("CMD", "config", "CXX17FLAGS"), stdout = TRUE),
  commit = system2("git", c("rev-parse", "HEAD"), stdout = TRUE),
  working_tree = system2("git", c("status", "--porcelain"), stdout = TRUE))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
prefix <- file.path(output_dir, paste0("cases-", format(metadata$timestamp, "%Y%m%d-%H%M%S")))
saveRDS(list(settings = list(N = N, cases = case_count, iterations = iterations,
  seed = seed, pattern = pattern_args, radius_bands = radius_bands,
  knn_bands = knn_bands), summary = summary, timings = timings,
  metadata = metadata), paste0(prefix, ".rds"))
write.csv(summary, paste0(prefix, ".csv"), row.names = FALSE)
writeLines(capture.output(print(metadata)), paste0(prefix, "-session.txt"))
message("Saved case benchmark to ", prefix)
