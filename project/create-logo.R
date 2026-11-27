# Run with source("project/create-logo.R") in Positron/the R console,
# or Rscript project/create-logo.R. Saves the image used by the GitHub README.

library(ggplot2)
library(dplyr)
library(sf)
library(hexSticker)
library(systemfonts)

# Design parameters -----------------------------------------------------------
word <- "seg\n  cetera"
point_layout <- "halton" # "random", "halton", "sobol", or "hexagonal"
point_count <- 2^15L # Exact for random/quasi-random; approximate for hexagonal
point_size <- 0.25
seed <- 357
noise <- 0.01
letter_strength <- 100
center_strength <- 3 # Weight on background log-distance, outside letters
radius <- 0.94
colour_quantile <- 0.9
viridis_option <- "D"
word_bounds <- c(left = -0.65, right = 0.65, bottom = -0.5, top = 0.5)
font_family <- "Helvetica Neue"
font_style <- "Bold"
outline_tolerance <- 0.1 # Curve approximation, in font points
background <- "#101519"
border_colour <- "#303940"
border_size <- 0.65
output_dpi <- 300
output_file <- file.path("man", "figures", "logo.png")
show_logo <- TRUE # Display in interactive Positron/R sessions
beta_a <- 0.2
beta_b <- 0.7
line_height <- 0.7

# Find the output location -----------------------------------------------------
script <- grep("^--file=", commandArgs(), value = TRUE)
root <- if (length(script)) {
  file.path(dirname(sub("^--file=", "", script[1L])), "..")
} else if (file.exists("DESCRIPTION")) {
  "."
} else {
  ".."
}
output <- file.path(normalizePath(root), output_file)
dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)

# Use a temporary device for command-line sticker rendering.
# Interactive runs retain the session's normal plot device for print(logo).
if (!interactive()) {
  drawing_file <- tempfile(fileext = ".pdf")
  grDevices::pdf(drawing_file)
  drawing_device <- grDevices::dev.cur()
}

# Turn the lettering into polygon shapes --------------------------------------
fonts <- systemfonts::system_fonts()
font <- fonts[fonts$family == font_family & fonts$style == font_style, ]
if (!nrow(font)) {
  stop(
    "Font unavailable: ",
    font_family,
    " ",
    font_style,
    ". Choose a family/style from systemfonts::system_fonts()."
  )
}
font <- font[1L, ]

shaped <- systemfonts::shape_string(
  word,
  path = font$path,
  index = font$index,
  size = 100,
  lineheight = line_height
)$shape
outlines <- systemfonts::glyph_outline(
  shaped$index,
  path = font$path,
  index = font$index,
  size = 100,
  tolerance = outline_tolerance
)
# Layout controls and spaces advance the text without contributing shapes.
outlines <- outlines[!grepl("^[[:space:]]$", shaped$glyph[outlines$glyph]), ]
if (!nrow(outlines)) {
  stop("The lettering has no drawable glyph outlines.")
}
outlines$x <- outlines$x + shaped$x_offset[outlines$glyph]
outlines$y <- outlines$y + shaped$y_offset[outlines$glyph]

letter_geometry <- lapply(split(outlines, outlines$glyph), function(glyph) {
  rings <- lapply(split(glyph, glyph$contour), function(contour) {
    xy <- as.matrix(contour[c("x", "y")])
    if (any(xy[1L, ] != xy[nrow(xy), ])) {
      xy <- rbind(xy, xy[1L, ])
    }
    st_sfc(st_polygon(list(xy)))
  })
  # Even-odd filling retains counters (letter holes) and detached components.
  Reduce(st_sym_difference, rings)[[1L]]
})
letter_geometry <- st_sfc(letter_geometry)

# Fit the complete word into its bounds while preserving font proportions.
bounds <- st_bbox(letter_geometry)
text_center <- c(
  mean(bounds[c("xmin", "xmax")]),
  mean(bounds[c("ymin", "ymax")])
)
target_center <- c(
  mean(word_bounds[c("left", "right")]),
  mean(word_bounds[c("bottom", "top")])
)
text_scale <- min(
  (word_bounds["right"] - word_bounds["left"]) /
    (bounds["xmax"] - bounds["xmin"]),
  (word_bounds["top"] - word_bounds["bottom"]) /
    (bounds["ymax"] - bounds["ymin"])
)
letter_shapes <- st_sf(
  geometry = (letter_geometry - text_center) *
    as.numeric(text_scale) +
    target_center
)

# Choose point locations -------------------------------------------------------
set.seed(seed)
point_layout <- match.arg(
  point_layout,
  c("random", "halton", "sobol", "hexagonal")
)
if (point_layout == "random") {
  dots <- tibble(
    x = runif(
      2 * point_count,
      -sqrt(3) * radius / 2,
      sqrt(3) * radius / 2
    ),
    y = runif(2 * point_count, -radius, radius)
  ) |>
    filter(abs(y) <= radius - abs(x) / sqrt(3)) |>
    slice_head(n = point_count)
} else if (point_layout %in% c("halton", "sobol")) {
  if (!requireNamespace("randtoolbox", quietly = TRUE)) {
    stop("Halton and Sobol layouts require the 'randtoolbox' package.")
  }
  # Fixed low-discrepancy sequences cover the bounding box evenly.
  # Clip to the same hexagon as the other layouts, then retain the requested N.
  sequence_generator <- switch(
    point_layout,
    halton = randtoolbox::halton,
    sobol = randtoolbox::sobol
  )
  candidate_count <- 2 * point_count
  repeat {
    unit_points <- sequence_generator(candidate_count, dim = 2, init = TRUE)
    dots <- tibble(
      x = (unit_points[, 1] - 0.5) * sqrt(3) * radius,
      y = (unit_points[, 2] - 0.5) * 2 * radius
    ) |>
      filter(abs(y) <= radius - abs(x) / sqrt(3))
    if (nrow(dots) >= point_count) {
      break
    }
    candidate_count <- 2 * candidate_count
  }
  dots <- slice_head(dots, n = point_count)
} else {
  domain <- st_as_sfc(st_bbox(c(
    xmin = -sqrt(3) * radius / 2,
    ymin = -radius,
    xmax = sqrt(3) * radius / 2,
    ymax = radius
  )))
  grid_spacing <- radius * sqrt(3 / point_count)
  grid <- st_make_grid(
    domain,
    cellsize = grid_spacing,
    square = FALSE,
    what = "centers"
  )
  xy <- st_coordinates(grid)
  dots <- tibble(x = xy[, 1], y = xy[, 2]) |>
    filter(
      abs(x) <= sqrt(3) * radius / 2,
      abs(y) <= radius - abs(x) / sqrt(3)
    )
}


# Assign continuous resources by noisy spatial utility -------------------------
n_points <- nrow(dots)
point_shapes <- st_as_sf(dots, coords = c("x", "y"), remove = FALSE)
in_word <- lengths(st_intersects(point_shapes, letter_shapes)) > 0L

dots <- dots |>
  mutate(
    in_word = lengths(st_intersects(point_shapes, letter_shapes)) > 0L,
    log_distance = log(sqrt(x^2 + y^2 + 0.1) / radius),
    center_score = ifelse(in_word, 0, log_distance),
    letter_score = as.numeric(in_word),
    score = letter_strength *
      scale(letter_score) +
      center_strength * scale(center_score),
    score = as.numeric(scale(score)) + noise * rnorm(n_points)
  ) |>
  arrange(score) |>
  mutate(z = sort(rbeta(n_points, beta_a, beta_b)))

# Cap the display scale only; retain the original resource values.
colour_cap <- unname(quantile(dots$z, colour_quantile))
dots <- dots |> mutate(colour_value = pmin(z, colour_cap))
dot_plot <- ggplot(dots, aes(x, y, colour = z)) +
  geom_point(size = point_size, stroke = 0) +
  scale_colour_viridis_c() +
  coord_fixed(xlim = c(-1, 1), ylim = c(-1, 1), expand = FALSE) +
  theme_void() +
  theme(legend.position = "none", plot.margin = margin(0, 0, 0, 0))

# Save, then show the logo in an interactive plot pane --------------------------
# hexSticker places the zero-centered subplot in its own frame at (1, 1).
sticker_file <- tempfile(fileext = ".png")
logo <- hexSticker::sticker(
  dot_plot,
  package = "",
  p_family = font_family,
  u_family = font_family,
  s_x = 1,
  s_y = 1,
  s_width = 2,
  s_height = 2,
  h_fill = background,
  h_color = border_colour,
  h_size = border_size,
  filename = sticker_file,
  dpi = output_dpi
)
logo <- logo + theme(plot.margin = margin(2, 2, 2, 2, unit = "mm"))
hexSticker::save_sticker(output, logo, dpi = output_dpi)
unlink(sticker_file)
if (!interactive()) {
  invisible(grDevices::dev.off(drawing_device))
  unlink(drawing_file)
}
if (interactive() && show_logo) {
  print(logo)
}
message("Saved ", n_points, " ", point_layout, " points to: ", output)
