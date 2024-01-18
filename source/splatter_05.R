# derived from splatter_01

Rcpp::sourceCpp(here::here("source", "automaton_01.cpp"))

sample_set <- function(set, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  sample(set, size = 1)
}

sample_range <- function(min, max, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  sample(min:max, size = 1)
}

sample_uniform <- function(min, max, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  stats::runif(1, min, max)
}

sample_palette <- function(n, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  #if (runif(1) < .5) {
    palette_file <- here::here("source", "palette_02.csv")
    palette_data <- readr::read_csv(palette_file, show_col_types = FALSE)
    palette_base <- palette_data |>
      dplyr::slice_sample(n = 1) |>
      unlist()
  # } else {
  #   palette_file <- here::here("source", "palette_01.csv")
  #   palette_data <- readr::read_csv(palette_file, show_col_types = FALSE)
  #   palette_base <- palette_data |>
  #     dplyr::slice_sample(n = 1) |>
  #     dplyr::select(-source) |>
  #     unlist()
  # }
  palette_func <- colorRampPalette(palette_base)
  return(palette_func(n))
}

normalise <- function(values, min = 0, max = 1) {
  value_min <- min(values)
  value_max <- max(values)
  unit_values <- (values - value_min) / (value_max - value_min)
  scaled_values <- (unit_values + min) * (max - min)
  return(scaled_values)
}

automaton_values <- function(rows, cols, iterations, span) {
  mat <- automaton(rows, cols, iterations, span)
  val <- mat |>
    as.vector() |>
    normalise(min = .0001, max = 1) # scale to a "left-open unit interval"
  return(val)
}

create_base_image <- function(seed, rows, cols, shades, iterations, span) {
  set.seed(seed)
  if (is.null(cols)) cols <- rows
  img <- tidyr::expand_grid(
    row = seq(0, 50, length.out = rows),
    col = seq(0, 50, length.out = cols),
    size = 1
  )
  img$value <- automaton_values(rows, cols, iterations, span)
  return(img)
}

# helper function used only for debugging
show_base_image <- function(img) {
  ggplot2::ggplot(
    data = img,
    mapping = ggplot2::aes(x = col, y = row, fill = value)
  ) +
    ggplot2::geom_raster(show.legend = FALSE) +
    ggplot2::coord_equal() +
    ggplot2::theme_void()
}

is_within_circle <- function(x, y, size, origin_x = .5, origin_y = .5) {
  x <- normalise(x) - origin_x
  y <- normalise(y) - origin_y
  return(x^2 + y^2 < (size / 2)^2)
}

is_within_square <- function(x, y, size, origin_x = .5, origin_y = .5) {
  x <- normalise(x) - origin_x
  y <- normalise(y) - origin_y
  return(abs(x) < (size / 2) & abs(y) < (size / 2))
}

is_within_stripe <- function(x, y, size, origin_x = .5, origin_y = .5) {
  x <- normalise(x) - origin_x
  y <- normalise(y) - origin_y
  return(abs(x) < (size / 2))
}


filter_with_shape <- function(data, shape_type, shape_size, shape_contrast) {
  is_within_shape <- NULL
  if (shape_type == "circle") is_within_shape <- is_within_circle
  if (shape_type == "square") is_within_shape <- is_within_square
  if (shape_type == "stripe") is_within_shape <- is_within_stripe
  if (is.null(is_within_shape)) rlang::abort("unknown shape_type")
  
  data |>
    dplyr::mutate(
      shape = is_within_shape(
        x_coord,
        y_coord,
        shape_size,
        origin_x = sample_uniform(min = 0, max = 1),
        origin_y = sample_uniform(min = 0, max = 1)
      ),
      value = value + shape * shape_contrast * sample_set(c(-1, 1))
    )
}

create_layered_image <- function(
    data,
    shape_size,
    shape_contrast,
    shape_type,
    n_shapes
) {
  
  for(s in 1:n_shapes) {
    data <- filter_with_shape(data, shape_type, shape_size, shape_contrast)
  }
  data |> dplyr::mutate(value = normalise(value, min = .0001, max = 1))
}

curl_step <- function(data, iteration, scale, octaves, seed) {
  noise_data <- ambient::curl_noise(
    x = data$x_coord,
    y = data$y_coord,
    z = data$z_coord,
    seed = seed,
    generator = ambient::fracture,
    noise = ambient::gen_simplex,
    fractal = ambient::billow,
    octaves = octaves
  )
  data$iteration <- iteration
  data$x_coord <- data$x_coord + noise_data$x * scale
  data$y_coord <- data$y_coord + noise_data$y * scale
  data$z_coord <- data$z_coord + noise_data$z * scale
  return(data)
}

curl_loop <- function(data, seed, iterations, scale, octaves) {
  data$z_coord <- data$iteration <- 1
  state <- purrr::accumulate(
    .x = (1:iterations) + 1,
    .f = curl_step,
    .init = data,
    scale = scale,
    octaves = octaves,
    seed = seed
  )
  state <- dplyr::bind_rows(state)
  return(state)
}

create_curl_image <- function(data, seed, iterations, scale, octaves) {
  data |>
    dplyr::mutate(x_coord = col * .01, y_coord = row * .01) |>
    curl_loop(
      seed = seed,
      iterations = iterations,
      scale = scale,
      octaves = octaves
    ) |>
    dplyr::mutate(
      prop_complete = iteration / iterations,
      curl_strength = 1 - prop_complete
    )
}

compute_axis_limit <- function(data, column, trim = .04) {
  value <- data[[column]][data$iteration == 1]
  limit <- range(value) + c(1, -1) * trim
  return(limit)
}

compute_dot_size <- function(raw_size, max_size, curl_strength) {
  raw_size * max_size * curl_strength
}

compute_dot_shade <- function(value, shades) {
  shades[ceiling(value * length(shades))]
}

create_curl_plot <- function(data, shades, dot_scale) {
  data$dot_colour <- compute_dot_shade(data$value, shades)
  data$dot_size <- compute_dot_size(data$size, dot_scale, data$curl_strength)
  data |>
    ggplot2::ggplot(ggplot2::aes(
      x = x_coord,
      y = y_coord,
      colour = dot_colour,
      size = dot_size
    )) +
    ggplot2::geom_point(alpha = 1, stroke = 0, show.legend = FALSE) +
    ggplot2::coord_cartesian(
      xlim = compute_axis_limit(data, "x_coord"),
      ylim = compute_axis_limit(data, "y_coord")
    ) +
    ggplot2::scale_x_continuous(name = NULL, expand = c(0, 0), breaks = NULL) +
    ggplot2::scale_y_continuous(name = NULL, expand = c(0, 0), breaks = NULL) +
    ggplot2::scale_size_identity() +
    ggplot2::scale_colour_identity() +
    ggplot2::scale_alpha_continuous(range = c(0, 1)) +
    ggplot2::theme_void() +
    ggplot2::theme(plot.background = ggplot2::element_rect(fill = "white"))
}

create_directories <- function(dir, sizes, types) {
  if (!fs::dir_exists(dir)) fs::dir_create(dir)
  for (size in sizes) {
    for (type in types) {
      subdir <- fs::path(dir, type, size)
      if (!fs::dir_exists(subdir)) fs::dir_create(subdir)
    }
  }
}

file_path <- function(dir, type, size, file) {
  fs::path(dir, type, size, file)
}

file_name <- function(name, type) {
  paste0(name, ".", type)
}

create_image_file <- function(plot, size, type, dir, name) {
  file <- file_name(name, type)
  path <- file_path(dir, type, size, file)
  scaling <- 40 / 3
  ggplot2::ggsave(
    filename = path,
    plot = plot,
    width = scaling,
    height = scaling,
    dpi = size / scaling
  )
}

create_extra_images <- function(image_path, sizes, types, dir, name) {
  image_size <- max(sizes)
  image_type <- "png"
  for (size in sizes) {
    for (type in types) {
      is_original_image <- size == image_size & type == image_type
      if (!is_original_image) {
        file <- file_name(name, type)
        path <- file_path(dir, type, size, file)
        img <- magick::image_read(image_path)
        img <- magick::image_resize(img, paste0(size, "x", size))
        magick::image_write(img, path)
        rm(img)
        gc()
      }
    }
  }
}

create_image_parameters <- function(seed) {
  set.seed(seed)
  list(
    base_rows = 100,
    base_cols = NULL,
    base_iterations = 400000,
    base_span = sample_range(min = 1, max = 8),
    layer_shape_size = sample_uniform(min = 0, max = .3),
    layer_shape_contrast = sample_uniform(min = 0, max = .2),
    layer_shape_type = "circle",
    layer_n_shapes = sample_range(min = 10, max = 20),
    curl_scale = .0002,  #sample_uniform(min = .00003, max = .0001),
    curl_octaves = 5,
    curl_iterations = 100, #sample_range(min = 150, max = 400),
    plot_dot_scale = sample_range(min = 30, max = 40),
    plot_shades = sample_palette(n = 1024, seed = seed)
  )
}

art <- function(seed) {
  
  # system identification
  sys_id <- "05"
  sys_name <- "splatter"
  
  # administrative set up
  message_stem <- paste0("image seed ", seed, ":")
  output_dir <- here::here("output", sys_id)
  output_name <- paste0(sys_name, "_", sys_id, "_", seed)
  output_sizes <- c(500, 2000)
  output_types <- c("png", "jpg")
  create_directories(output_dir, output_sizes, output_types)
  
  message(paste(message_stem, "making image parameters"))
  params <- create_image_parameters(seed)
  
  message(paste(message_stem, "making data for base image"))
  image_data <- create_base_image(
    seed = seed,
    rows = params$base_rows,
    cols = params$base_cols,
    iterations = params$base_iterations,
    span = params$base_span
  )
  
  message(paste(message_stem, "making data for curled image"))
  image_data <- create_curl_image(
    data = image_data,
    seed = seed,
    iterations = params$curl_iterations,
    scale = params$curl_scale,
    octaves = params$curl_octaves
  )
  
  # message(paste(message_stem, "making foreground shape"))
  # image_data <- create_layered_image(
  #   data = image_data,
  #   shape_size = params$layer_shape_size,
  #   shape_contrast = params$layer_shape_contrast,
  #   shape_type = params$layer_shape_type,
  #   n_shapes = params$layer_n_shapes
  # )
  
  message(paste(message_stem, "making plot specification"))
  curl_plot <- create_curl_plot(
    data = image_data,
    shades = params$plot_shades,
    dot_scale = params$plot_dot_scale
  )
  
  message(paste(message_stem, "making primary image file"))
  image_path <- create_image_file(
    plot = curl_plot,
    size = max(output_sizes),
    type = "png",
    dir = output_dir,
    name = output_name
  )
  
  message(paste(message_stem, "making additional image files"))
  create_extra_images(
    image_path = image_path,
    sizes = output_sizes,
    types = output_types,
    dir = output_dir,
    name = output_name
  )
}

if (TRUE) for (s in 500:599) art(s)
