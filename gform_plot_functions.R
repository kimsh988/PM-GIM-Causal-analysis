############################################################
# Plot 
############################################################

.check_ggplot2 <- function() {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("The ggplot2 package is required for plotting.")
  }
}

# Extract variable-specific labels from a named vector/list
.get_label <- function(labels, variable, default) {
  if (is.null(labels)) {
    return(default)
  }
  
  if (!is.null(names(labels)) && variable %in% names(labels)) {
    return(labels[[variable]])
  }
  
  default
}

# Check whether a variable in dt_cov_plot has a continuous plotting structure
.is_continuous_plot <- function(x) {
  if (is.null(x)) return(FALSE)
  all(c("t_yov", "cov", "legend") %in% names(x))
}

# Check whether a variable in dt_cov_plot has a categorical plotting structure
.is_categorical_plot <- function(x, variable_name) {
  if (is.null(x)) return(FALSE)
  all(c("t0", "V1", "legend", variable_name) %in% names(x))
}

# Convert a continuous variable to a data frame for plotting
.as_continuous_plot_df <- function(x, variable_name) {
  data.frame(
    variable = variable_name,
    time = x$t_yov,
    value = x$cov,
    source = x$legend,
    stringsAsFactors = FALSE
  )
}

# Convert a categorical variable to a data frame for plotting
.as_categorical_plot_df <- function(x, variable_name) {
  data.frame(
    variable = variable_name,
    time = x$t0,
    level = x[[variable_name]],
    value = x$V1,
    source = x$legend,
    stringsAsFactors = FALSE
  )
}

# Continuous variables: plot selected variables as line graphs in a single figure
plot_gform_continuous_covariates <- function(
    gform,
    variables,
    dt_cov_plot = gform$dt_cov_plot,
    x_label = "Time",
    y_label = "Mean",
    variable_labels = NULL,
    source_labels = NULL,
    free_y = TRUE,
    point_size = 2,
    line_width = 0.7) {
  
  .check_ggplot2()
  
  plot_data <- do.call(
    rbind,
    lapply(
      variables,
      function(v) .as_continuous_plot_df(dt_cov_plot[[v]], v)
    )
  )
  
  plot_data$variable <- factor(plot_data$variable, levels = variables)
  
  if (!is.null(variable_labels)) {
    plot_data$variable <- dplyr::recode(
      plot_data$variable,
      !!!variable_labels
    )
  }
  
  if (!is.null(source_labels)) {
    plot_data$source <- dplyr::recode(
      plot_data$source,
      !!!source_labels
    )
  }
  
  scale_y <- if (free_y) "free_y" else "fixed"
  
  ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = time,
      y = value,
      color = source,
      linetype = source,
      shape = source,
      group = source
    )
  ) +
    ggplot2::geom_line(linewidth = line_width) +
    ggplot2::geom_point(size = point_size) +
    ggplot2::facet_wrap(~ variable, scales = scale_y) +
    ggplot2::labs(
      x = x_label,
      y = y_label,
      color = "Estimate",
      linetype = "Estimate",
      shape = "Estimate"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      strip.background = ggplot2::element_rect(fill = "white"),
      strip.text = ggplot2::element_text(face = "bold"),
      legend.position = "bottom"
    )
}

# Categorical variables: plot selected variables as bar graphs in a single figure
plot_gform_categorical_covariates <- function(
    gform,
    variables,
    dt_cov_plot = gform$dt_cov_plot,
    x_label = "Category",
    y_label = "Proportion",
    variable_labels = NULL,
    source_labels = NULL,
    level_labels = NULL,
    position = c("dodge", "stack"),
    free_y = TRUE,
    bar_width = 0.7) {
  
  .check_ggplot2()
  
  if (!requireNamespace("ggpattern", quietly = TRUE)) {
    stop("The ggpattern package is required for patterned bar plots.")
  }
  
  position <- match.arg(position)
  
  plot_data <- do.call(
    rbind,
    lapply(
      variables,
      function(v) .as_categorical_plot_df(dt_cov_plot[[v]], v)
    )
  )
  
  plot_data$variable <- factor(plot_data$variable, levels = variables)
  plot_data$level <- factor(plot_data$level)
  
  if (!is.null(variable_labels)) {
    plot_data$variable <- dplyr::recode(
      plot_data$variable,
      !!!variable_labels
    )
  }
  
  if (!is.null(source_labels)) {
    plot_data$source <- dplyr::recode(
      plot_data$source,
      !!!source_labels
    )
  }
  
  if (!is.null(level_labels)) {
    plot_data$level <- dplyr::recode(
      plot_data$level,
      !!!level_labels
    )
  }
  
  bar_position <- if (position == "dodge") {
    ggplot2::position_dodge(width = 0.8)
  } else {
    "stack"
  }
  
  scale_y <- if (free_y) "free_y" else "fixed"
  
  ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = level,
      y = value,
      fill = source,
      pattern = source,
      group = source
    )
  ) +
    ggpattern::geom_col_pattern(
      color = "black",
      position = bar_position,
      width = bar_width,
      pattern_fill = "black",
      pattern_color = "black",
      pattern_angle = 45,
      pattern_density = 0.08,
      pattern_spacing = 0.03,
      pattern_key_scale_factor = 0.6
    ) +
    ggpattern::scale_pattern_manual(
      values = c(
        "Parametric g-formula estimate" = "stripe",
        "Nonparametric estimate" = "none",
        "parametric g-formula estimate" = "stripe",
        "parametric g-formula estimates" = "stripe",
        "nonparametric estimate" = "none",
        "nonparametric estimates" = "none"
      )
    ) +
    ggplot2::facet_grid(
      variable ~ time,
      scales = scale_y
    ) +
    ggplot2::labs(
      x = x_label,
      y = y_label,
      fill = "Estimate",
      alpha = "Estimate"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      strip.background = ggplot2::element_rect(fill = "white"),
      strip.text = ggplot2::element_text(face = "bold"),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      legend.position = "bottom"
    )
}

# Wrapper function: automatically use line plots for continuous variables and bar plots for categorical variables
plot_gform_covariates <- function(
    gform = NULL,
    variables = NULL,
    dt_cov_plot = gform$dt_cov_plot,
    x_label_continuous = "Time",
    y_label_continuous = "Mean",
    x_label_categorical = "Category",
    y_label_categorical = "Proportion",
    variable_labels = NULL,
    source_labels = NULL,
    level_labels = NULL,
    categorical_position = c("dodge", "stack"),
    free_y_continuous = TRUE,
    free_y_categorical = TRUE,
    point_size = 2,
    line_width = 0.7,
    bar_width = 0.7) {
  
  categorical_position <- match.arg(categorical_position)
  
  if (is.null(dt_cov_plot)) {
    stop("dt_cov_plot is NULL. Please check gform$dt_cov_plot.")
  }
  
  available_variables <- names(dt_cov_plot)
  
  if (is.null(variables)) {
    variables <- available_variables
  }
  
  missing_variables <- setdiff(variables, available_variables)
  
  if (length(missing_variables) > 0) {
    warning(
      "These variable(s) were not found in dt_cov_plot: ",
      paste(missing_variables, collapse = ", "),
      call. = FALSE
    )
  }
  
  selected_variables <- intersect(variables, available_variables)
  
  if (length(selected_variables) == 0) {
    stop("No eligible variable was selected for plotting.")
  }
  
  continuous_vars <- selected_variables[
    sapply(
      selected_variables,
      function(v) .is_continuous_plot(dt_cov_plot[[v]])
    )
  ]
  
  categorical_vars <- selected_variables[
    sapply(
      selected_variables,
      function(v) .is_categorical_plot(dt_cov_plot[[v]], v)
    )
  ]
  
  unidentified_vars <- setdiff(
    selected_variables,
    union(continuous_vars, categorical_vars)
  )
  
  if (length(unidentified_vars) > 0) {
    warning(
      "These variable(s) were not recognized as either continuous or categorical plot data: ",
      paste(unidentified_vars, collapse = ", "),
      call. = FALSE
    )
  }
  
  plots <- list()
  
  if (length(continuous_vars) > 0) {
    plots$continuous <- plot_gform_continuous_covariates(
      gform = gform,
      variables = continuous_vars,
      dt_cov_plot = dt_cov_plot,
      x_label = x_label_continuous,
      y_label = y_label_continuous,
      variable_labels = variable_labels,
      source_labels = source_labels,
      free_y = free_y_continuous,
      point_size = point_size,
      line_width = line_width
    )
  }
  
  if (length(categorical_vars) > 0) {
    plots$categorical <- plot_gform_categorical_covariates(
      gform = gform,
      variables = categorical_vars,
      dt_cov_plot = dt_cov_plot,
      x_label = x_label_categorical,
      y_label = y_label_categorical,
      variable_labels = variable_labels,
      source_labels = source_labels,
      level_labels = level_labels,
      position = categorical_position,
      free_y = free_y_categorical,
      bar_width = bar_width
    )
  }
  
  if (length(plots) == 0) {
    stop("No eligible variable was selected for plotting.")
  }
  
  if (length(plots) == 1) {
    return(plots[[1]])
  }
  
  plots
}

# Short alias
plot_selected_variables <- plot_gform_covariates

plot_gform_outcome <- function(
    gform = NULL,
    dt_out_plot = NULL,
    outcome_type = c("risk", "survival"),
    x_label = "Years since baseline",
    y_label = NULL,
    source_labels = NULL,
    y_limits = NULL,
    point_size = 2,
    line_width = 0.7) {
  
  outcome_type <- match.arg(outcome_type)
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("The ggplot2 package is required for plotting.")
  }
  
  if (is.null(dt_out_plot)) {
    if (!is.null(gform) && "dt_out_plot" %in% names(gform)) {
      dt_out_plot <- gform$dt_out_plot
    }
  }
  
  if (is.null(dt_out_plot)) {
    stop("dt_out_plot is NULL. Pass `gform = gform` or `dt_out_plot = gform$dt_out_plot`.")
  }
  
  required_names <- c("t_yov", outcome_type, "legend")
  
  if (!all(required_names %in% names(dt_out_plot))) {
    stop(
      "dt_out_plot must contain: ",
      paste(required_names, collapse = ", "),
      ". Available elements are: ",
      paste(names(dt_out_plot), collapse = ", ")
    )
  }
  
  if (is.null(y_label)) {
    y_label <- if (outcome_type == "risk") "Cumulative risk" else "Survival probability"
  }
  
  plot_data <- data.frame(
    time = dt_out_plot$t_yov,
    value = dt_out_plot[[outcome_type]],
    source = dt_out_plot$legend,
    stringsAsFactors = FALSE
  )
  
  if (!is.null(source_labels)) {
    plot_data$source <- dplyr::recode(
      plot_data$source,
      !!!source_labels
    )
  }
  
  p <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = time,
      y = value,
      color = source,
      linetype = source,
      shape = source,
      group = source
    )
  ) +
    ggplot2::geom_line(linewidth = line_width) +
    ggplot2::geom_point(size = point_size) +
    ggplot2::labs(
      x = x_label,
      y = y_label,
      color = "Estimate",
      linetype = "Estimate",
      shape = "Estimate"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      legend.position = "bottom"
    )
  
  if (!is.null(y_limits)) {
    p <- p + ggplot2::coord_cartesian(ylim = y_limits)
  }
  
  p
}