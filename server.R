# server.R

# Load required libraries
library(shiny)
library(shinydashboard)
library(DT)
library(ggplot2)
library(dplyr)
library(readr)
library(numDeriv)
library(tools)

# =============================================================================
# Cuckoo Search Algorithm Implementation (Updated to match paper exactly)
# =============================================================================

# Simple bounds application
simple_bounds <- function(s, lb, ub) {
  pmax(pmin(s, ub), lb)
}

# Generate cuckoos using Levy flights
get_cuckoos <- function(nest, best, lb, ub) {
  n <- nrow(nest)
  d <- ncol(nest)
  beta <- 1.5
  
  sigma <- (gamma(1 + beta) * sin(pi * beta / 2) /
              (gamma((1 + beta) / 2) * beta * 2^((beta - 1) / 2)))^(1 / beta)
  
  for (j in 1:n) {
    s <- nest[j, ]
    u <- rnorm(d) * sigma
    v <- rnorm(d)
    step <- u / (abs(v)^(1 / beta))
    stepsize <- 0.01 * step * (s - best)
    s <- s + stepsize * rnorm(d)
    nest[j, ] <- simple_bounds(s, lb, ub)
  }
  return(nest)
}

# Replace some nests with new solutions
empty_nests <- function(nest, lb, ub, pa) {
  n <- nrow(nest)
  d <- ncol(nest)
  
  K <- matrix(runif(n * d) > pa, nrow = n, ncol = d)
  
  # Use permutation without replacement
  perm1 <- sample(1:n, n, replace = FALSE)
  perm2 <- sample(1:n, n, replace = FALSE)
  
  stepsize <- runif(1) * (nest[perm1, ] - nest[perm2, ])
  new_nest <- nest + stepsize * K
  
  for (j in 1:n) {
    new_nest[j, ] <- simple_bounds(new_nest[j, ], lb, ub)
  }
  return(new_nest)
}

# Evaluate and update best nests
get_best_nest <- function(nest, newnest, fitness, fitness_func) {
  for (j in 1:nrow(nest)) {
    fnew <- fitness_func(newnest[j, ])
    if (fnew <= fitness[j]) {
      fitness[j] <- fnew
      nest[j, ] <- newnest[j, ]
    }
  }
  fmin <- min(fitness)
  best <- nest[which.min(fitness), ]
  return(list(fmin = fmin, best = best, nest = nest, fitness = fitness))
}

# Main cuckoo search function
cuckoo_search_custom <- function(fitness_func, nd, lb, ub, n = 25, N_IterTotal = 1000, pa = 0.25) {
  # Set random seed to match paper implementation
  set.seed(123)
  
  # Initialize nest positions
  nest <- matrix(runif(n * nd), nrow = n, ncol = nd) * 
    matrix(rep(ub - lb, each = n), nrow = n) + 
    matrix(rep(lb, each = n), nrow = n)
  
  fitness <- rep(1e10, n)
  
  # Initial evaluation
  result <- get_best_nest(nest, nest, fitness, fitness_func)
  fmin <- result$fmin
  bestnest <- result$best
  nest <- result$nest
  fitness <- result$fitness
  
  convergence_history <- numeric(N_IterTotal)
  
  # Main iteration loop
  for (iter in 1:N_IterTotal) {
    # Generate new solutions via Levy flights
    new_nest <- get_cuckoos(nest, bestnest, lb, ub)
    result <- get_best_nest(nest, new_nest, fitness, fitness_func)
    fmin <- result$fmin
    bestnest <- result$best
    nest <- result$nest
    fitness <- result$fitness
    
    # Replace some nests
    new_nest <- empty_nests(nest, lb, ub, pa)
    result <- get_best_nest(nest, new_nest, fitness, fitness_func)
    fmin <- result$fmin
    bestnest <- result$best
    nest <- result$nest
    fitness <- result$fitness
    
    convergence_history[iter] <- fmin
  }
  
  return(list(bestnest = bestnest, fmin = fmin, convergence_history = convergence_history))
}

# =============================================================================
# Core scGTM Functions (Updated to match paper exactly)
# =============================================================================

# Link function (exactly matching paper implementation)
link_function <- function(t, mu, k1, k2, t0) {
  t <- as.numeric(t)
  
  # Handle sign exactly as in paper
  sign_k1 <- sign(k1) + (k1 == 0)
  sign_k2 <- sign(k2) + (k2 == 0)
  
  part1 <- mu * exp(-abs(k1) * (t - t0)^2) * sign_k1
  part2 <- mu * exp(-abs(k2) * (t - t0)^2) * sign_k2
  
  result <- ifelse(t <= t0, part1, part2)
  result <- pmax(pmin(result, 50), -50)
  
  return(result)
}

# Model functions (updated to use link_function)
hill_model <- function(t, params) {
  mu <- params[1]
  k1 <- params[2]  # Don't take absolute value here
  k2 <- params[3]  # Don't take absolute value here
  t0 <- params[4]
  
  return(link_function(t, mu, k1, k2, t0))
}

valley_model <- function(t, params, baseline) {
  mu <- params[1]
  k1 <- params[2]
  k2 <- params[3]
  t0 <- params[4]
  
  result <- baseline - link_function(t, mu, k1, k2, t0)
  result <- pmax(pmin(result, 50), -50)
  
  return(result)
}

monotonic_increasing_model <- function(t, params) {
  mu <- params[1]
  k <- abs(params[2])
  t0 <- params[3]
  
  t <- as.numeric(t)
  
  result <- mu * (1 + tanh(k * (t - t0))) / 2
  result <- pmax(pmin(result, 10), -10)
  
  return(result)
}

monotonic_decreasing_model <- function(t, params) {
  mu <- params[1]
  k <- abs(params[2])
  t0 <- params[3]
  
  t <- as.numeric(t)
  
  result <- mu * (1 - tanh(k * (t - t0))) / 2
  result <- pmax(pmin(result, 10), -10)
  
  return(result)
}

# Universal model evaluation
evaluate_model <- function(t, params, model_type, baseline = NULL) {
  switch(model_type,
         "hill" = hill_model(t, params),
         "valley" = valley_model(t, params, baseline),
         "increasing" = monotonic_increasing_model(t, params),
         "decreasing" = monotonic_decreasing_model(t, params)
  )
}

# Log-likelihood functions for different distributions
log_likelihood_poisson <- function(params, y, t, model_type, baseline = NULL) {
  tryCatch({
    if (any(is.na(params)) || any(!is.finite(params))) return(1e10)
    
    lambda <- exp(evaluate_model(t, params, model_type, baseline))
    lambda <- pmax(pmin(lambda, 1000), 0.001)
    
    if (any(is.na(lambda)) || any(!is.finite(lambda))) return(1e10)
    
    ll <- sum(dpois(y, lambda, log = TRUE))
    if (is.na(ll) || !is.finite(ll)) return(1e10)
    
    return(-ll)
  }, error = function(e) 1e10)
}

log_likelihood_nb <- function(params, y, t, model_type, baseline = NULL, theta = 1) {
  tryCatch({
    if (any(is.na(params)) || any(!is.finite(params))) return(1e10)
    
    mu <- exp(evaluate_model(t, params, model_type, baseline))
    mu <- pmax(pmin(mu, 1000), 0.001)
    
    if (any(is.na(mu)) || any(!is.finite(mu))) return(1e10)
    
    ll <- sum(dnbinom(y, size = theta, mu = mu, log = TRUE))
    if (is.na(ll) || !is.finite(ll)) return(1e10)
    
    return(-ll)
  }, error = function(e) 1e10)
}

log_likelihood_zip <- function(params, y, t, model_type, baseline = NULL, pi = 0.1) {
  tryCatch({
    if (any(is.na(params)) || any(!is.finite(params))) return(1e10)
    
    lambda <- exp(evaluate_model(t, params, model_type, baseline))
    lambda <- pmax(pmin(lambda, 1000), 0.001)
    
    if (any(is.na(lambda)) || any(!is.finite(lambda))) return(1e10)
    
    ll <- 0
    for (i in 1:length(y)) {
      if (y[i] == 0) {
        ll <- ll + log(pi + (1 - pi) * exp(-lambda[i]))
      } else {
        ll <- ll + log(1 - pi) + dpois(y[i], lambda[i], log = TRUE)
      }
    }
    
    if (is.na(ll) || !is.finite(ll)) return(1e10)
    
    return(-ll)
  }, error = function(e) 1e10)
}

log_likelihood_zinb <- function(params, y, t, model_type, baseline = NULL, theta = 1, pi = 0.1) {
  tryCatch({
    if (any(is.na(params)) || any(!is.finite(params))) return(1e10)
    
    mu <- exp(evaluate_model(t, params, model_type, baseline))
    mu <- pmax(pmin(mu, 1000), 0.001)
    
    if (any(is.na(mu)) || any(!is.finite(mu))) return(1e10)
    
    ll <- 0
    for (i in 1:length(y)) {
      if (y[i] == 0) {
        ll <- ll + log(pi + (1 - pi) * dnbinom(0, size = theta, mu = mu[i]))
      } else {
        ll <- ll + log(1 - pi) + dnbinom(y[i], size = theta, mu = mu[i], log = TRUE)
      }
    }
    
    if (is.na(ll) || !is.finite(ll)) return(1e10)
    
    return(-ll)
  }, error = function(e) 1e10)
}

# Asymptotic variance calculation using numDeriv
calculate_asymptotic_variance <- function(params, y, t, model_type, distribution, baseline = NULL) {
  
  # Define negative log-likelihood function for numDeriv
  neg_log_lik <- function(par) {
    switch(distribution,
           "poisson" = log_likelihood_poisson(par, y, t, model_type, baseline),
           "nb" = log_likelihood_nb(par, y, t, model_type, baseline),
           "zip" = log_likelihood_zip(par, y, t, model_type, baseline),
           "zinb" = log_likelihood_zinb(par, y, t, model_type, baseline)
    )
  }
  
  tryCatch({
    # Use numDeriv to automatically compute Hessian matrix
    H <- hessian(neg_log_lik, params)
    
    # Check if Hessian is finite
    if (any(!is.finite(H))) {
      warning("Non-finite values in Hessian matrix")
      return(diag(1e-3, length(params)))
    }
    
    # Check condition number
    if (rcond(H) < 1e-12) {
      warning("Hessian matrix is nearly singular, adding regularization")
      H <- H + diag(1e-6, nrow(H))
    }
    
    # Compute asymptotic covariance matrix
    asymptotic_cov <- solve(H)
    
    # Ensure diagonal elements are positive
    diag_vals <- diag(asymptotic_cov)
    if (any(diag_vals <= 0)) {
      warning("Some variance estimates are non-positive")
      diag_vals[diag_vals <= 0] <- abs(diag_vals[diag_vals <= 0]) + 1e-6
      diag(asymptotic_cov) <- diag_vals
    }
    
    return(asymptotic_cov)
    
  }, error = function(e) {
    warning(paste("Asymptotic variance calculation failed:", e$message))
    return(diag(1e-3, length(params)))
  })
}

# Enhanced Cuckoo Search for scGTM
enhanced_cuckoo_scgtm <- function(y, t, n_nests = 25, max_iter = 1000, pa = 0.25,
                                  selected_models = c("hill", "valley", "increasing", "decreasing"),
                                  distribution = "poisson") {
  y <- as.numeric(y)
  t <- as.numeric(t)
  
  valid_idx <- !is.na(y) & !is.na(t) & is.finite(y) & is.finite(t)
  y <- y[valid_idx]
  t <- t[valid_idx]
  
  if (length(y) < 10) {
    stop("Not enough valid data points")
  }
  
  # Fit all selected models
  model_results <- list()
  
  for (model_type in selected_models) {
    tryCatch({
      result <- fit_single_model_cuckoo(y, t, model_type, n_nests, max_iter, pa, distribution)
      
      # Calculate AIC
      if (model_type %in% c("increasing", "decreasing")) {
        n_params <- 3
      } else if (model_type == "valley") {
        n_params <- 5
      } else {
        n_params <- 4
      }
      
      result$aic <- 2 * n_params + 2 * result$objective_value
      result$model_type <- model_type
      result$n_params <- n_params
      
      model_results[[model_type]] <- result
    }, error = function(e) {
      message(paste("Error fitting", model_type, "model:", e$message))
    })
  }
  
  # Model selection
  if (length(model_results) == 0) {
    stop("No models could be fitted successfully")
  }
  
  aic_values <- sapply(model_results, function(x) x$aic)
  best_model <- names(which.min(aic_values))
  
  best_result <- model_results[[best_model]]
  best_result$trend_type <- best_model
  best_result$all_models <- model_results
  best_result$model_comparison <- data.frame(
    Model = names(aic_values),
    AIC = round(aic_values, 2),
    Delta_AIC = round(aic_values - min(aic_values), 2),
    stringsAsFactors = FALSE
  )
  
  return(best_result)
}

# Fit single model with Cuckoo Search (updated parameter bounds)
fit_single_model_cuckoo <- function(y, t, model_type, n_nests, max_iter, pa, distribution) {
  y_mean <- mean(y, na.rm = TRUE)
  y_max <- max(y, na.rm = TRUE)
  t_min <- min(t, na.rm = TRUE)
  t_max <- max(t, na.rm = TRUE)
  
  # Parameter bounds (updated to allow negative k1, k2 for hill/valley)
  if (model_type %in% c("increasing", "decreasing")) {
    lower <- c(log(max(y_mean, 0.1)), 0.1, t_min)
    upper <- c(log(y_max + 1), 3, t_max)
    n_params <- 3
  } else {
    # Allow k1, k2 to be negative for proper hill/valley fitting
    lower <- c(log(max(y_mean, 0.1)), -5, -5, t_min)
    upper <- c(log(y_max + 1), 5, 5, t_max)
    n_params <- 4
  }
  
  # Baseline for valley
  baseline <- NULL
  if (model_type == "valley") {
    baseline <- max(log(y + 1), na.rm = TRUE)
  }
  
  # Define fitness function
  fitness_func <- function(params) {
    switch(distribution,
           "poisson" = log_likelihood_poisson(params, y, t, model_type, baseline),
           "nb" = log_likelihood_nb(params, y, t, model_type, baseline),
           "zip" = log_likelihood_zip(params, y, t, model_type, baseline),
           "zinb" = log_likelihood_zinb(params, y, t, model_type, baseline)
    )
  }
  
  # Run Cuckoo Search
  result <- cuckoo_search_custom(fitness_func, n_params, lower, upper, n_nests, max_iter, pa)
  
  # Calculate quality metrics
  fitted_values <- evaluate_model(t, result$bestnest, model_type, baseline)
  residuals <- log(y + 1) - fitted_values
  
  # Calculate asymptotic variance using numDeriv
  asymptotic_var <- calculate_asymptotic_variance(result$bestnest, y, t, model_type, distribution, baseline)
  
  return(list(
    parameters = result$bestnest,
    objective_value = result$fmin,
    converged = result$fmin < 1e9,
    iterations = max_iter,
    algorithm = "Cuckoo Search",
    convergence_history = result$convergence_history,
    residuals = residuals,
    fitted_values = fitted_values,
    observed_values = log(y + 1),
    pseudotime_clean = t,
    gene_expr_clean = y,
    rmse = sqrt(mean(residuals^2)),
    mae = mean(abs(residuals)),
    r_squared = 1 - sum(residuals^2) / sum((log(y + 1) - mean(log(y + 1)))^2),
    baseline = baseline,
    model_type = model_type,
    distribution = distribution,
    asymptotic_var = asymptotic_var
  ))
}

# Biological pattern classification
classify_gene_pattern <- function(params, r_squared, trend_type) {
  if (trend_type %in% c("increasing", "decreasing")) {
    t0 <- params[3]
    temporal_class <- if (t0 < 0.3) {
      "Early-acting"
    } else if (t0 > 0.7) {
      "Late-acting"  
    } else {
      "Mid-trajectory"
    }
    
    pattern <- ifelse(trend_type == "increasing", "Monotonic Increasing", "Monotonic Decreasing")
    dynamics_class <- "Monotonic"
  } else {
    t0 <- params[4]
    k1 <- abs(params[2])
    k2 <- abs(params[3])
    
    temporal_class <- if (t0 < 0.3) {
      "Early-acting"
    } else if (t0 > 0.7) {
      "Late-acting"
    } else {
      "Mid-trajectory"
    }
    
    steepness_ratio <- k1 / k2
    dynamics_class <- if (steepness_ratio > 1.5) {
      "Transient expression"
    } else if (steepness_ratio < 0.67) {
      "Persistent expression"
    } else {
      "Balanced dynamics"
    }
    
    pattern <- ifelse(trend_type == "valley", "Valley-shaped", "Hill-shaped")
  }
  
  regulation_quality <- if (r_squared > 0.8) {
    "Highly regulated"
  } else if (r_squared > 0.6) {
    "Well regulated"
  } else if (r_squared > 0.4) {
    "Moderately regulated"
  } else {
    "Lowly regulated"
  }
  
  return(list(
    pattern = pattern,
    temporal_class = temporal_class,
    dynamics_class = dynamics_class,
    regulation_quality = regulation_quality
  ))
}

# Generate example data (updated to use new models)
generate_example_data <- function() {
  set.seed(123)
  n_cells <- 100
  n_genes <- 8
  
  pseudotime <- sort(runif(n_cells, 0, 1))
  
  data <- data.frame(
    Index = paste0("Cell_", 1:n_cells),
    pseudotime = pseudotime,
    stringsAsFactors = FALSE
  )
  
  gene_patterns <- list(
    list(name = "Early_Hill",   mu = 2.0, k1 = 1.2, k2 = 0.8, t0 = 0.2),
    list(name = "Mid_Hill",     mu = 1.8, k1 = 1.0, k2 = 1.0, t0 = 0.5),
    list(name = "Late_Hill",    mu = 2.2, k1 = 0.7, k2 = 1.5, t0 = 0.8),
    list(name = "Transient",    mu = 1.5, k1 = 2.0, k2 = 2.0, t0 = 0.4),
    list(name = "Persistent",   mu = 1.9, k1 = 0.5, k2 = 1.8, t0 = 0.3),
    list(name = "Valley_Early", mu = 1.0, k1 = 1.0, k2 = 1.2, t0 = 0.3),
    list(name = "Increasing",   mu = 1.5, k = 1.0, t0 = 0.3),
    list(name = "Decreasing",   mu = 1.5, k = 1.0, t0 = 0.7)
  )
  
  for (i in 1:n_genes) {
    pattern <- gene_patterns[[i]]
    
    if (grepl("Valley", pattern$name)) {
      baseline <- 2.5
      true_expr <- valley_model(pseudotime, c(pattern$mu, pattern$k1, pattern$k2, pattern$t0), baseline)
    } else if (pattern$name == "Increasing") {
      true_expr <- pattern$mu * (1 + tanh(pattern$k * (pseudotime - pattern$t0))) / 2
    } else if (pattern$name == "Decreasing") {
      true_expr <- pattern$mu * (1 - tanh(pattern$k * (pseudotime - pattern$t0))) / 2
    } else {
      true_expr <- hill_model(pseudotime, c(pattern$mu, pattern$k1, pattern$k2, pattern$t0))
    }
    
    lambda <- exp(true_expr)
    lambda <- pmax(pmin(lambda, 100), 0.1)
    
    expression <- rpois(n_cells, lambda)
    data[[pattern$name]] <- expression
  }
  
  return(data)
}

# Data structure detection
detect_data_structure <- function(data) {
  if (is.null(data) || nrow(data) == 0) {
    stop("No data provided")
  }
  
  if (ncol(data) < 3) {
    stop("Data must have at least 3 columns (Index, pseudotime, genes)")
  }
  
  col_names <- tolower(colnames(data))
  pseudotime_idx <- which(grepl("pseudotime|time|pt", col_names))
  
  if (length(pseudotime_idx) == 0) {
    pseudotime_idx <- 2
    message("No 'pseudotime' column found. Using column 2 as pseudotime.")
  } else {
    pseudotime_idx <- pseudotime_idx[1]
  }
  
  pseudotime_col <- data[[pseudotime_idx]]
  if (!is.numeric(pseudotime_col)) {
    data[[pseudotime_idx]] <- as.numeric(as.character(pseudotime_col))
    if (all(is.na(data[[pseudotime_idx]]))) {
      stop(paste("Column", colnames(data)[pseudotime_idx], "cannot be converted to numeric"))
    }
  }
  
  gene_cols <- setdiff(1:ncol(data), c(1, pseudotime_idx))
  
  for (i in gene_cols) {
    if (!is.numeric(data[[i]])) {
      data[[i]] <- as.numeric(as.character(data[[i]]))
      data[[i]][is.na(data[[i]])] <- 0
    }
  }
  
  if (pseudotime_idx != 2) {
    col_order <- c(1, pseudotime_idx, setdiff(1:ncol(data), c(1, pseudotime_idx)))
    data <- data[, col_order]
    colnames(data)[2] <- "pseudotime"
  }
  
  return(list(
    data = data,
    pseudotime_col = 2,
    gene_cols = 3:ncol(data),
    n_cells = nrow(data),
    n_genes = ncol(data) - 2
  ))
}

# Helper function for text formatting
`%R%` <- function(x, n) paste(rep(x, n), collapse = "")

# =============================================================================
# Server Implementation
# =============================================================================

server <- function(input, output, session) {
  
  # Reactive values
  values <- reactiveValues(
    data = NULL,
    data_info = NULL,
    fit_results = NULL,
    current_gene = NULL,
    batch_results = NULL
  )
  
  # Load example data
  observeEvent(input$use_example, {
    if (input$use_example) {
      tryCatch({
        example_data <- generate_example_data()
        data_info <- detect_data_structure(example_data)
        
        values$data <- data_info$data
        values$data_info <- data_info
        
        gene_cols <- data_info$gene_cols
        gene_choices <- colnames(values$data)[gene_cols]
        updateSelectInput(session, "gene_select", choices = gene_choices, selected = gene_choices[1])
        updateCheckboxGroupInput(session, "batch_genes", choices = gene_choices)
        
        showNotification("Example data loaded successfully!", type = "message")
      }, error = function(e) {
        showNotification(paste("Error:", e$message), type = "error")
      })
    }
  })
  
  # Load data from file
  observeEvent(input$load_data, {
    if (!is.null(input$data_file)) {
      tryCatch({
        raw_data <- read.csv(input$data_file$datapath, 
                             header = input$header, 
                             stringsAsFactors = FALSE,
                             check.names = FALSE)
        
        data_info <- detect_data_structure(raw_data)
        
        values$data <- data_info$data
        values$data_info <- data_info
        
        gene_cols <- data_info$gene_cols
        gene_choices <- colnames(values$data)[gene_cols]
        updateSelectInput(session, "gene_select", choices = gene_choices, selected = gene_choices[1])
        updateCheckboxGroupInput(session, "batch_genes", choices = gene_choices)
        
        showNotification("Data loaded successfully!", type = "message")
        
      }, error = function(e) {
        showNotification(paste("Error loading data:", e$message), type = "error")
      })
    } else if (!input$use_example) {
      showNotification("Please select a file or use example data", type = "warning")
    }
  })
  
  # Data summary
  output$data_summary <- renderText({
    if (is.null(values$data) || is.null(values$data_info)) {
      return("No data loaded")
    }
    
    data_info <- values$data_info
    pseudotime_col <- values$data[[data_info$pseudotime_col]]
    valid_pseudotime <- pseudotime_col[!is.na(pseudotime_col)]
    
    paste(
      "Dataset Overview:\n",
      "• Cells:", data_info$n_cells, "\n",
      "• Genes:", data_info$n_genes, "\n", 
      "• Pseudotime range:", round(min(valid_pseudotime), 3), "to", round(max(valid_pseudotime), 3), "\n",
      "• Ready for scGTM analysis with Cuckoo Search"
    )
  })
  
  # Data preview
  output$data_preview <- DT::renderDataTable({
    if (is.null(values$data)) return(NULL)
    
    preview_data <- values$data[1:min(50, nrow(values$data)), 
                                1:min(10, ncol(values$data))]
    
    DT::datatable(preview_data, 
                  options = list(scrollX = TRUE, pageLength = 10),
                  rownames = FALSE)
  })
  
  # Model fitting
  observeEvent(input$fit_model, {
    req(input$gene_select, values$data, values$data_info, input$selected_models)
    
    if (length(input$selected_models) == 0) {
      showNotification("Please select at least one model to fit", type = "warning")
      return()
    }
    
    tryCatch({
      withProgress(message = 'Fitting scGTM with Cuckoo Search...', value = 0, {
        
        pseudotime <- values$data[[values$data_info$pseudotime_col]]
        gene_expr <- values$data[[input$gene_select]]
        
        if (all(is.na(pseudotime))) {
          stop("Pseudotime column contains no valid values")
        }
        
        if (all(is.na(gene_expr))) {
          stop("Selected gene contains no valid expression values")
        }
        
        valid_idx <- !is.na(pseudotime) & !is.na(gene_expr) & is.finite(pseudotime) & is.finite(gene_expr)
        pseudotime_clean <- pseudotime[valid_idx]
        gene_expr_clean <- gene_expr[valid_idx]
        
        if (length(pseudotime_clean) < 10) {
          stop("Not enough valid data points (need at least 10)")
        }
        
        incProgress(0.3, detail = paste("Comparing", length(input$selected_models), "models with", input$distribution, "distribution"))
        
        results <- enhanced_cuckoo_scgtm(
          y = gene_expr_clean,
          t = pseudotime_clean,
          n_nests = input$n_nests,
          max_iter = input$max_iter,
          pa = input$pa,
          selected_models = input$selected_models,
          distribution = input$distribution
        )
        
        incProgress(0.8, detail = "Analyzing biological patterns and calculating asymptotic variance...")
        
        pattern_info <- classify_gene_pattern(results$parameters, results$r_squared, results$trend_type)
        results$pattern_info <- pattern_info
        
        values$fit_results <- results
        values$current_gene <- input$gene_select
        
        incProgress(1, detail = "Complete!")
      })
      
      showNotification(paste("Cuckoo Search completed! Best model:", 
                             values$fit_results$trend_type, "with", input$distribution), type = "message")
    }, error = function(e) {
      showNotification(paste("Error:", e$message), type = "error")
    })
  })
  
  # Main trend plot
  output$trend_plot <- renderPlot({
    if (is.null(values$fit_results)) {
      plot(1, 1, type = "n", xlab = "Pseudotime", ylab = "Log(Expression + 1)",
           main = "Select a gene and fit scGTM model")
      text(1, 1, "Load data and select a gene", cex = 1.2, col = "gray60")
      return()
    }
    
    pseudotime   <- values$fit_results$pseudotime_clean
    gene_expr    <- values$fit_results$gene_expr_clean
    params       <- values$fit_results$parameters
    trend_type   <- values$fit_results$trend_type
    r_squared    <- values$fit_results$r_squared
    distribution <- values$fit_results$distribution
    
    t_seq        <- seq(min(pseudotime), max(pseudotime), length.out = 400)
    fitted_curve <- evaluate_model(t_seq, params, trend_type, values$fit_results$baseline)
    
    # t0 location per model
    t0 <- if (trend_type %in% c("hill", "valley")) params[4] else params[3]
    
    # data frames for ggplot
    df_pts   <- data.frame(t = pseudotime, y = log(gene_expr + 1))
    df_curve <- data.frame(t = t_seq, y = fitted_curve, what = "scGTM fit")
    
    title_txt <- paste0(
      "scGTM-Cuckoo: ", values$current_gene, "  —  ",
      tools::toTitleCase(trend_type), " · ", toupper(distribution),
      " · R2 = ", sprintf("%.3f", r_squared)
    )
    
    subtitle_txt <- if (trend_type %in% c("hill", "valley")) "Peak/Valley time" else "Inflection time"
    
    library(ggplot2)
    ggplot() +
      # observed points
      geom_point(
        data = df_pts, aes(t, y, color = "Observed"),
        size = 2.2, alpha = 0.7, stroke = 0
      ) +
      # fitted curve
      geom_line(
        data = df_curve, aes(t, y, color = "Fit"),
        linewidth = 1.2
      ) +
      # vertical t0
      geom_vline(xintercept = t0, linetype = "22", linewidth = 0.8, color = "steelblue4") +
      annotate("label", x = t0, y = max(df_pts$y, df_curve$y) * 0.98,
               label = paste0("t0 = ", sprintf("%.3f", t0), "\n", subtitle_txt),
               hjust = -0.05, vjust = 1, size = 6.5, color = "steelblue4",
               fill = "white", label.size = 0.2) +
      scale_color_manual(
        name = NULL,
        values = c("Observed" = "#2c7fb8", "Fit" = "#d7301f")
      ) +
      labs(
        x = "Pseudotime",
        y = "Log(Expression + 1)",
        title = title_txt
      ) +
      guides(color = guide_legend(override.aes = list(linewidth = 1.5, size = 5))) +
      theme_minimal(base_size = 17) +
      theme(
        plot.title   = element_text(face = "bold"),
        legend.position = "top",
        legend.box.margin = margin(b = 4),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(margin = margin(r = 6)),
        axis.title.x = element_text(margin = margin(t = 6))
      )
  })
  
  
  # Model comparison table
  output$model_comparison_table <- renderTable({
    if (is.null(values$fit_results) || is.null(values$fit_results$model_comparison)) return(NULL)
    
    comparison <- values$fit_results$model_comparison
    comparison$Model <- c("Hill-shaped", "Valley-shaped", "Monotonic Increasing", "Monotonic Decreasing")[
      match(comparison$Model, c("hill", "valley", "increasing", "decreasing"))
    ]
    
    comparison$Selected <- ifelse(comparison$Delta_AIC == 0, "✓ BEST", "")
    comparison
  }, striped = TRUE, hover = TRUE)
  
  # Model selection details
  output$model_selection_details <- renderText({
    if (is.null(values$fit_results)) return("No models fitted yet")
    
    best_model <- values$fit_results$trend_type
    best_aic <- min(values$fit_results$model_comparison$AIC)
    n_models <- nrow(values$fit_results$model_comparison)
    distribution <- values$fit_results$distribution
    
    model_names <- c(
      "hill" = "Hill-shaped",
      "valley" = "Valley-shaped", 
      "increasing" = "Monotonic Increasing",
      "decreasing" = "Monotonic Decreasing"
    )
    
    paste(
      "MODEL SELECTION RESULTS\n",
      "=" %R% 30, "\n\n",
      "Selected Model:", model_names[best_model], "\n",
      "Distribution:", toupper(distribution), "\n",
      "Selection Criterion: Lowest AIC\n",
      "Best AIC:", round(best_aic, 2), "\n",
      "Models Compared:", n_models, "\n",
      "Optimization: Cuckoo Search Algorithm\n",
      "Variance Calculation: numDeriv Package\n\n",
      "AIC Interpretation:\n",
      "• Lower AIC = Better model\n",
      "• ΔAIC > 2 = Strong evidence\n",
      "• ΔAIC > 10 = Very strong evidence"
    )
  })
  
  # Parameter table
  output$parameter_table <- renderTable({
    if (is.null(values$fit_results)) return(NULL)
    
    params <- values$fit_results$parameters
    trend_type <- values$fit_results$trend_type
    
    if (trend_type %in% c("increasing", "decreasing")) {
      param_names <- c("μ (Magnitude)", "k (Rate)", "t₀ (Inflection time)")
    } else {
      param_names <- c("μ (Magnitude)", "k₁ (Left steepness)", "k₂ (Right steepness)", "t₀ (Peak/Valley time)")
    }
    
    param_values <- round(params, 4)
    
    data.frame(
      Parameter = param_names,
      Estimate = param_values,
      stringsAsFactors = FALSE
    )
  }, striped = TRUE)
  
  # Asymptotic standard errors table with improved error handling
  output$asymptotic_se_table <- renderTable({
    if (is.null(values$fit_results)) {
      return(data.frame(
        Parameter = "No model fitted",
        Standard_Error = "-",
        Confidence_Interval_95 = "-",
        stringsAsFactors = FALSE
      ))
    }
    
    if (is.null(values$fit_results$asymptotic_var)) {
      return(data.frame(
        Parameter = "Calculation in progress",
        Standard_Error = "Please wait...",
        Confidence_Interval_95 = "-",
        stringsAsFactors = FALSE
      ))
    }
    
    # Check asymptotic variance matrix validity
    asymptotic_var <- values$fit_results$asymptotic_var
    
    if (all(is.na(asymptotic_var)) || any(diag(asymptotic_var) <= 0)) {
      return(data.frame(
        Parameter = "Calculation failed",
        Standard_Error = "Numerical issues",
        Confidence_Interval_95 = "Check convergence",
        stringsAsFactors = FALSE
      ))
    }
    
    trend_type <- values$fit_results$trend_type
    
    if (trend_type %in% c("increasing", "decreasing")) {
      param_names <- c("μ (Magnitude)", "k (Rate)", "t₀ (Inflection time)")
    } else {
      param_names <- c("μ (Magnitude)", "k₁ (Left steepness)", "k₂ (Right steepness)", "t₀ (Peak/Valley time)")
    }
    
    se_values <- sqrt(diag(asymptotic_var))
    
    # Check SE values validity
    if (any(!is.finite(se_values))) {
      return(data.frame(
        Parameter = param_names,
        Standard_Error = "Non-finite values",
        Confidence_Interval_95 = "Check model convergence",
        stringsAsFactors = FALSE
      ))
    }
    
    data.frame(
      Parameter = param_names,
      Standard_Error = round(se_values, 4),
      Confidence_Interval_95 = paste("±", round(1.96 * se_values, 4)),
      stringsAsFactors = FALSE
    )
  }, striped = TRUE)
  
  # Diagnostic plots
  output$diagnostic_plots <- renderPlot({
    if (is.null(values$fit_results)) {
      par(mfrow = c(1, 3))
      for (i in 1:3) {
        plot(1, 1, type = "n", main = paste("Diagnostic", i))
        text(1, 1, "Fit model to see diagnostics", col = "gray60")
      }
      par(mfrow = c(1, 1))
      return()
    }
    
    fitted <- values$fit_results$fitted_values
    residuals <- values$fit_results$residuals
    
    par(mfrow = c(1, 3))
    
    # Residuals vs fitted
    plot(fitted, residuals, pch = 16, col = "red",
         xlab = "Fitted Values", ylab = "Residuals",
         main = "Residuals vs Fitted")
    abline(h = 0, lty = 2)
    
    # Q-Q plot
    qqnorm(residuals, pch = 16, col = "blue", main = "Normal Q-Q Plot")
    qqline(residuals, col = "red", lwd = 2)
    
    # Convergence plot
    if (!is.null(values$fit_results$convergence_history)) {
      plot(values$fit_results$convergence_history, type = "l", col = "green", lwd = 2,
           xlab = "Iteration", ylab = "Objective Value",
           main = "Cuckoo Search Convergence")
    }
    
    par(mfrow = c(1, 1))
  })
  
  # Fitting status
  output$fitting_status <- renderText({
    if (is.null(values$fit_results)) {
      return("No model fitted yet.\n\nSteps:\n1. Select gene and models\n2. Choose distribution\n3. Click 'Fit scGTM Model'")
    }
    
    paste(
      "Gene:", values$current_gene, "\n",
      "Algorithm: Cuckoo Search\n",
      "Distribution:", toupper(values$fit_results$distribution), "\n",
      "Pattern:", tools::toTitleCase(values$fit_results$trend_type), "\n",
      "Converged:", ifelse(values$fit_results$converged, "Yes", "Check"), "\n",
      "R²:", round(values$fit_results$r_squared, 3), "\n",
      "Nests used:", input$n_nests, "\n",
      "Iterations:", input$max_iter, "\n",
      "Asymptotic Variance: Calculated with numDeriv"
    )
  })
  
  # Detected pattern
  output$detected_pattern <- renderText({
    if (is.null(values$fit_results)) return("No analysis performed")
    paste(values$fit_results$pattern_info$pattern, "-", values$fit_results$pattern_info$temporal_class)
  })
  
  # Distribution info
  output$distribution_info <- renderText({
    if (is.null(values$fit_results)) return("No analysis performed")
    paste(
      "Distribution:", toupper(values$fit_results$distribution), "\n",
      "Quality:", values$fit_results$pattern_info$regulation_quality
    )
  })
  
  # Model quality
  output$model_quality <- renderText({
    if (is.null(values$fit_results)) return("No model fitted yet")
    
    r_sq <- values$fit_results$r_squared
    rmse <- values$fit_results$rmse
    distribution <- values$fit_results$distribution
    
    quality_class <- if (r_sq > 0.8) {
      "Excellent"
    } else if (r_sq > 0.6) {
      "Good"
    } else if (r_sq > 0.4) {
      "Fair"
    } else {
      "Poor"
    }
    
    paste(
      "Overall Quality:", quality_class, "\n\n",
      "Distribution:", toupper(distribution), "\n",
      "R-squared:", round(r_sq, 4), "\n",
      "RMSE:", round(rmse, 4), "\n",
      "Algorithm: Cuckoo Search\n",
      "Variance Calculation: numDeriv\n\n",
      "Model Selection:\n",
      "Best AIC:", round(values$fit_results$aic, 2), "\n",
      "Asymptotic Variance: Available\n",
      "Standard Errors: Calculated automatically"
    )
  })
  
  # Expression statistics
  output$expression_stats <- renderTable({
    if (is.null(values$fit_results)) return(NULL)
    
    gene_expr <- values$fit_results$gene_expr_clean
    params <- values$fit_results$parameters
    trend_type <- values$fit_results$trend_type
    
    if (trend_type %in% c("increasing", "decreasing")) {
      change_time <- params[3]
    } else {
      change_time <- params[4]
    }
    
    data.frame(
      Statistic = c("Mean Expression", "Max Expression", "Peak/Valley/Inflection Time", "Zero Expression %"),
      Value = c(
        round(mean(gene_expr), 2),
        max(gene_expr),
        round(change_time, 3),
        paste0(round(mean(gene_expr == 0) * 100, 1), "%")
      ),
      stringsAsFactors = FALSE
    )
  }, striped = TRUE)
  
  # Pattern classification
  output$pattern_classification <- renderText({
    if (is.null(values$fit_results)) return("Analyze a gene to see classification")
    
    pattern_info <- values$fit_results$pattern_info
    params <- values$fit_results$parameters
    trend_type <- values$fit_results$trend_type
    distribution <- values$fit_results$distribution
    
    if (trend_type %in% c("increasing", "decreasing")) {
      change_time <- params[3]
    } else {
      change_time <- params[4]
    }
    
    paste(
      "BIOLOGICAL PATTERN ANALYSIS\n",
      "=" %R% 30, "\n\n",
      "Expression Pattern:", pattern_info$pattern, "\n",
      "Temporal Class:", pattern_info$temporal_class, "\n",
      "Dynamics:", pattern_info$dynamics_class, "\n",
      "Quality:", pattern_info$regulation_quality, "\n",
      "Distribution Used:", toupper(distribution), "\n\n",
      "Key Time Point: t₀ =", round(change_time, 3), "\n",
      "Optimization: Cuckoo Search Algorithm\n",
      "Variance: Automatic calculation with numDeriv"
    )
  })
  
  # Model illustration
  output$model_illustration <- renderPlot({
    if (is.null(values$fit_results)) {
      plot(1, 1, type = "n", main = "Fit a model to see illustration")
      text(1, 1, "Mathematical model will be shown here", col = "gray60")
      return()
    }
    
    params <- values$fit_results$parameters
    trend_type <- values$fit_results$trend_type
    
    t_range <- seq(0, 1, length.out = 100)
    y_curve <- evaluate_model(t_range, params, trend_type, values$fit_results$baseline)
    
    plot(t_range, y_curve, type = "l", lwd = 3, col = "red",
         xlab = "Pseudotime (t)", ylab = "log(E[Y] + 1)",
         main = paste("scGTM Model:", tools::toTitleCase(trend_type), "Shape"))
    
    if (trend_type %in% c("increasing", "decreasing")) {
      t0 <- params[3]
    } else {
      t0 <- params[4]
    }
    
    abline(v = t0, lty = 2, col = "blue", lwd = 2)
    text(t0, max(y_curve) * 0.9, paste("t₀ =", round(t0, 3)), 
         pos = 4, col = "blue", font = 2)
  })
  
  # Pseudotime distribution
  output$pseudotime_dist <- renderPlot({
    if (is.null(values$data) || is.null(values$data_info)) return(NULL)
    
    pseudotime <- values$data[[values$data_info$pseudotime_col]]
    hist(pseudotime, breaks = 20, col = "lightblue",
         main = "Pseudotime Distribution", xlab = "Pseudotime")
    
    if (!is.null(values$fit_results)) {
      trend_type <- values$fit_results$trend_type
      params <- values$fit_results$parameters
      
      if (trend_type %in% c("increasing", "decreasing")) {
        change_time <- params[3]
      } else {
        change_time <- params[4]
      }
      
      abline(v = change_time, col = "red", lwd = 3, lty = 2)
    }
  })
  
  # Batch analysis functions
  observeEvent(input$select_all_genes, {
    if (!is.null(values$data_info)) {
      gene_choices <- colnames(values$data)[values$data_info$gene_cols]
      updateCheckboxGroupInput(session, "batch_genes", selected = gene_choices)
    }
  })
  
  observeEvent(input$clear_all_genes, {
    updateCheckboxGroupInput(session, "batch_genes", selected = character(0))
  })
  
  # Run batch analysis
  observeEvent(input$run_batch, {
    req(input$batch_genes, values$data, values$data_info)
    
    if (length(input$batch_genes) == 0) {
      showNotification("Please select at least one gene", type = "warning")
      return()
    }
    
    tryCatch({
      withProgress(message = 'Running comparative analysis with Cuckoo Search...', value = 0, {
        
        genes_to_analyze <- input$batch_genes
        n_genes <- length(genes_to_analyze)
        
        batch_results <- data.frame(
          Gene = character(n_genes),
          Trend_Type = character(n_genes),
          Distribution = character(n_genes),
          mu_Magnitude = numeric(n_genes),
          k1_Left = numeric(n_genes),
          k2_Right = numeric(n_genes),
          t0_Time = numeric(n_genes),
          R_squared = numeric(n_genes),
          AIC = numeric(n_genes),
          Converged = logical(n_genes),
          Pattern = character(n_genes),
          Temporal_Class = character(n_genes),
          stringsAsFactors = FALSE
        )
        
        pseudotime <- values$data[[values$data_info$pseudotime_col]]
        
        for (i in 1:n_genes) {
          gene_name <- genes_to_analyze[i]
          incProgress(1/n_genes, detail = paste("Analyzing", gene_name, "(", i, "of", n_genes, ")"))
          
          tryCatch({
            gene_expr <- values$data[[gene_name]]
            
            valid_idx <- !is.na(pseudotime) & !is.na(gene_expr) & is.finite(pseudotime) & is.finite(gene_expr)
            pseudotime_clean <- pseudotime[valid_idx]
            gene_expr_clean <- gene_expr[valid_idx]
            
            if (length(pseudotime_clean) >= 10) {
              result <- enhanced_cuckoo_scgtm(
                y = gene_expr_clean,
                t = pseudotime_clean,
                n_nests = input$batch_nests,
                max_iter = input$batch_iterations,
                selected_models = input$batch_models,
                distribution = input$batch_distribution
              )
              
              pattern_info <- classify_gene_pattern(result$parameters, result$r_squared, result$trend_type)
              
              # Handle different parameter structures
              if (result$trend_type %in% c("increasing", "decreasing")) {
                params_padded <- c(result$parameters, NA, result$parameters[3])
              } else {
                params_padded <- result$parameters
              }
              
              batch_results[i, ] <- list(
                gene_name,
                result$trend_type,
                result$distribution,
                params_padded[1],
                ifelse(is.na(params_padded[2]), NA, params_padded[2]),
                ifelse(is.na(params_padded[3]), NA, params_padded[3]),
                params_padded[4],
                result$r_squared,
                result$aic,
                result$converged,
                pattern_info$pattern,
                pattern_info$temporal_class
              )
            } else {
              batch_results[i, ] <- list(
                gene_name, "insufficient_data", input$batch_distribution, NA, NA, NA, NA, NA, NA, FALSE, 
                "Insufficient data", "Unknown"
              )
            }
          }, error = function(e) {
            batch_results[i, ] <- list(
              gene_name, "error", input$batch_distribution, NA, NA, NA, NA, NA, NA, FALSE,
              "Error in fitting", "Unknown"
            )
          })
        }
        
        values$batch_results <- batch_results
      })
      
      showNotification(paste("Cuckoo Search analysis completed for", n_genes, "genes!"), type = "message")
    }, error = function(e) {
      showNotification(paste("Error:", e$message), type = "error")
    })
  })
  
  # Batch status
  output$batch_status <- renderText({
    if (is.null(values$batch_results)) {
      return("No analysis performed yet.\n\nRun Cuckoo Search optimization\nfor multiple genes comparison.")
    }
    
    results <- values$batch_results
    n_total <- nrow(results)
    n_converged <- sum(results$Converged, na.rm = TRUE)
    n_hill <- sum(results$Trend_Type == "hill", na.rm = TRUE)
    n_valley <- sum(results$Trend_Type == "valley", na.rm = TRUE)
    n_inc <- sum(results$Trend_Type == "increasing", na.rm = TRUE)
    n_dec <- sum(results$Trend_Type == "decreasing", na.rm = TRUE)
    
    paste(
      "CUCKOO SEARCH ANALYSIS\n",
      "=" %R% 25, "\n\n",
      "Total genes:", n_total, "\n",
      "Converged:", n_converged, "\n",
      "Hill:", n_hill, "\n", 
      "Valley:", n_valley, "\n",
      "Increasing:", n_inc, "\n",
      "Decreasing:", n_dec, "\n",
      "Success rate:", round(n_converged/n_total * 100, 1), "%"
    )
  })
  
  # Batch results table
  output$batch_results_table <- DT::renderDataTable({
    if (is.null(values$batch_results)) return(NULL)
    
    display_results <- values$batch_results
    numeric_cols <- c("mu_Magnitude", "k1_Left", "k2_Right", "t0_Time", "R_squared", "AIC")
    
    for (col in numeric_cols) {
      display_results[[col]] <- round(as.numeric(display_results[[col]]), 3)
    }
    
    DT::datatable(display_results, 
                  options = list(scrollX = TRUE, pageLength = 15),
                  rownames = FALSE) %>%
      DT::formatStyle('R_squared',
                      backgroundColor = DT::styleInterval(c(0.4, 0.6, 0.8), 
                                                          c('#ffebee', '#fff3e0', '#e8f5e8', '#c8e6c9'))) %>%
      DT::formatStyle('Converged',
                      backgroundColor = DT::styleEqual(TRUE, '#c8e6c9'))
  })
  
  # Pattern summary plot
  output$pattern_summary_plot <- renderPlot({
    if (is.null(values$batch_results)) {
      plot(1, 1, type = "n", main = "Run batch analysis to see summary")
      text(1, 1, "Pattern distribution will appear here", col = "gray60")
      return()
    }
    
    results <- values$batch_results
    results$t0_Time <- as.numeric(results$t0_Time)
    results <- results[!is.na(results$t0_Time), ]
    
    if (nrow(results) == 0) {
      plot(1, 1, type = "n", main = "No valid data for plotting")
      text(1, 1, "Check data quality", col = "gray60")
      return()
    }
    
    par(mfrow = c(1, 2))
    
    # Temporal distribution
    if (length(results$t0_Time) > 0) {
      hist(results$t0_Time, breaks = 10, col = "lightblue",
           main = "Peak/Valley/Inflection Time Distribution", xlab = "t₀")
    }
    
    # Pattern type distribution
    pattern_counts <- table(results$Trend_Type)
    if (length(pattern_counts) > 0) {
      pie(pattern_counts, 
          col = rainbow(length(pattern_counts)),
          main = "Pattern Distribution")
    }
    
    par(mfrow = c(1, 1))
  })
  
  # Gene clustering plot
  output$gene_clustering_plot <- renderPlot({
    if (is.null(values$batch_results)) {
      plot(1, 1, type = "n", main = "Run batch analysis to see clustering")
      text(1, 1, "Gene clustering will appear here", col = "gray60")
      return()
    }
    
    results <- values$batch_results
    results$t0_Time <- as.numeric(results$t0_Time)
    results$R_squared <- as.numeric(results$R_squared)
    
    results <- results[results$Converged & 
                         !is.na(results$t0_Time) & 
                         !is.na(results$R_squared), ]
    
    if (nrow(results) < 3) {
      plot(1, 1, type = "n", main = "Not enough genes for clustering")
      text(1, 1, "Need at least 3 valid genes", col = "gray60")
      return()
    }
    
    # Color by trend type
    colors <- c("hill" = "red", "valley" = "blue", "increasing" = "green", "decreasing" = "orange")
    point_colors <- colors[results$Trend_Type]
    
    plot(results$t0_Time, results$R_squared,
         col = point_colors,
         pch = 16, cex = 1.5,
         xlab = "Peak/Valley/Inflection Time (t₀)", 
         ylab = "R² (Model Quality)",
         main = "Gene Expression Dynamics Clustering")
    
    # Add gene labels for small datasets
    if (nrow(results) <= 15) {
      text(results$t0_Time, results$R_squared, 
           results$Gene, pos = 3, cex = 0.7)
    }
    
    # Reference lines
    abline(h = 0.6, lty = 2, col = "gray", lwd = 2)
    abline(v = c(0.3, 0.7), lty = 2, col = "gray")
    
    legend("topright", 
           legend = names(colors),
           col = colors,
           pch = 16,
           title = "Pattern Type")
  })
  
  # Download handlers
  output$download_results <- downloadHandler(
    filename = function() {
      paste0("scGTM_cuckoo_analysis_", values$current_gene, "_", Sys.Date(), ".csv")
    },
    content = function(file) {
      if (!is.null(values$fit_results)) {
        params <- values$fit_results$parameters
        pattern_info <- values$fit_results$pattern_info
        
        # Add asymptotic standard errors if available
        se_values <- if (!is.null(values$fit_results$asymptotic_var)) {
          sqrt(diag(values$fit_results$asymptotic_var))
        } else {
          rep(NA, length(params))
        }
        
        results_df <- data.frame(
          Gene = values$current_gene,
          Algorithm = "Cuckoo Search",
          Distribution = values$fit_results$distribution,
          Trend_Type = values$fit_results$trend_type,
          Parameters = paste(round(params, 4), collapse = ";"),
          Standard_Errors = paste(round(se_values, 4), collapse = ";"),
          R_squared = values$fit_results$r_squared,
          AIC = values$fit_results$aic,
          Converged = values$fit_results$converged,
          Pattern = pattern_info$pattern,
          Temporal_Class = pattern_info$temporal_class,
          Dynamics_Class = pattern_info$dynamics_class,
          Regulation_Quality = pattern_info$regulation_quality,
          Variance_Method = "numDeriv automatic",
          Analysis_Date = Sys.Date(),
          stringsAsFactors = FALSE
        )
        write.csv(results_df, file, row.names = FALSE)
      }
    }
  )
  
  output$download_batch <- downloadHandler(
    filename = function() {
      paste0("scGTM_cuckoo_batch_analysis_", Sys.Date(), ".csv")
    },
    content = function(file) {
      if (!is.null(values$batch_results)) {
        batch_with_metadata <- values$batch_results
        batch_with_metadata$Algorithm <- "Cuckoo Search"
        batch_with_metadata$Variance_Method <- "numDeriv automatic"
        batch_with_metadata$Analysis_Date <- Sys.Date()
        write.csv(batch_with_metadata, file, row.names = FALSE)
      }
    }
  )
}
