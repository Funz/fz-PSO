#title: Particle Swarm Optimization
#author: Yann Richet, Claus Bendtsen
#type: optimization
#options: yminimization=true;max_iterations=30;nparticles=NA;seed=123;w=0.7213;c_p=1.1931;c_g=1.1931
#options.help: yminimization='Minimize output value ?';max_iterations='Maximum number of iterations';nparticles='Number of particles (default: auto = 10+2*sqrt(d))';seed='Random seed';w='Inertia weight';c_p='Cognitive (personal best) coefficient';c_g='Social (global best) coefficient'

# Constructor for PSO S3 class
PSO <- function(...) {
  # Get options from ... arguments
  opts <- list(...)

  # Create object with initial state
  # Use an environment for mutable state (idiomatic S3 pattern)
  state <- new.env(parent = emptyenv())
  state$i <- 0
  state$input <- NULL
  state$d <- NULL
  state$s <- NULL    # swarm size (set in get_initial_design)
  state$X <- NULL    # current particle positions (matrix: d x s)
  state$V <- NULL    # current particle velocities (matrix: d x s)
  state$P <- NULL    # personal best positions (matrix: d x s)
  state$f_p <- NULL  # personal best fitness values (vector: s)
  state$i_best <- NULL  # index of global best particle
  state$converged <- FALSE

  obj <- list(
    options = list(
      yminimization = isTRUE(as.logical(
        ifelse(is.null(opts$yminimization), TRUE, opts$yminimization)
      )),
      max_iterations = as.integer(
        ifelse(is.null(opts$max_iterations), 30, opts$max_iterations)
      ),
      nparticles = if (is.null(opts$nparticles) || is.na(opts$nparticles)) NA else as.integer(opts$nparticles),
      seed = as.integer(
        ifelse(is.null(opts$seed), 123, opts$seed)
      ),
      w = as.numeric(
        ifelse(is.null(opts$w), 1 / (2 * log(2)), opts$w)
      ),
      c_p = as.numeric(
        ifelse(is.null(opts$c_p), 0.5 + log(2), opts$c_p)
      ),
      c_g = as.numeric(
        ifelse(is.null(opts$c_g), 0.5 + log(2), opts$c_g)
      )
    ),
    state = state  # Environment for mutable state
  )

  # Set S3 class
  class(obj) <- "PSO"

  return(obj)
}

# Generic function definitions (if not already defined)
if (!exists("get_initial_design")) {
  get_initial_design <- function(obj, ...) UseMethod("get_initial_design")
}

if (!exists("get_next_design")) {
  get_next_design <- function(obj, ...) UseMethod("get_next_design")
}

if (!exists("get_analysis")) {
  get_analysis <- function(obj, ...) UseMethod("get_analysis")
}

if (!exists("get_analysis_tmp")) {
  get_analysis_tmp <- function(obj, ...) UseMethod("get_analysis_tmp")
}

# Method: get_initial_design
get_initial_design.PSO <- function(obj, input_variables, output_variables) {
  # Store input variables in mutable state
  obj$state$input <- input_variables
  obj$state$i <- 0
  obj$state$converged <- FALSE

  d <- length(input_variables)
  obj$state$d <- d

  # Set swarm size
  if (is.na(obj$options$nparticles)) {
    obj$state$s <- as.integer(floor(10 + 2 * sqrt(d)))
  } else {
    obj$state$s <- obj$options$nparticles
  }

  set.seed(obj$options$seed)

  s <- obj$state$s

  # Initialize particle positions randomly in [0,1]^d
  X <- matrix(runif(d * s), nrow = d, ncol = s)
  colnames(X) <- NULL

  # Initialize velocities
  X_alt <- matrix(runif(d * s), nrow = d, ncol = s)
  V <- (X_alt - X) / 2

  obj$state$X <- X
  obj$state$V <- V

  # Convert positions to real space and return as list of named lists
  design <- matrix_to_points(X, names(input_variables), input_variables)

  return(design)
}

# Method: get_next_design
get_next_design.PSO <- function(obj, X, Y) {
  # Check max iterations
  if (obj$state$i >= obj$options$max_iterations) {
    obj$state$converged <- FALSE
    return(list())
  }

  s <- obj$state$s
  d <- obj$state$d
  input_names <- names(obj$state$input)

  # Convert Y to numeric vector
  Y_vec <- unlist(Y)

  # Handle minimization/maximization by flipping sign
  if (!obj$options$yminimization) {
    Y_eval <- -Y_vec
  } else {
    Y_eval <- Y_vec
  }

  # Get the latest s evaluations (the most recent batch)
  n <- length(Y_eval)
  f_x <- Y_eval[(n - s + 1):n]

  if (obj$state$i == 0) {
    # First iteration: initialize personal bests
    obj$state$P <- obj$state$X  # personal best positions = initial positions
    obj$state$f_p <- f_x        # personal best fitness = initial fitness
    obj$state$i_best <- which.min(obj$state$f_p)
  } else {
    # Update personal bests
    improved <- f_x < obj$state$f_p
    # Handle NAs
    improved[is.na(improved)] <- FALSE

    if (any(improved)) {
      obj$state$P[, improved] <- obj$state$X[, improved]
      obj$state$f_p[improved] <- f_x[improved]
      obj$state$i_best <- which.min(obj$state$f_p)
    }
  }

  obj$state$i <- obj$state$i + 1

  # Update velocities and positions (SPSO 2007 style)
  w <- obj$options$w
  c_p <- obj$options$c_p
  c_g <- obj$options$c_g
  i_best <- obj$state$i_best

  for (i in 1:s) {
    # Inertia
    obj$state$V[, i] <- w * obj$state$V[, i]

    # Cognitive component (personal best)
    r_p <- runif(d)
    obj$state$V[, i] <- obj$state$V[, i] +
      c_p * r_p * (obj$state$P[, i] - obj$state$X[, i])

    # Social component (global best)
    if (i != i_best) {
      r_g <- runif(d)
      obj$state$V[, i] <- obj$state$V[, i] +
        c_g * r_g * (obj$state$P[, i_best] - obj$state$X[, i])
    }

    # Update position
    obj$state$X[, i] <- obj$state$X[, i] + obj$state$V[, i]

    # Enforce bounds [0, 1] with velocity reset
    below <- obj$state$X[, i] < 0
    above <- obj$state$X[, i] > 1
    if (any(below)) {
      obj$state$X[below, i] <- 0
      obj$state$V[below, i] <- 0
    }
    if (any(above)) {
      obj$state$X[above, i] <- 1
      obj$state$V[above, i] <- 0
    }
  }

  # Convert new positions to real space and return
  design <- matrix_to_points(obj$state$X, input_names, obj$state$input)

  return(design)
}

# Method: get_analysis
get_analysis.PSO <- function(obj, X, Y) {
  analysis_dict <- list(text = "", data = list())

  # Convert Y to numeric
  Y_vec <- unlist(Y)

  # Filter valid values
  valid <- !is.na(Y_vec) & !is.null(Y_vec)
  Y_valid <- Y_vec[valid]

  if (length(Y_valid) < 1) {
    analysis_dict$text <- "No valid results to analyze"
    analysis_dict$data <- list(valid_samples = 0)
    return(analysis_dict)
  }

  # Convert X to matrix
  input_names <- names(obj$state$input)
  X_mat <- points_to_matrix(X, input_names)

  # Find optimum
  if (obj$options$yminimization) {
    m <- min(Y_valid)
    m_ix <- which.min(Y_vec)
  } else {
    m <- max(Y_valid)
    m_ix <- which.max(Y_vec)
  }
  m_ix <- m_ix[1]
  x_opt <- X_mat[m_ix, ]

  opt_type <- if (obj$options$yminimization) "minimum" else "maximum"

  # Store data
  analysis_dict$data <- list(
    optimum = m,
    optimum_point = as.list(x_opt),
    n_evaluations = length(Y_valid),
    iterations = obj$state$i,
    nparticles = obj$state$s
  )

  # Create text summary
  analysis_dict$text <- sprintf(
"PSO Optimization Results:
  Iterations: %d
  Swarm size: %d
  Total evaluations: %d
  %s: %.6f
  Found at: %s
",
    obj$state$i,
    obj$state$s,
    length(Y_valid),
    opt_type,
    m,
    paste(paste(names(x_opt), "=", sprintf("%.6f", x_opt), sep = ""), collapse = "; ")
  )

  # Try to create HTML with plot
  tryCatch({
    # Color by fitness
    if (max(Y_valid) > min(Y_valid)) {
      red <- (Y_valid - min(Y_valid)) / (max(Y_valid) - min(Y_valid))
    } else {
      red <- rep(0.5, length(Y_valid))
    }

    png_file <- tempfile(fileext = ".png")
    png(png_file, width = 600, height = 600, bg = "transparent")

    d <- ncol(X_mat)
    if (d > 1) {
      pairs(cbind(X_mat[valid, , drop = FALSE], Y = Y_valid),
            col = rgb(r = red, g = 0, b = 1 - red),
            pch = 20,
            labels = c(input_names, "Y"))
    } else {
      plot(x = X_mat[valid, 1], y = Y_valid,
           xlab = input_names[1], ylab = "Y",
           col = rgb(r = red, g = 0, b = 1 - red),
           pch = 20,
           main = "PSO Optimization")
      points(x_opt[1], m, pch = 4, col = "red", cex = 2, lwd = 2)
    }

    dev.off()

    if (requireNamespace("base64enc", quietly = TRUE)) {
      img_base64 <- base64enc::base64encode(png_file)

      html_output <- sprintf(
'<div>
  <h3>PSO Optimization Results</h3>
  <p><strong>%s:</strong> %.6f</p>
  <p><strong>Found at:</strong> %s</p>
  <p><strong>Iterations:</strong> %d (swarm size: %d)</p>
  <p><strong>Total evaluations:</strong> %d</p>
  <img src="data:image/png;base64,%s" alt="PSO plot" style="max-width:600px;"/>
</div>',
        opt_type, m,
        paste(paste(names(x_opt), "=", sprintf("%.6f", x_opt), sep = ""), collapse = "; "),
        obj$state$i, obj$state$s,
        length(Y_valid),
        img_base64
      )
      analysis_dict$html <- html_output
    }

    unlink(png_file)
  }, error = function(e) {
    # If plotting fails, just skip it
  })

  return(analysis_dict)
}

# Method: get_analysis_tmp
get_analysis_tmp.PSO <- function(obj, X, Y) {
  Y_vec <- unlist(Y)
  Y_valid <- Y_vec[!is.na(Y_vec)]

  if (length(Y_valid) < 1) {
    return(list(
      text = sprintf("  Progress: Iteration %d, no valid samples yet", obj$state$i),
      data = list(iteration = obj$state$i, valid_samples = 0)
    ))
  }

  if (obj$options$yminimization) {
    current_best <- min(Y_valid)
  } else {
    current_best <- max(Y_valid)
  }

  opt_type <- if (obj$options$yminimization) "min" else "max"

  return(list(
    text = sprintf(
      "  Progress: Iteration %d, %d evaluations, current %s=%.6f",
      obj$state$i,
      length(Y_valid),
      opt_type,
      current_best
    ),
    data = list(
      iteration = obj$state$i,
      n_evaluations = length(Y_valid),
      current_best = current_best
    )
  ))
}

# Helper function: matrix_to_points
# Convert d x s matrix in [0,1] space to list of named lists in real space
matrix_to_points <- function(X_01, input_names, input_variables) {
  d <- nrow(X_01)
  s <- ncol(X_01)

  points <- list()
  for (j in 1:s) {
    point <- list()
    for (i in 1:d) {
      name <- input_names[i]
      bounds <- input_variables[[name]]
      min_val <- bounds[1]
      max_val <- bounds[2]
      point[[name]] <- X_01[i, j] * (max_val - min_val) + min_val
    }
    points[[j]] <- point
  }

  return(points)
}

# Helper function: points_to_matrix
# Convert list of named lists to n x d matrix
points_to_matrix <- function(X, input_names) {
  X_list <- lapply(X, function(point) {
    unlist(point[input_names])
  })
  X_mat <- do.call(rbind, X_list)
  if (is.null(dim(X_mat))) {
    X_mat <- matrix(X_mat, nrow = 1)
  }
  colnames(X_mat) <- input_names
  return(X_mat)
}
