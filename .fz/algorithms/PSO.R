#title: Particle Swarm Optimization
#author: Yann Richet <yann.richet@irsn.fr>, Claus Bendtsen <papyrus.bendtsen@gmail.com>
#type: optimization
#options: maxit=30;seed=123
#require: future;templr;base64enc
#references: Clerc, M. et al. (2010) http://www.particleswarm.info/standard_pso_2007.c

# Constructor for PSO S3 class
PSO <- function(...) {
  # Get options from ... arguments
  opts <- list(...)
  
  # Create object with initial state
  # Use an environment for mutable state (idiomatic S3 pattern)
  state <- new.env(parent = emptyenv())
  state$i <- 0
  state$s <- NULL
  state$id <- floor(1000 * runif(1))
  state$input_variables <- list()
  state$output_variables <- list()
  state$d <- 0
  state$job <- NULL
  state$files <- NULL
  
  obj <- list(
    options = list(
      maxit = as.integer(
        ifelse(is.null(opts$maxit), 30, opts$maxit)
      ),
      seed = as.integer(
        ifelse(is.null(opts$seed), 123, opts$seed)
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

# Helper functions to extract min/max from input variables
min_input <- function(input_vars) {
  sapply(input_vars, function(v) v[1])
}

max_input <- function(input_vars) {
  sapply(input_vars, function(v) v[2])
}

# Method: get_initial_design
get_initial_design.PSO <- function(obj, input_variables, output_variables) {
  set.seed(obj$options$seed)
  
  # Store variable information in mutable state
  # input_variables is a named list: list(var1 = c(min, max), var2 = c(min, max))
  obj$state$input_variables <- input_variables
  obj$state$output_variables <- output_variables
  obj$state$d <- length(input_variables)
  
  library(future)
  library(templr)
  wd <- getwd()
  
  # Capture variables for future
  id <- obj$state$id
  d <- obj$state$d
  seed <- obj$options$seed
  maxit <- obj$options$maxit
  
  # Start async PSO optimization job
  obj$state$job <- future(evaluator = plan("multisession"), lazy = FALSE, {
    sink(file.path(wd, paste0('PSO_', id, '.out')), type = 'output')
    print("Starting psoptim()")
    set.seed(seed)
    o <- psoptim(
      par = (min_input(input_variables) + max_input(input_variables)) / 2,
      fn = function(x) {
        ask_Y(id = id, x = matrix(x, ncol = d))
      },
      lower = min_input(input_variables),
      upper = max_input(input_variables),
      control = list(maxit = maxit, vectorized = TRUE)
    )
    print("psoptim() ended")
    print(o)
    sink(type = 'output')
    ask_Y(id = id, matrix(NaN, ncol = d))
  })
  
  obj$state$i <- 0
  
  Sys.sleep(.1)
  
  Xn <- ask_X(id = obj$state$id)
  obj$state$s <- nrow(Xn)
  
  # Convert matrix to list of named lists
  samples <- list()
  for (i in 1:nrow(Xn)) {
    sample <- list()
    for (j in 1:length(input_variables)) {
      sample[[names(input_variables)[j]]] <- Xn[i, j]
    }
    samples[[i]] <- sample
  }
  
  return(samples)
}

# Method: get_next_design
get_next_design.PSO <- function(obj, X, Y) {
  obj$state$i <- obj$state$i + 1
  
  # Extract last s results
  n <- length(Y)
  y <- Y[(n - obj$state$s + 1):n]
  y <- unlist(y)
  
  if (all(is.na(y))) {
    tell_Y(id = obj$state$id, NULL)
    return(list())
  } else {
    tell_Y(id = obj$state$id, y)
  }
  
  Sys.sleep(.1)
  
  Xn <- ask_X(id = obj$state$id)
  if (is.null(Xn)) return(list())
  obj$state$s <- nrow(Xn)
  
  # Convert matrix to list of named lists
  samples <- list()
  for (i in 1:nrow(Xn)) {
    sample <- list()
    for (j in 1:length(obj$state$input_variables)) {
      sample[[names(obj$state$input_variables)[j]]] <- Xn[i, j]
    }
    samples[[i]] <- sample
  }
  
  return(samples)
}

# Method: get_analysis
get_analysis.PSO <- function(obj, X, Y) {
  analysis_dict <- list(text = "", data = list())
  
  # Filter out NULL/NA values
  Y_valid <- Y[!sapply(Y, is.null) & !is.na(Y)]
  Y_valid <- unlist(Y_valid)
  
  if (length(Y_valid) == 0) {
    analysis_dict$text <- "No valid results to analyze"
    return(analysis_dict)
  }
  
  # Convert X to matrix
  X_matrix <- matrix(unlist(X), ncol = obj$state$d, byrow = TRUE)
  colnames(X_matrix) <- names(obj$state$input_variables)
  
  # Find minimum
  min_idx <- which.min(Y_valid)
  min_y <- Y_valid[min_idx]
  argmin_x <- X_matrix[min_idx, ]
  
  # Store data
  analysis_dict$data <- list(
    min = min_y,
    argmin = as.list(argmin_x),
    n_evaluations = length(Y_valid)
  )
  
  # Create text summary
  argmin_str <- paste(names(argmin_x), "=", round(argmin_x, 6), collapse = ", ")
  analysis_dict$text <- sprintf(
"PSO Optimization Results:
  Total evaluations: %d
  Minimum: %.6f
  At: %s
",
    length(Y_valid),
    min_y,
    argmin_str
  )
  
  # Try to create visualization
  tryCatch({
    obj$state$files <- paste0("pso_plot_", obj$state$i, ".png")
    png(file = obj$state$files, height = 600, width = 600)
    
    if (ncol(X_matrix) > 1) {
      red <- (Y_valid - min(Y_valid)) / (max(Y_valid) - min(Y_valid))
      pairs(X_matrix, col = rgb(r = red, g = 0, b = 1 - red))
    } else {
      plot(X_matrix, Y_valid, xlab = names(obj$state$input_variables)[1], 
           ylab = "Objective")
    }
    
    dev.off()
    
    # Convert to base64 for HTML
    if (requireNamespace("base64enc", quietly = TRUE)) {
      img_base64 <- base64enc::base64encode(obj$state$files)
      
      html_output <- sprintf(
'<div>
  <p><strong>Minimum:</strong> %.6f</p>
  <p><strong>At:</strong> %s</p>
  <img src="data:image/png;base64,%s" alt="PSO Results" style="max-width:600px;"/>
</div>',
        min_y,
        argmin_str,
        img_base64
      )
      analysis_dict$html <- html_output
    }
  }, error = function(e) {
    # If plotting fails, just skip it
  })
  
  return(analysis_dict)
}

# Method: get_analysis_tmp
get_analysis_tmp.PSO <- function(obj, X, Y) {
  # Filter out NULL/NA values
  Y_valid <- Y[!sapply(Y, is.null) & !is.na(Y)]
  Y_valid <- unlist(Y_valid)
  
  if (length(Y_valid) == 0) {
    return(list(
      text = sprintf("  Progress: iteration %d, 0 evaluations", obj$state$i),
      data = list(iteration = obj$state$i, n_evaluations = 0)
    ))
  }
  
  min_y <- min(Y_valid)
  
  return(list(
    text = sprintf(
      "  Progress: iteration %d, %d evaluations, current min=%.6f",
      obj$state$i,
      length(Y_valid),
      min_y
    ),
    data = list(
      iteration = obj$state$i,
      n_evaluations = length(Y_valid),
      current_min = min_y
    )
  ))
}

#' @ref https://cran.r-project.org/web/packages/pso/ , with slight mods for vectorized objective fun
psoptim <- function (par, fn, gr = NULL, ..., lower = -1, upper = 1,
                     control = list()) {
    
    fn1 <- function(par) fn(par, ...)/p.fnscale
    mrunif <- function(n, m, lower, upper) {
        return(matrix(runif(n * m, 0, 1), nrow = n, ncol = m) * (upper - lower) + lower)
    }
    norm <- function(x) sqrt(sum(x * x))
    rsphere.unif <- function(n, r) {
        temp <- runif(n)
        return((runif(1, min = 0, max = r) / norm(temp)) * temp)
    }
    svect <- function(a, b, n, k) {
        temp <- rep(a, n)
        temp[k] <- b
        return(temp)
    }
    mrsphere.unif <- function(n, r) {
        m <- length(r)
        temp <- matrix(runif(n * m), n, m)
        return(temp %*% diag(runif(m, min = 0, max = r) / apply(temp, 2, norm)))
    }
    npar <- length(par)
    lower <- as.double(rep(lower, , npar))
    upper <- as.double(rep(upper, , npar))
    con <- list(trace = 0, fnscale = 1, maxit = 1000L, maxf = Inf,
                abstol = -Inf, reltol = 0, REPORT = 10,
                s = NA, k = 3, p = NA, w = 1/(2*log(2)),
                c.p = .5+log(2), c.g = .5+log(2), d = NA,
                v.max = NA, rand.order = TRUE, max.restart = Inf,
                maxit.stagnate = Inf,
                vectorized = FALSE, hybrid = FALSE, hybrid.control = NULL,
                trace.stats = FALSE, type = "SPSO2007")
    nmsC <- names(con)
    con[(namc <- names(control))] <- control
    if (length(noNms <- namc[!namc %in% nmsC])) 
        warning("unknown names in control: ", paste(noNms, collapse = ", "))
    ## Argument error checks
    if (any(upper == Inf | lower == -Inf))
        stop("fixed bounds must be provided")
    
    p.type <- pmatch(con[["type"]], c("SPSO2007", "SPSO2011")) - 1
    if (is.na(p.type)) stop("type should be one of \"SPSO2007\", \"SPSO2011\"")
    
    p.trace <- con[["trace"]] > 0L # provide output on progress?
    p.fnscale <- con[["fnscale"]] # scale funcion by 1/fnscale
    p.maxit <- con[["maxit"]] # maximal number of iterations
    p.maxf <- con[["maxf"]] # maximal number of function evaluations
    p.abstol <- con[["abstol"]] # absolute tolerance for convergence
    p.reltol <- con[["reltol"]] # relative minimal tolerance for restarting
    p.report <- as.integer(con[["REPORT"]]) # output every REPORT iterations
    p.s <- ifelse(is.na(con[["s"]]), ifelse(p.type == 0, floor(10 + 2 * sqrt(npar)), 40),
                  con[["s"]]) # swarm size
    p.p <- ifelse(is.na(con[["p"]]), 1 - (1 - 1 / p.s)^con[["k"]], con[["p"]]) # average % of informants
    p.w0 <- con[["w"]] # exploitation constant
    if (length(p.w0) > 1) {
        p.w1 <- p.w0[2]
        p.w0 <- p.w0[1]
    } else {
        p.w1 <- p.w0
    }
    p.c.p <- con[["c.p"]] # local exploration constant
    p.c.g <- con[["c.g"]] # global exploration constant
    p.d <- ifelse(is.na(con[["d"]]), norm(upper - lower), con[["d"]]) # domain diameter
    p.vmax <- con[["v.max"]] * p.d # maximal velocity
    p.randorder <- as.logical(con[["rand.order"]]) # process particles in random order?
    p.maxrestart <- con[["max.restart"]] # maximal number of restarts
    p.maxstagnate <- con[["maxit.stagnate"]] # maximal number of iterations without improvement
    p.vectorized <- as.logical(con[["vectorized"]]) # vectorized?
    if (is.character(con[["hybrid"]])) {
        p.hybrid <- pmatch(con[["hybrid"]], c("off", "on", "improved")) - 1
        if (is.na(p.hybrid)) stop("hybrid should be one of \"off\", \"on\", \"improved\"")
    } else {
        p.hybrid <- as.integer(as.logical(con[["hybrid"]])) # use local BFGS search
    }
    p.hcontrol <- con[["hybrid.control"]] # control parameters for hybrid optim
    if ("fnscale" %in% names(p.hcontrol))
        p.hcontrol["fnscale"] <- p.hcontrol["fnscale"] * p.fnscale
    else
        p.hcontrol["fnscale"] <- p.fnscale
    p.trace.stats <- as.logical(con[["trace.stats"]]) # collect detailed stats?
    
    if (p.trace) {
        message("S=", p.s, ", K=", con[["k"]], ", p=", signif(p.p, 4), ", w0=",
                signif(p.w0, 4), ", w1=",
                signif(p.w1, 4), ", c.p=", signif(p.c.p, 4),
                ", c.g=", signif(p.c.g, 4))
        message("v.max=", signif(con[["v.max"]], 4),
                ", d=", signif(p.d, 4), ", vectorized=", p.vectorized,
                ", hybrid=", c("off", "on", "improved")[p.hybrid + 1])
        if (p.trace.stats) {
            stats.trace.it <- c()
            stats.trace.error <- c()
            stats.trace.f <- NULL
            stats.trace.x <- NULL
        }
    }
    ## Initialization
    if (p.reltol != 0) p.reltol <- p.reltol * p.d
    if (p.vectorized) {
        lowerM <- matrix(lower, nrow = npar, ncol = p.s)
        upperM <- matrix(upper, nrow = npar, ncol = p.s)
    }
    X <- mrunif(npar, p.s, lower, upper)
    if (!any(is.na(par)) && all(par >= lower) && all(par <= upper)) X[, 1] <- par
    if (p.type == 0) {
        V <- (mrunif(npar, p.s, lower, upper) - X) / 2
    } else { ## p.type==1
        V <- matrix(runif(npar * p.s, min = as.vector(lower - X), max = as.vector(upper - X)), npar, p.s)
        p.c.p2 <- p.c.p / 2 # precompute constants
        p.c.p3 <- p.c.p / 3
        p.c.g3 <- p.c.g / 3
        p.c.pg3 <- p.c.p3 + p.c.g3
    }
    if (!is.na(p.vmax)) { # scale to maximal velocity
        temp <- apply(V, 2, norm)
        temp <- pmin.int(temp, p.vmax) / temp
        V <- V %*% diag(temp)
    }
    f.x <- fn1(t(X))
    stats.feval <- p.s
    P <- X
    f.p <- f.x
    P.improved <- rep(FALSE, p.s)
    i.best <- which.min(f.p)
    error <- f.p[i.best]
    init.links <- TRUE
    if (p.trace && p.report == 1) {
        message("It 1: fitness=", signif(error, 4))
        if (p.trace.stats) {
            stats.trace.it <- c(stats.trace.it, 1)
            stats.trace.error <- c(stats.trace.error, error)
            stats.trace.f <- c(stats.trace.f, list(f.x))
            stats.trace.x <- c(stats.trace.x, list(X))
        }
    }
    ## Iterations
    stats.iter <- 1
    stats.restart <- 0
    stats.stagnate <- 0
    while (stats.iter < p.maxit && stats.feval < p.maxf && error > p.abstol &&
           stats.restart < p.maxrestart && stats.stagnate < p.maxstagnate) {
        stats.iter <- stats.iter + 1
        if (p.p != 1 && init.links) {
            links <- matrix(runif(p.s * p.s, 0, 1) <= p.p, p.s, p.s)
            diag(links) <- TRUE
        }
        ## The swarm moves
        if (!p.vectorized) {
            if (p.randorder) {
                index <- sample(p.s)
            } else {
                index <- 1:p.s
            }
            for (i in index) {
                if (p.p == 1)
                    j <- i.best
                else
                    j <- which(links[, i])[which.min(f.p[links[, i]])] # best informant
                temp <- (p.w0 + (p.w1 - p.w0) * max(stats.iter / p.maxit, stats.feval / p.maxf))
                V[, i] <- temp * V[, i] # exploration tendency
                if (p.type == 0) {
                    V[, i] <- V[, i] + runif(npar, 0, p.c.p) * (P[, i] - X[, i]) # exploitation
                    if (i != j) V[, i] <- V[, i] + runif(npar, 0, p.c.g) * (P[, j] - X[, i])
                } else { # SPSO 2011
                    if (i != j)
                        temp <- p.c.p3 * P[, i] + p.c.g3 * P[, j] - p.c.pg3 * X[, i] # Gi-Xi
                    else
                        temp <- p.c.p2 * P[, i] - p.c.p2 * X[, i] # Gi-Xi for local=best
                    V[, i] <- V[, i] + temp + rsphere.unif(npar, norm(temp))
                }
                if (!is.na(p.vmax)) {
                    temp <- norm(V[, i])
                    if (temp > p.vmax) V[, i] <- (p.vmax / temp) * V[, i]
                }
                X[, i] <- X[, i] + V[, i]
                ## Check bounds
                temp <- X[, i] < lower
                if (any(temp)) {
                    X[temp, i] <- lower[temp]
                    V[temp, i] <- 0
                }
                temp <- X[, i] > upper
                if (any(temp)) {
                    X[temp, i] <- upper[temp]
                    V[temp, i] <- 0
                }
                ## Evaluate function
                if (p.hybrid == 1) {
                    temp <- optim(X[, i], fn, gr, ..., method = "L-BFGS-B", lower = lower,
                                  upper = upper, control = p.hcontrol)
                    V[, i] <- V[, i] + temp$par - X[, i] # disregards any v.max imposed
                    X[, i] <- temp$par
                    f.x[i] <- temp$value
                    stats.feval <- stats.feval + as.integer(temp$counts[1])
                } else {
                    f.x[i] <- fn1(X[, i])
                    stats.feval <- stats.feval + 1
                }
                if (f.x[i] < f.p[i]) { # improvement
                    P[, i] <- X[, i]
                    f.p[i] <- f.x[i]
                    if (f.p[i] < f.p[i.best]) {
                        i.best <- i
                        if (p.hybrid == 2) {
                            temp <- optim(X[, i], fn, gr, ..., method = "L-BFGS-B", lower = lower,
                                          upper = upper, control = p.hcontrol)
                            V[, i] <- V[, i] + temp$par - X[, i] # disregards any v.max imposed
                            X[, i] <- temp$par
                            P[, i] <- temp$par
                            f.x[i] <- temp$value
                            f.p[i] <- temp$value
                            stats.feval <- stats.feval + as.integer(temp$counts[1])
                        }
                    }
                }
                if (stats.feval >= p.maxf) break
            }
        } else {
            if (p.p == 1)
                j <- rep(i.best, p.s)
            else # best informant
                j <- sapply(1:p.s, function(i)
                    which(links[, i])[which.min(f.p[links[, i]])]) 
            temp <- (p.w0 + (p.w1 - p.w0) * max(stats.iter / p.maxit, stats.feval / p.maxf))
            V <- temp * V # exploration tendency
            if (p.type == 0) {
                V <- V + mrunif(npar, p.s, 0, p.c.p) * (P - X) # exploitation
                temp <- j != (1:p.s)
                V[, temp] <- V[, temp] + mrunif(npar, sum(temp), 0, p.c.p) * (P[, j[temp]] - X[, temp])
            } else { # SPSO 2011
                temp <- j == (1:p.s)
                temp <- P %*% diag(svect(p.c.p3, p.c.p2, p.s, temp)) +
                    P[, j] %*% diag(svect(p.c.g3, 0, p.s, temp)) -
                    X %*% diag(svect(p.c.pg3, p.c.p2, p.s, temp)) # G-X
                V <- V + temp + mrsphere.unif(npar, apply(temp, 2, norm))
            }
            if (!is.na(p.vmax)) {
                temp <- apply(V, 2, norm)
                temp <- pmin.int(temp, p.vmax) / temp
                V <- V %*% diag(temp)
            }
            X <- X + V
            ## Check bounds
            temp <- X < lowerM
            if (any(temp)) {
                X[temp] <- lowerM[temp] 
                V[temp] <- 0
            }
            temp <- X > upperM
            if (any(temp)) {
                X[temp] <- upperM[temp]
                V[temp] <- 0
            }
            ## Evaluate function
            if (p.hybrid == 1) { # not really vectorizing
                for (i in 1:p.s) {
                    temp <- optim(X[, i], fn, gr, ..., method = "L-BFGS-B", lower = lower,
                                  upper = upper, control = p.hcontrol)
                    V[, i] <- V[, i] + temp$par - X[, i] # disregards any v.max imposed
                    X[, i] <- temp$par
                    f.x[i] <- temp$value
                    stats.feval <- stats.feval + as.integer(temp$counts[1])
                }
            } else {
                f.x <- fn1(t(X))
                stats.feval <- stats.feval + p.s
            }
            temp <- sapply(isTRUE, X = as.numeric(f.x) < f.p)
            if (any(temp)) { # improvement
                P[, temp] <- X[, temp]
                f.p[temp] <- f.x[temp]
                i.best <- which.min(f.p)
                if (temp[i.best] && p.hybrid == 2) { # overall improvement
                    temp <- optim(X[, i.best], fn, gr, ..., method = "L-BFGS-B", lower = lower,
                                  upper = upper, control = p.hcontrol)
                    V[, i.best] <- V[, i.best] + temp$par - X[, i.best] # disregards any v.max imposed
                    X[, i.best] <- temp$par
                    P[, i.best] <- temp$par
                    f.x[i.best] <- temp$value
                    f.p[i.best] <- temp$value
                    stats.feval <- stats.feval + as.integer(temp$counts[1])
                }
            }
            if (stats.feval >= p.maxf) break
        }
        if (p.reltol != 0) {
            d <- X - P[, i.best]
            d <- sqrt(max(colSums(d * d)))
            if (d < p.reltol) {
                X <- mrunif(npar, p.s, lower, upper)
                V <- (mrunif(npar, p.s, lower, upper) - X) / 2
                if (!is.na(p.vmax)) {
                    temp <- apply(V, 2, norm)
                    temp <- pmin.int(temp, p.vmax) / temp
                    V <- V %*% diag(temp)
                }
                stats.restart <- stats.restart + 1
                if (p.trace) message("It ", stats.iter, ": restarting")
            }
        }
        init.links <- f.p[i.best] == error # if no overall improvement
        stats.stagnate <- ifelse(init.links, stats.stagnate + 1, 0)
        error <- f.p[i.best]
        if (p.trace && stats.iter %% p.report == 0) {
            if (p.reltol != 0) 
                message("It ", stats.iter, ": fitness=", signif(error, 4),
                        ", swarm diam.=", signif(d, 4))
            else
                message("It ", stats.iter, ": fitness=", signif(error, 4))
            if (p.trace.stats) {
                stats.trace.it <- c(stats.trace.it, stats.iter)
                stats.trace.error <- c(stats.trace.error, error)
                stats.trace.f <- c(stats.trace.f, list(f.x))
                stats.trace.x <- c(stats.trace.x, list(X))
            }
        }
    }
    if (error <= p.abstol) {
        msg <- "Converged"
        msgcode <- 0
    } else if (stats.feval >= p.maxf) {
        msg <- "Maximal number of function evaluations reached"
        msgcode <- 1
    } else if (stats.iter >= p.maxit) {
        msg <- "Maximal number of iterations reached"
        msgcode <- 2
    } else if (stats.restart >= p.maxrestart) {
        msg <- "Maximal number of restarts reached"
        msgcode <- 3
    } else {
        msg <- "Maximal number of iterations without improvement reached"
        msgcode <- 4
    }
    if (p.trace) message(msg)
    o <- list(par = P[, i.best], value = f.p[i.best],
              counts = c("function" = stats.feval, "iteration" = stats.iter,
                       "restarts" = stats.restart),
              convergence = msgcode, message = msg)
    if (p.trace && p.trace.stats) o <- c(o, list(stats = list(it = stats.trace.it,
                                                           error = stats.trace.error,
                                                           f = stats.trace.f,
                                                           x = stats.trace.x)))
    
    return(o)
}
