# Example: Using PSO for Branin function optimization

# This example demonstrates how to use the PSO algorithm
# to optimize the Branin function, a common benchmark

# Load the algorithm
source(".fz/algorithms/PSO.R")

# Define the Branin test function
branin <- function(x) {
  x1 <- x$x1 * 15 - 5   
  x2 <- x$x2 * 15     
  (x2 - 5/(4*pi^2)*(x1^2) + 5/pi*x1 - 6)^2 + 10*(1 - 1/(8*pi))*cos(x1) + 10
}

# Create PSO optimizer
pso <- PSO(maxit = 30, seed = 123)

# Define input variables (normalized to [0,1])
input_vars <- list(
  x1 = c(0, 1),
  x2 = c(0, 1)
)

output_vars <- "branin"

cat("Starting PSO optimization of Branin function\n")
cat("Expected minimum: 0.3978874\n\n")

# Get initial design
cat("Getting initial design...\n")
X_initial <- get_initial_design(pso, input_vars, output_vars)

# In a real scenario, you would evaluate these points
# and collect results. For this example, we'll simulate:

# Evaluate initial design
Y_initial <- lapply(X_initial, branin)

X_all <- X_initial
Y_all <- Y_initial

# Iterate
iteration <- 1
repeat {
  cat(sprintf("Iteration %d: Current best = %.6f\n", 
              iteration, 
              min(unlist(Y_all))))
  
  # Get next design
  X_next <- get_next_design(pso, X_all, Y_all)
  
  # Check if finished
  if (length(X_next) == 0) {
    cat("Optimization complete!\n")
    break
  }
  
  # Evaluate new points
  Y_next <- lapply(X_next, branin)
  
  # Append to history
  X_all <- c(X_all, X_next)
  Y_all <- c(Y_all, Y_next)
  
  iteration <- iteration + 1
}

# Get final analysis
cat("\nFinal Analysis:\n")
results <- get_analysis(pso, X_all, Y_all)
cat(results$text)
cat("\nExpected minimum: 0.3978874\n")
cat(sprintf("Found minimum: %.6f\n", results$data$min))
cat(sprintf("Error: %.6f\n", abs(results$data$min - 0.3978874)))
