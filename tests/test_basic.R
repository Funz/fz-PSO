#!/usr/bin/env Rscript

# Basic test for PSO algorithm structure
# This test verifies the algorithm can be loaded and instantiated

cat("Testing PSO algorithm structure...\n")

# Source the algorithm
source(".fz/algorithms/PSO.R")

# Test 1: Constructor
cat("Test 1: Creating PSO object...\n")
pso <- PSO(maxit = 10, seed = 42)
stopifnot(class(pso) == "PSO")
stopifnot(pso$options$maxit == 10)
stopifnot(pso$options$seed == 42)
cat("  ✓ Constructor works\n")

# Test 2: Check methods exist
cat("Test 2: Checking methods exist...\n")
stopifnot(exists("get_initial_design.PSO"))
stopifnot(exists("get_next_design.PSO"))
stopifnot(exists("get_analysis.PSO"))
stopifnot(exists("get_analysis_tmp.PSO"))
cat("  ✓ All methods defined\n")

# Test 3: Check state environment
cat("Test 3: Checking state environment...\n")
stopifnot(is.environment(pso$state))
stopifnot(exists("i", envir = pso$state))
stopifnot(exists("id", envir = pso$state))
cat("  ✓ State environment initialized\n")

cat("\nAll basic tests passed! ✓\n")
