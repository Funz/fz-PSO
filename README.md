# fz-PSO

Particle Swarm Optimization (PSO) algorithm for the new fz framework.

## Description

This repository contains a port of the PSO optimization algorithm from the old Funz plugin format to the new fz format. PSO is a population-based stochastic optimization technique inspired by social behavior of bird flocking or fish schooling.

## Features

- **Algorithm**: Particle Swarm Optimization (SPSO2007/SPSO2011)
- **Type**: Optimization
- **Authors**: Yann Richet <yann.richet@irsn.fr>, Claus Bendtsen <papyrus.bendtsen@gmail.com>
- **Reference**: Clerc, M. et al. (2010) http://www.particleswarm.info/standard_pso_2007.c

## Installation

Place the `PSO.R` file in your `.fz/algorithms/` directory.

## Usage

The algorithm follows the new fz S3 class structure with the following methods:

- `PSO(...)`: Constructor that accepts options as named arguments
- `get_initial_design(obj, input_variables, output_variables)`: Returns initial design points
- `get_next_design(obj, X, Y)`: Returns next batch of design points based on previous evaluations
- `get_analysis(obj, X, Y)`: Returns final analysis with visualization
- `get_analysis_tmp(obj, X, Y)`: Returns intermediate progress information

## Options

- `maxit`: Maximum number of iterations (default: 30)
- `seed`: Random seed for reproducibility (default: 123)

## Requirements

- `future`: For asynchronous optimization
- `templr`: For communication between main process and optimizer
- `base64enc`: For image encoding in HTML output

## Example

```R
# Create PSO optimizer
pso <- PSO(maxit = 50, seed = 42)

# Define input variables with bounds
input_vars <- list(
  x1 = c(0, 1),
  x2 = c(0, 1)
)

# Get initial design
initial_design <- get_initial_design(pso, input_vars, "objective")

# Evaluate objective function on initial design and get next batch
# ... (iterate until get_next_design returns empty list)

# Get final analysis
results <- get_analysis(pso, all_X, all_Y)
```

## Changes from Original

This implementation has been adapted from the original algorithm-PSO plugin to work with the new fz framework:

1. Converted from function-based to S3 class-based structure
2. Renamed functions to match new conventions (`get_initial_design`, `get_next_design`, `get_analysis`)
3. Updated to use list-based input/output format instead of matrix-based
4. Added `get_analysis_tmp` method for intermediate progress reporting
5. Enhanced HTML visualization support with base64-encoded images

## License

Same as original Funz algorithms.