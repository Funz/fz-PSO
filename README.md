# fz-PSO

A [Funz](https://github.com/Funz/fz) algorithm plugin implementing **Particle Swarm Optimization (PSO)** in R.

This repository provides the PSO algorithm for the fz framework, ported from the original [algorithm-PSO](https://github.com/Funz/algorithm-PSO) Funz plugin.

**Algorithm reference:** [Clerc, M. et al. (2010) *Particle Swarm Optimization*](http://clerc.maurice.free.fr/pso/)

## Features

### Algorithm Interface (R S3 class)

The algorithm implements the fz R algorithm interface:

- `PSO(...)`: S3 constructor accepting algorithm-specific options
- `get_initial_design.PSO(obj, input_variables, output_variables)`: Return initial swarm of particles
- `get_next_design.PSO(obj, X, Y)`: Update velocities and positions, return next swarm, or `list()` when max iterations reached
- `get_analysis.PSO(obj, X, Y)`: Return optimum value, location, and optional visualization
- `get_analysis_tmp.PSO(obj, X, Y)`: Return intermediate progress (current iteration, best value)

### Algorithm Behavior

1. **Initialization**: Creates a swarm of `nparticles` particles randomly distributed in the input space, with random initial velocities.

2. **Iteration**: At each iteration:
   - Updates personal best positions (each particle remembers its best location)
   - Updates global best (the best position found by any particle)
   - Updates velocities using inertia, cognitive (personal best), and social (global best) components
   - Updates positions and clamps to bounds with velocity reset

3. **Convergence**: Stops when:
   - Maximum iterations reached

## Requirements

- **R** must be installed on your system
- **rpy2** Python package: `pip install rpy2`
- **fz** framework: `pip install git+https://github.com/Funz/fz.git`

## Installation

```bash
pip install git+https://github.com/Funz/fz.git
pip install rpy2
```

### Install the Algorithm Plugin

```python
import fz
fz.install_algorithm("PSO")
```

Or from a URL:
```python
fz.install_algorithm("https://github.com/Funz/fz-PSO")
```

Or using the CLI:
```bash
fz install PSO
```

## Algorithm Options

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `yminimization` | boolean | `true` | Minimize output value? Set to `false` for maximization |
| `max_iterations` | integer | `30` | Maximum number of iterations |
| `nparticles` | integer | `NA` (auto) | Number of particles in swarm (default: `10 + 2*sqrt(d)` where `d` = number of dimensions) |
| `seed` | integer | `123` | Random seed for reproducibility |
| `w` | numeric | `0.7213` | Inertia weight (`1/(2*log(2))`) |
| `c_p` | numeric | `1.1931` | Cognitive (personal best) coefficient (`0.5+log(2)`) |
| `c_g` | numeric | `1.1931` | Social (global best) coefficient (`0.5+log(2)`) |

## Usage

### Without fzd (standalone algorithm testing)

You can test the algorithm without any simulation code, using rpy2 directly:

```python
from rpy2 import robjects

# Source the R algorithm
robjects.r.source(".fz/algorithms/PSO.R")
r_globals = robjects.globalenv

# Create an instance
r_algo = robjects.r["PSO"](
    yminimization=True, max_iterations=20,
    nparticles=15, seed=123
)

# Define input variable ranges as R list
r_input_vars = robjects.r('list(x1 = c(0.0, 1.0), x2 = c(0.0, 1.0))')
r_output_vars = robjects.StrVector(["output"])

# Get initial design (swarm of particles)
r_design = r_globals['get_initial_design'](r_algo, r_input_vars, r_output_vars)
print(f"Initial design: {len(r_design)} particles")
```

Or via fz's automatic wrapper:

```python
from fz.algorithms import load_algorithm

# Load R algorithm (fz handles rpy2 wrapping automatically)
algo = load_algorithm("PSO",
                      yminimization=True, max_iterations=20, nparticles=15)

# Same Python interface as Python algorithms
design = algo.get_initial_design(
    {"x1": (0.0, 1.0), "x2": (0.0, 1.0)}, ["output"]
)
print(f"Initial design: {len(design)} particles")
```

### With fzd (coupled with a model)

Use `fz.fzd()` to run the algorithm coupled with a model and calculators:

```python
import fz

# Install model and algorithm plugins
fz.install("Model")
fz.install_algorithm("PSO")

# Run optimization
analysis = fz.fzd(
    input_path="examples/Model/input.txt",
    input_variables={"x": "[0;10]", "y": "[-5;5]"},
    model="Model",
    output_expression="result",
    algorithm="PSO",
    algorithm_options={
        "yminimization": True,
        "max_iterations": 30,
        "nparticles": 20,
        "seed": 123
    },
    calculators="localhost_Model",
    analysis_dir="analysis_results_pso"
)

print(analysis)
```

## Output

The algorithm provides:

- **Final analysis**:
  - Optimum value found
  - Location of optimum (input values)
  - Number of iterations and evaluations
  - Swarm size used
  - Visualization plot (pairs plot for multi-dimensional, scatter plot for 1D)

- **Intermediate progress**:
  - Current iteration number
  - Current best value
  - Number of evaluations so far

## Directory Structure

```
fz-PSO/
├── .fz/
│   └── algorithms/
│       └── PSO.R                     # R algorithm implementation (S3 class)
├── .github/
│   └── workflows/
│       └── test.yml                  # CI workflow (includes R setup)
├── tests/
│   └── test_plugin.py                # Test suite (uses rpy2)
├── example_standalone.ipynb          # Notebook: algorithm without fzd
├── example_with_fzd.ipynb            # Notebook: algorithm with fzd
├── LICENSE
└── README.md
```

## Technical Details

### Algorithm Type
- **Type**: Global optimization
- **Order**: Zero-order (derivative-free)
- **Method**: Particle Swarm Optimization (SPSO 2007-style velocity update)

### Dependencies
- Base R (no special packages required)
- Optional: `base64enc` for HTML visualization output

### Porting Notes

This algorithm has been ported from the original Funz plugin format to the new fz format:
- Original: [algorithm-PSO](https://github.com/Funz/algorithm-PSO)

Key changes from old format:
- Constructor pattern using S3 classes instead of environments
- Methods renamed: `getInitialDesign` → `get_initial_design`, `getNextDesign` → `get_next_design`, etc.
- Removed `future`/`templr` async dependency — now uses synchronous step-by-step interface
- Input/output format adapted to new fz expectations (list of points instead of matrices)
- Return empty `list()` to signal completion instead of `NULL`
- State management using environments for mutable state
- Simplified PSO core: inertia + cognitive + social velocity update with boundary clamping

## Running Tests

```bash
# Run all tests
python -m pytest tests/ -v

# Or directly
python tests/test_plugin.py
```

## License

BSD 3-Clause License (same as original Funz project)

## Authors

Yann Richet, Claus Bendtsen (ported to new fz format)