# Porting Guide: Old PSO to New fz Format

This document details the changes made when porting from the old Funz algorithm format to the new fz format.

## Key Structural Changes

### 1. From Function-Based to S3 Class-Based

**Old Format:**
```r
PSO <- function(options) {
    algorithm = new.env()
    # ... initialization
    return(algorithm)
}
```

**New Format:**
```r
PSO <- function(...) {
    opts <- list(...)
    state <- new.env(parent = emptyenv())
    obj <- list(options = ..., state = state)
    class(obj) <- "PSO"
    return(obj)
}
```

### 2. Function Naming Convention

| Old Format | New Format |
|------------|------------|
| `PSO()` | `PSO()` (constructor remains same) |
| `getInitialDesign()` | `get_initial_design.PSO()` (S3 method) |
| `getNextDesign()` | `get_next_design.PSO()` (S3 method) |
| `displayResults()` | `get_analysis.PSO()` (S3 method) |
| `displayResultsTmp()` | `get_analysis_tmp.PSO()` (S3 method) |

### 3. Input/Output Format

**Old Format (Matrix-based):**
```r
# Input: matrix with column names
Xn <- matrix(...)
colnames(Xn) <- names(algorithm$input)

# Output: matrix or vector
Y <- matrix(...)
```

**New Format (List-based):**
```r
# Input: list of named lists
samples <- list()
for (i in 1:n) {
    sample <- list(x1 = value1, x2 = value2, ...)
    samples[[i]] <- sample
}

# Output: list of values
Y <- list(value1, value2, ...)
```

### 4. Variable Definition

**Old Format:**
```r
input <- list(
    x1 = list(min = 0, max = 1),
    x2 = list(min = 0, max = 1)
)
```

**New Format:**
```r
input_variables <- list(
    x1 = c(0, 1),  # c(min, max)
    x2 = c(0, 1)
)
```

### 5. Return Values for Analysis

**Old Format (HTML string):**
```r
displayResults <- function(...) {
    return(paste(sep = "",
                 "<HTML name='points'>",
                 "<img src='", file, "'/>",
                 "</HTML>",
                 "<min>", min(Y), "</min>"))
}
```

**New Format (Structured list):**
```r
get_analysis.PSO <- function(...) {
    analysis_dict <- list(
        text = "Human-readable summary\n...",
        data = list(min = ..., argmin = ...),
        html = "<div>...</div>"  # Optional
    )
    return(analysis_dict)
}
```

### 6. Termination Signal

**Old Format:**
```r
getNextDesign <- function(...) {
    if (finished) return(NULL)
    # ...
}
```

**New Format:**
```r
get_next_design.PSO <- function(...) {
    if (finished) return(list())  # Empty list
    # ...
}
```

## Implementation Details

### Metadata Header

Both formats support metadata in comments at the top of the file:

```r
#title: Particle Swarm Optimization
#author: Author Name <email>
#type: optimization
#options: maxit=30;seed=123
#require: future;templr;base64enc
#references: Reference URL or citation
```

### State Management

**Old Format:**
- Used an environment directly as the algorithm object
- All state stored in this environment

**New Format:**
- Separates options (immutable) from state (mutable)
- State stored in nested environment for clarity
- More idiomatic S3 pattern

### Async Execution

Both formats use the same `future` and `templr` packages for async PSO execution:
- The core `psoptim` function remains unchanged
- Communication via `ask_X`, `tell_Y`, `ask_Y` remains the same
- Only the wrapping and variable capturing changed slightly

## Files Created

```
.fz/
└── algorithms/
    └── PSO.R               # Main algorithm (ported)

tests/
├── test_basic.R            # Basic structure tests
└── validate_syntax.py      # Syntax validation

examples/
└── branin_example.R        # Example usage

.github/
└── workflows/
    └── test.yml            # CI/CD workflow

README.md                    # Documentation
.gitignore                  # Ignore patterns
```

## Testing Strategy

1. **Structure Tests** (`test_basic.R`): Verify S3 class structure and methods exist
2. **Syntax Validation** (`validate_syntax.py`): Check file structure and balanced braces
3. **Example** (`branin_example.R`): Functional test with Branin optimization
4. **CI/CD** (`.github/workflows/test.yml`): Automated testing on GitHub

## Migration Checklist for Other Algorithms

- [ ] Convert main function to S3 constructor
- [ ] Rename functions to S3 method convention
- [ ] Update input format from matrix to list of named lists
- [ ] Update output format from vectors to lists
- [ ] Update return format for analysis methods
- [ ] Change termination signal from NULL to empty list
- [ ] Add generic function definitions
- [ ] Update metadata header
- [ ] Create tests
- [ ] Add examples
- [ ] Document changes
