#!/usr/bin/env python3
"""
Simple validation checks for the PSO.R file
"""

import re
import sys

def check_r_file(filepath):
    """Perform basic validation checks on R file"""
    errors = []
    warnings = []
    
    with open(filepath, 'r') as f:
        content = f.read()
        lines = content.split('\n')
    
    # Check 1: Metadata header
    if not content.startswith('#title:'):
        errors.append("Missing #title: header")
    
    # Check 2: Required metadata fields
    required_fields = ['#title:', '#author:', '#type:', '#options:', '#require:']
    for field in required_fields:
        if field not in content:
            warnings.append(f"Missing metadata field: {field}")
    
    # Check 3: S3 class definition
    if 'PSO <- function(' not in content:
        errors.append("Missing PSO constructor function")
    
    # Check 4: Required methods
    required_methods = [
        'get_initial_design.PSO',
        'get_next_design.PSO',
        'get_analysis.PSO',
        'get_analysis_tmp.PSO'
    ]
    for method in required_methods:
        if method not in content:
            errors.append(f"Missing required method: {method}")
    
    # Check 5: Balanced braces
    open_braces = content.count('{')
    close_braces = content.count('}')
    if open_braces != close_braces:
        errors.append(f"Unbalanced braces: {open_braces} open, {close_braces} close")
    
    # Check 6: Balanced parentheses (rough check)
    open_parens = content.count('(')
    close_parens = content.count(')')
    if open_parens != close_parens:
        warnings.append(f"Possibly unbalanced parentheses: {open_parens} open, {close_parens} close")
    
    # Check 7: psoptim function exists
    if 'psoptim <- function' not in content:
        errors.append("Missing psoptim function definition")
    
    # Print results
    if errors:
        print("❌ ERRORS:")
        for error in errors:
            print(f"  - {error}")
    
    if warnings:
        print("⚠️  WARNINGS:")
        for warning in warnings:
            print(f"  - {warning}")
    
    if not errors and not warnings:
        print("✅ All validation checks passed!")
        return 0
    elif not errors:
        print("\n✅ No critical errors found (only warnings)")
        return 0
    else:
        return 1

if __name__ == '__main__':
    filepath = '.fz/algorithms/PSO.R'
    sys.exit(check_r_file(filepath))
