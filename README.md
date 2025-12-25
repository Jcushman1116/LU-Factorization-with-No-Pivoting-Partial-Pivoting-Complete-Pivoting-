# LU Factorization with No Pivoting, Partial Pivoting, and Complete Pivoting

## Overview
This project implements a **dense LU factorization algorithm** using **Gaussian elimination** under three pivoting strategies:

- **No pivoting:** `A = LU`
- **Partial pivoting:** `PA = LU`
- **Complete pivoting:** `PAQ = LU`

The algorithm is tested empirically on a wide range of **structured matrices** to evaluate **numerical stability**, **accuracy**, and **growth behavior** as matrix size increases. Performance is assessed using **factorization error**, **residual error**, **growth factor**, and **condition number** metrics.

All routines are implemented in **MATLAB**.

---

## Methods Implemented

### 1. LU Factorization Without Pivoting
Assumes the matrix is well-conditioned (e.g., diagonally dominant).

- Fastest method
- No row or column swaps
- Fails when pivot elements are zero or very small
- Least numerically stable

---

### 2. LU Factorization with Partial Pivoting
Computes

    P A = L U

by swapping rows to move the largest element in the active column into the pivot position.

- Improves numerical stability
- Adds an `O(n²)` pivot search cost
- Standard method used in practice

---

### 3. LU Factorization with Complete Pivoting
Computes

    P A Q = L U

by swapping both rows and columns to select the largest element in the active submatrix.

- Most numerically stable
- Highest computational cost
- Best control of growth factor

---

## Implementation Details

- All methods are implemented in **one function** using a routine-selection flag
- `L` and `U` are stored **in-place** within a single matrix for memory efficiency
- Row and column permutations are tracked using **permutation vectors**
- Rank-one updates are used to efficiently update trailing submatrices
- Error guards prevent division by near-zero pivots (`< 1e-12`)

### Computational Complexity
- No pivoting: `O(n³)`
- Partial pivoting: `O(n³)`
- Complete pivoting: `O(n³)` with higher constants due to pivot search

---

## Experimental Design

Matrices of increasing size were tested under various **structured configurations**, with multiple trials per size. The following metrics were computed using the matrix 2-norm.

### Metrics

**Factorization Error**

    || P_r A P_c − L U ||₂ / ||A||₂

**Residual Error**

    || b − A x ||₂ / ||b||₂

**Growth Factor**

    || |L| |U| ||₂ / ||A||₂

**Condition Number**

    κ₂(A) = ||A||₂ · ||A⁻¹||₂

Mean and maximum values were recorded for each metric across trials.

---

## Structured Matrix Tests

The algorithm was tested on the following matrix classes:

- Diagonal matrices
- Anti-diagonal matrices
- Diagonal + anti-diagonal matrices
- Unit lower triangular matrices
- General lower triangular matrices
- Tridiagonal, diagonally dominant matrices
- Growth-factor stress matrices
- Symmetric positive definite matrices

---

## Results and Observations

- Diagonal and triangular matrices produced near-zero errors and stable growth factors
- Anti-diagonal matrices fail without pivoting but succeed with partial and complete pivoting
- Growth-factor stress matrices demonstrate the necessity of **complete pivoting**
- Ill-conditioned matrices trigger error guards or yield large errors
- Complete pivoting consistently offers the best numerical stability at the cost of speed
- Conditioning strongly correlates with stability and error behavior

---

## Files Included

- `LUdense.m`  
  In-place LU factorization with selectable pivoting strategy

- `Lvsolve.m`  
  Forward substitution solver

- `Uvsolve.m`  
  Backward substitution solver

- `TestDriverProgram3.m`  
  Runs all structured matrix tests and reports metrics

All tests complete in under 30 seconds.

---

## Key Takeaways

- Pivoting is essential for numerical stability in LU factorization
- Complete pivoting controls growth factor most effectively
- Matrix structure and conditioning strongly influence algorithm performance
- In-place storage significantly improves memory efficiency
