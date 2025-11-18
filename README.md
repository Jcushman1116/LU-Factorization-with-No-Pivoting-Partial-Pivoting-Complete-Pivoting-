✅ README Template: Program 3 — LU Factorization (No, Partial, Complete Pivoting)
📌 Overview

This project develops a unified LU factorization function supporting:

No pivoting

Partial pivoting

Complete pivoting

The function performs in-place factorization and stores permutation information in Prow and Pcol. Comprehensive empirical testing evaluates accuracy, growth factor, and stability.

🧠 Key Concepts

Gaussian elimination

Pivoting strategies

Permutation matrices (row & column)

Growth factor analysis

Error guards for near-singular matrices

Backsubstitution and reconstruction

⚙️ Implemented Methods
1. No Pivoting (A = LU)

Fastest but numerically risky.

2. Partial Pivoting (PA = LU)

Stable for most matrices; industry standard.

3. Complete Pivoting (PAQ = LU)

Most stable but slowest.

🧪 Structured Matrix Tests

The algorithm is tested on:

Diagonal matrices

Anti-diagonal matrices

Diagonal + anti-diagonal

Unit lower triangular

General lower triangular

Growth factor matrices

Metrics:

Factorization error

Residual error

Growth factor γ

Conditioning number κ2(A)

📈 Results Summary

No pivoting fails on singular / anti-diagonal cases (as expected).

Partial and complete pivoting handle all structured tests.

Growth factor aligns with theory:

Explodes for LU without pivoting

Controlled for partial

Smallest for complete

Error metrics near machine precision except in ill-conditioned constructions.
