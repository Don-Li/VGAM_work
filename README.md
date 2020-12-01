# VGAM work

## Introduction

We have been given a quest by the VGAM Lord Thomas Yee to improve some of the matrix multiplication in VGAM.

The particular problem is transforming a band-format matrix, doing Cholesky decomposition for each of the resulting square matrices, and multiplying it with a design matrix.

## Results

We have achieved ~10x speedup by rewriting the matrix multiplication function with `Rcpp` and `armadillo`.

`Unit: microseconds`

| ## | expr  | min  | lq   | mean     | median | uq   | max    | neval |
|----|-------|------|------|----------|--------|------|--------|-------|
| ## | VGAM  | 36.2 | 41.1 | **50.76489** | 42.5   | 45.2 | 4353.5 | 10000 |
| ## | cpp1  | 3.5  | 4.3  | **5.92894**  | 4.7    | 5.5  | 3508.4 | 1000  |
| ## | cpp2  | 6.1  | 7.0  | 9.41595  | 7.3    | 8.4  | 3564.6 | 1000  |
| ## | cpp3  | 7.4  | 8.4  | 10.79717 | 8.9    | 10.2 | 2273.1 | 1000  |
| ## | cpp4  | 7.4  | 8.3  | 10.64852 | 8.8    | 10.0 | 2840.7 | 1000  |


*  **Old**: 50.76489

*  **New**: 5.92894

## Implementation

For full implementation details, please see the documentation and testing results at https://don-li.github.io/VGAM_work .