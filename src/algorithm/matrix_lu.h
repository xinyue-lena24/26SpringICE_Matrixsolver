/*******************************************************************************
 * File Name   : matrix_lu.h
 * Course      : C Program and Algorithm Analysis Lab -- Matrix Solver Library
 * Author      : Gong Helin
 * Year        : 2026
 * Description : LU decomposition and triangular solve interfaces.
 ******************************************************************************/

#ifndef MATRIX_LU_H
#define MATRIX_LU_H

#include "matrix_core.h"

/*
 * No-pivot Doolittle LU:
 *
 *   A = L U.
 */
MatrixError LUDecomposeNoPivot(const Matrix *A, Matrix *L, Matrix *U, REAL tol);

/*
 * No-pivot LU by recording Gaussian elimination multipliers:
 *
 *   A = L U.
 */
MatrixError LUDecomposeByEliminationNoPivot(const Matrix *A, Matrix *L, Matrix *U, REAL tol);

/*
 * LU with partial pivoting:
 *
 *   P A = L U.
 *
 * The permutation matrix P is represented by pivot.
 */
MatrixError LUDecomposePartialPivot(const Matrix *A, Matrix *L, Matrix *U, int *pivot, int *swap_count, REAL tol);

/*
 * LU with partial pivoting by recording Gaussian elimination multipliers:
 *
 *   P A = L U.
 */
MatrixError LUDecomposeByEliminationPartialPivot(const Matrix *A, Matrix *L, Matrix *U, int *pivot, int *swap_count, REAL tol);

/*
 * Solve a lower triangular system:
 *
 *   L y = b.
 */
MatrixError ForwardSubstitution(const Matrix *L, const Matrix *b, Matrix *y, REAL tol);

/*
 * Solve a lower triangular system with multiple right-hand sides:
 *
 *   L Y = B.
 */
MatrixError ForwardSubstitutionMultiple(const Matrix *L, const Matrix *B, Matrix *Y, REAL tol);

/*
 * Solve an upper triangular system:
 *
 *   U x = y.
 */
MatrixError BackSubstitution(const Matrix *U, const Matrix *y, Matrix *x, REAL tol);

/*
 * Solve an upper triangular system with multiple right-hand sides:
 *
 *   U X = Y.
 */
MatrixError BackSubstitutionMultiple(const Matrix *U, const Matrix *Y, Matrix *X, REAL tol);

/*
 * Compute determinant from the U factor of a no-pivot LU decomposition.
 */
MatrixError LUDeterminant(const Matrix *U, REAL *det);

#endif