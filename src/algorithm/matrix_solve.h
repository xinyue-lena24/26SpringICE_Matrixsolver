/*******************************************************************************
 * File Name   : matrix_solve.h
 * Course      : C Program and Algorithm Analysis Lab -- Matrix Solver Library
 * Author      : Gong Helin
 * Year        : 2026
 * Description : Gaussian elimination and LU-based solver interfaces.
 ******************************************************************************/

#ifndef MATRIX_SOLVE_H
#define MATRIX_SOLVE_H

#include "matrix_core.h"

/*
 * Swap two rows of a matrix.
 */
MatrixError MatrixSwapRows(Matrix *A, int r1, int r2);

/*
 * Solve A x = b using Gaussian elimination with partial pivoting.
 */
MatrixError GaussianSolvePartialPivot(const Matrix *A, const Matrix *b, Matrix *x, REAL tol);

/*
 * Solve A X = B column by column using Gaussian elimination.
 */
MatrixError GaussianSolveMultiple(const Matrix *A, const Matrix *B, Matrix *X, REAL tol);

/*
 * Solve A X = B using batch Gaussian elimination.
 */
MatrixError GaussianSolveMultipleBatch(const Matrix *A, const Matrix *B, Matrix *X, REAL tol);

/*
 * Solve L U x = b from a no-pivot LU decomposition.
 */
MatrixError LUSolve(const Matrix *L, const Matrix *U, const Matrix *b, Matrix *x, REAL tol);

/*
 * Solve L U X = B from a no-pivot LU decomposition.
 */
MatrixError LUSolveMultiple(const Matrix *L, const Matrix *U, const Matrix *B, Matrix *X, REAL tol);

/*
 * Solve L U x = P b from a pivoted LU decomposition.
 */
MatrixError LUSolveWithPivot(const Matrix *L, const Matrix *U, const int *pivot, const Matrix *b, Matrix *x, REAL tol);

/*
 * Solve L U X = P B from a pivoted LU decomposition.
 */
MatrixError LUSolveMultipleWithPivot(const Matrix *L, const Matrix *U, const int *pivot, const Matrix *B, Matrix *X, REAL tol);

/*
 * Solve A X = B by first computing no-pivot LU:
 *
 *   A = L U.
 */
MatrixError LUDecomposeSolveMultiple(const Matrix *A, const Matrix *B, Matrix *X, REAL tol);

#endif
