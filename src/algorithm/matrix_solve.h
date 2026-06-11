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
MatrixError GaussianSolveVector(const Matrix *A, const Matrix *b, Matrix *x, REAL tol);

/*
 * Solve A X = B column by column using Gaussian elimination.
 */
MatrixError GaussianSolveMatrix(const Matrix *A, const Matrix *B, Matrix *X, REAL tol);

/*
 * Solve A X = B using batch Gaussian elimination.
 */
MatrixError GaussianSolveMatrixBatch(const Matrix *A, const Matrix *B, Matrix *X, REAL tol);

/*
 * Solve L U x = b from a no-pivot LU decomposition.
 */
MatrixError LUSolveVector(const Matrix *L, const Matrix *U, const Matrix *b, Matrix *x, REAL tol);

/*
 * Solve L U X = B from a no-pivot LU decomposition.
 */
MatrixError LUSolveMatrix(const Matrix *L, const Matrix *U, const Matrix *B, Matrix *X, REAL tol);

/*
 * Solve L U x = P b from a pivoted LU decomposition.
 */
MatrixError LUPivotSolveVector(const Matrix *L, const Matrix *U, const int *pivot, const Matrix *b, Matrix *x, REAL tol);

/*
 * Solve L U X = P B from a pivoted LU decomposition.
 */
MatrixError LUPivotSolveMatrix(const Matrix *L, const Matrix *U, const int *pivot, const Matrix *B, Matrix *X, REAL tol);

/*
 * Solve A X = B by first computing no-pivot LU:
 *
 *   A = L U.
 */
MatrixError LUFactorSolvMatrix(const Matrix *A, const Matrix *B, Matrix *X, REAL tol);

#endif
