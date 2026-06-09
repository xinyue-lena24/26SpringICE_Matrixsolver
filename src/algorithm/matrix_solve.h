/*******************************************************************************
* File Name   : matrix_solve.h
* Course      : C Program and Algorithm Analysis Lab -- Matrix Solver Library
* Author      : Gong Helin
* Year        : 2026
* Description : Gaussian elimination 求解接口声明。
*******************************************************************************/

#ifndef MATRIX_SOLVE_H
#define MATRIX_SOLVE_H

#include "matrix_core.h"

MatrixError MatrixSwapRows(Matrix *A, int r1, int r2);
MatrixError GaussianSolvePartialPivot(const Matrix *A, const Matrix *b, Matrix *x, REAL tol);
MatrixError GaussianSolveMultiple(const Matrix *A, const Matrix *B, Matrix *X, REAL tol);
MatrixError LUSolve(const Matrix *L, const Matrix *U, const Matrix *b, Matrix *x, REAL tol);
MatrixError LUDecomposeSolveMultiple(const Matrix *A, const Matrix *B, Matrix *X, REAL tol);

#endif
