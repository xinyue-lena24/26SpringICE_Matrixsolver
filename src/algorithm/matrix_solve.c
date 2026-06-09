/*******************************************************************************
* File Name   : matrix_solve.c
* Course      : C Program and Algorithm Analysis Lab -- Matrix Solver Library
* Author      : Gong Helin
* Year        : 2026
* Description : Gaussian elimination 线性方程组求解函数实现。
*******************************************************************************/

#include "matrix_solve.h"
#include "matrix_lu.h"
#include <math.h>

MatrixError MatrixSwapRows(Matrix *A, int r1, int r2)
{
    if (!MatrixIsValid(A)) {
        return MATRIX_ERROR_NULL_POINTER;
    }
    if (r1 < 0 || r1 >= A->row || r2 < 0 || r2 >= A->row) {
        return MATRIX_ERROR_INDEX_OUT_OF_RANGE;
    }
    if (r1 == r2) {
        return MATRIX_SUCCESS;
    }
    for (int j = 0; j < A->column; ++j) {
        int idx1 = MatrixIndex(A, r1, j);
        int idx2 = MatrixIndex(A, r2, j);
        REAL tmp = A->data[idx1];
        A->data[idx1] = A->data[idx2];
        A->data[idx2] = tmp;
    }
    return MATRIX_SUCCESS;
}

static MatrixError CheckLinearSystem(const Matrix *A, const Matrix *b, Matrix *x)
{
    if (!MatrixIsValid(A) || !MatrixIsValid(b) || !MatrixIsValid(x)) {
        return MATRIX_ERROR_NULL_POINTER;
    }
    if (A->row != A->column) {
        return MATRIX_ERROR_NOT_SQUARE;
    }
    if (b->row != A->row || b->column != 1) {
        return MATRIX_ERROR_SIZE_MISMATCH;
    }
    if (x->row != A->column || x->column != 1) {
        return MATRIX_ERROR_SIZE_MISMATCH;
    }
    return MATRIX_SUCCESS;
}

static MatrixError CheckLinearSystemMultiple(const Matrix *A, const Matrix *B, Matrix *X)
{
    if (!MatrixIsValid(A) || !MatrixIsValid(B) || !MatrixIsValid(X)) {
        return MATRIX_ERROR_NULL_POINTER;
    }
    if (A->row != A->column) {
        return MATRIX_ERROR_NOT_SQUARE;
    }
    if (B->row != A->row) {
        return MATRIX_ERROR_SIZE_MISMATCH;
    }
    if (X->row != A->column || X->column != B->column) {
        return MATRIX_ERROR_SIZE_MISMATCH;
    }
    return MATRIX_SUCCESS;
}

MatrixError GaussianSolvePartialPivot(const Matrix *A, const Matrix *b, Matrix *x, REAL tol)
{
    MatrixError error = CheckLinearSystem(A, b, x);
    if (error != MATRIX_SUCCESS) {
        return error;
    }
    int n = A->row;
    Matrix U, rhs;
    MatrixInit(&U); MatrixInit(&rhs);
    MatrixCreate(&U, n, n);
    MatrixCreate(&rhs, n, 1);
    MatrixCopy(A, &U);
    MatrixCopy(b, &rhs);

    // Turn U into an upper triangular matrix with partial pivoting
    for (int k = 0; k < n - 1; ++k) {
        // Find the pivot row with the largest absolute value in the current column k
        int pivot_row = k;
        REAL pivot_abs = fabs(U.data[MatrixIndex(&U, k, k)]);
        for (int i = k + 1; i < n; ++i) {
            REAL value_abs = fabs(U.data[MatrixIndex(&U, i, k)]);
            if (value_abs > pivot_abs) {
                pivot_abs = value_abs;
                pivot_row = i;
            }
        }
        if (pivot_abs < tol) {
            MatrixFree(&U); MatrixFree(&rhs);
            return MATRIX_ERROR_SINGULAR;
        }
        MatrixSwapRows(&U, k, pivot_row);
        MatrixSwapRows(&rhs, k, pivot_row);

        // Perform elimination below the pivot
        REAL pivot = U.data[MatrixIndex(&U, k, k)];
        for (int i = k + 1; i < n; ++i) {
            REAL factor = U.data[MatrixIndex(&U, i, k)] / pivot;
            U.data[MatrixIndex(&U, i, k)] = 0.0;
            for (int j = k + 1; j < n; ++j) {
                U.data[MatrixIndex(&U, i, j)] -= factor * U.data[MatrixIndex(&U, k, j)];
            }
            rhs.data[i] -= factor * rhs.data[k];
        }
    }

    // Back substitution to solve Ux = rhs
    for (int i = n - 1; i >= 0; --i) {
        REAL sum = rhs.data[i];
        for (int j = i + 1; j < n; ++j) {
            sum -= U.data[MatrixIndex(&U, i, j)] * x->data[j];
        }
        REAL diag = U.data[MatrixIndex(&U, i, i)];
        if (fabs(diag) < tol) {
            MatrixFree(&U); MatrixFree(&rhs);
            return MATRIX_ERROR_SINGULAR;
        }
        x->data[i] = sum / diag;
    }

    MatrixFree(&U); MatrixFree(&rhs);
    return MATRIX_SUCCESS;
}


// Solves Ax = b using the LU decomposition of A (L and U) without pivoting
MatrixError LUSolve(const Matrix *L, const Matrix *U, const Matrix *b, Matrix *x, REAL tol)
{
    Matrix y;
    MatrixInit(&y);
    MatrixError error = MatrixCreate(&y, b->row, 1);
    if (error != MATRIX_SUCCESS) {
        return error;
    }
    error = ForwardSubstitution(L, b, &y, tol);
    if (error == MATRIX_SUCCESS) {
        error = BackSubstitution(U, &y, x, tol);
    }
    MatrixFree(&y);
    return error;
}

// Solves AX = B for multiple right-hand sides
MatrixError LUDecomposeSolveMultiple(const Matrix *A, const Matrix *B, Matrix *X, REAL tol)
{
    MatrixError error = CheckLinearSystemMultiple(A, B, X);
    if (error != MATRIX_SUCCESS) {
        return error;
    }
    int n = A->row;
    Matrix L, U;
    MatrixInit(&L); MatrixInit(&U);
    MatrixCreate(&L, n, n);
    MatrixCreate(&U, n, n);

    error = LUDecomposeNoPivot(A, &L, &U, tol);
    if (error != MATRIX_SUCCESS) {
        MatrixFree(&L); MatrixFree(&U);
        return error;
    }
    
    // Solve columes of X one by one
    for(int col = 0; col < B->column; ++col) {
        Matrix b_col, x_col;
        MatrixInit(&b_col); MatrixInit(&x_col);
        MatrixCreate(&b_col, n, 1);
        MatrixCreate(&x_col, n, 1);
        
        // Extract the current column of B as the right-hand side vector
        for (int i = 0; i < n; ++i) {
            b_col.data[i] = B->data[MatrixIndex(B, i, col)];
        }

        // Solve L(Ux) = b_col
        error = LUSolve(&L, &U, &b_col, &x_col, tol);
        if (error != MATRIX_SUCCESS) {
            MatrixFree(&L); MatrixFree(&U);
            MatrixFree(&b_col); MatrixFree(&x_col);
            return error;
        }

        // Copy solution to X
        for (int i = 0; i < n; ++i) {
            X->data[MatrixIndex(X, i, col)] = x_col.data[i];
        }
        MatrixFree(&b_col); MatrixFree(&x_col);
    }
    MatrixFree(&L); MatrixFree(&U);
    return error;
}