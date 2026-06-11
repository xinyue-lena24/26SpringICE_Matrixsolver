#include "matrix_solve.h"
#include "matrix_lu.h"

#include <math.h>
#include <stdlib.h>

/*
 * Swap two rows of a matrix.
 */
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

/*
 * Check A x = b.
 */
static MatrixError CheckLinearSystem(const Matrix *A, const Matrix *b, Matrix *x)
{
    if (!MatrixIsValid(A) || !MatrixIsValid(b) || !MatrixIsValid(x)) {
        return MATRIX_ERROR_NULL_POINTER;
    }

    if (A->row != A->column) {
        return MATRIX_ERROR_NOT_SQUARE;
    }

    int n = A->row;
    if (!MatrixHasShape(b, n, 1) || !MatrixHasShape(x, n, 1)) {
        return MATRIX_ERROR_SIZE_MISMATCH;
    }

    return MATRIX_SUCCESS;
}

/*
 * Check A X = B.
 */
static MatrixError CheckLinearSystemMultiple(const Matrix *A, const Matrix *B, Matrix *X)
{
    if (!MatrixIsValid(A) || !MatrixIsValid(B) || !MatrixIsValid(X)) {
        return MATRIX_ERROR_NULL_POINTER;
    }

    if (A->row != A->column) {
        return MATRIX_ERROR_NOT_SQUARE;
    }

    int n = A->row;
    if (B->row != n || X->row != n || X->column != B->column) {
        return MATRIX_ERROR_SIZE_MISMATCH;
    }

    return MATRIX_SUCCESS;
}

/*
 * Check L U x = b.
 */
static MatrixError CheckLUSystem(const Matrix *L, const Matrix *U, const Matrix *b, Matrix *x)
{
    if (!MatrixIsValid(L) || !MatrixIsValid(U) || !MatrixIsValid(b) || !MatrixIsValid(x)) {
        return MATRIX_ERROR_NULL_POINTER;
    }

    if (L->row != L->column || U->row != U->column) {
        return MATRIX_ERROR_NOT_SQUARE;
    }

    if (L->row != U->row || L->column != U->column) {
        return MATRIX_ERROR_SIZE_MISMATCH;
    }

    int n = L->row;
    if (!MatrixHasShape(b, n, 1) || !MatrixHasShape(x, n, 1)) {
        return MATRIX_ERROR_SIZE_MISMATCH;
    }

    return MATRIX_SUCCESS;
}

/*
 * Check L U X = B.
 */
static MatrixError CheckLUSystemMultiple(const Matrix *L, const Matrix *U, const Matrix *B, Matrix *X)
{
    if (!MatrixIsValid(L) || !MatrixIsValid(U) || !MatrixIsValid(B) || !MatrixIsValid(X)) {
        return MATRIX_ERROR_NULL_POINTER;
    }

    if (L->row != L->column || U->row != U->column) {
        return MATRIX_ERROR_NOT_SQUARE;
    }

    if (L->row != U->row || L->column != U->column) {
        return MATRIX_ERROR_SIZE_MISMATCH;
    }

    int n = L->row;
    if (B->row != n || X->row != n || X->column != B->column) {
        return MATRIX_ERROR_SIZE_MISMATCH;
    }

    return MATRIX_SUCCESS;
}

/*
 * Check pivot array.
 */
static MatrixError CheckPivotArray(const int *pivot, int n)
{
    if (pivot == NULL) {
        return MATRIX_ERROR_NULL_POINTER;
    }

    for (int i = 0; i < n; ++i) {
        if (pivot[i] < 0 || pivot[i] >= n) {
            return MATRIX_ERROR_INDEX_OUT_OF_RANGE;
        }
    }

    return MATRIX_SUCCESS;
}

/*
 * Solve A x = b using Gaussian elimination with partial pivoting.
 */
MatrixError GaussianSolveVector(const Matrix *A, const Matrix *b, Matrix *x, REAL tol)
{
    MatrixError error = CheckLinearSystem(A, b, x);
    if (error != MATRIX_SUCCESS) {
        return error;
    }

    int n = A->row;

    Matrix U, rhs;
    MatrixInit(&U);
    MatrixInit(&rhs);

    error = MatrixCreate(&U, n, n);
    if (error != MATRIX_SUCCESS) {
        return error;
    }

    error = MatrixCreate(&rhs, n, 1);
    if (error != MATRIX_SUCCESS) {
        MatrixFree(&U);
        return error;
    }

    error = MatrixCopy(A, &U);
    if (error != MATRIX_SUCCESS) {
        MatrixFree(&U);
        MatrixFree(&rhs);
        return error;
    }

    error = MatrixCopy(b, &rhs);
    if (error != MATRIX_SUCCESS) {
        MatrixFree(&U);
        MatrixFree(&rhs);
        return error;
    }

    for (int k = 0; k < n - 1; ++k) {
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
            MatrixFree(&U);
            MatrixFree(&rhs);
            return MATRIX_ERROR_SINGULAR;
        }

        error = MatrixSwapRows(&U, k, pivot_row);
        if (error != MATRIX_SUCCESS) {
            MatrixFree(&U);
            MatrixFree(&rhs);
            return error;
        }

        error = MatrixSwapRows(&rhs, k, pivot_row);
        if (error != MATRIX_SUCCESS) {
            MatrixFree(&U);
            MatrixFree(&rhs);
            return error;
        }

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

    error = BackSubstitution(&U, &rhs, x, tol);

    MatrixFree(&U);
    MatrixFree(&rhs);
    return error;
}

/*
 * Solve A X = B by solving each column separately.
 */
MatrixError GaussianSolveMatrix(const Matrix *A, const Matrix *B, Matrix *X, REAL tol)
{
    MatrixError error = CheckLinearSystemMultiple(A, B, X);
    if (error != MATRIX_SUCCESS) {
        return error;
    }

    Matrix b_col, x_col;
    MatrixInit(&b_col);
    MatrixInit(&x_col);

    error = MatrixCreate(&b_col, B->row, 1);
    if (error != MATRIX_SUCCESS) {
        return error;
    }

    error = MatrixCreate(&x_col, X->row, 1);
    if (error != MATRIX_SUCCESS) {
        MatrixFree(&b_col);
        return error;
    }

    for (int col = 0; col < B->column; ++col) {
        for (int i = 0; i < B->row; ++i) {
            b_col.data[i] = B->data[MatrixIndex(B, i, col)];
        }

        error = GaussianSolveVector(A, &b_col, &x_col, tol);
        if (error != MATRIX_SUCCESS) {
            MatrixFree(&b_col);
            MatrixFree(&x_col);
            return error;
        }

        for (int i = 0; i < X->row; ++i) {
            X->data[MatrixIndex(X, i, col)] = x_col.data[i];
        }
    }

    MatrixFree(&b_col);
    MatrixFree(&x_col);
    return MATRIX_SUCCESS;
}

/*
 * Solve A X = B by batch Gaussian elimination.
 */
MatrixError GaussianSolveMatrixBatch(const Matrix *A, const Matrix *B, Matrix *X, REAL tol)
{
    MatrixError error = CheckLinearSystemMultiple(A, B, X);
    if (error != MATRIX_SUCCESS) {
        return error;
    }

    int n = A->row;
    int nrhs = B->column;

    Matrix U, RHS;
    MatrixInit(&U);
    MatrixInit(&RHS);

    error = MatrixCreate(&U, n, n);
    if (error != MATRIX_SUCCESS) {
        return error;
    }

    error = MatrixCreate(&RHS, n, nrhs);
    if (error != MATRIX_SUCCESS) {
        MatrixFree(&U);
        return error;
    }

    error = MatrixCopy(A, &U);
    if (error != MATRIX_SUCCESS) {
        MatrixFree(&U);
        MatrixFree(&RHS);
        return error;
    }

    error = MatrixCopy(B, &RHS);
    if (error != MATRIX_SUCCESS) {
        MatrixFree(&U);
        MatrixFree(&RHS);
        return error;
    }

    for (int k = 0; k < n - 1; ++k) {
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
            MatrixFree(&U);
            MatrixFree(&RHS);
            return MATRIX_ERROR_SINGULAR;
        }

        error = MatrixSwapRows(&U, k, pivot_row);
        if (error != MATRIX_SUCCESS) {
            MatrixFree(&U);
            MatrixFree(&RHS);
            return error;
        }

        error = MatrixSwapRows(&RHS, k, pivot_row);
        if (error != MATRIX_SUCCESS) {
            MatrixFree(&U);
            MatrixFree(&RHS);
            return error;
        }

        REAL pivot = U.data[MatrixIndex(&U, k, k)];

        for (int i = k + 1; i < n; ++i) {
            REAL factor = U.data[MatrixIndex(&U, i, k)] / pivot;

            U.data[MatrixIndex(&U, i, k)] = 0.0;

            for (int j = k + 1; j < n; ++j) {
                U.data[MatrixIndex(&U, i, j)] -= factor * U.data[MatrixIndex(&U, k, j)];
            }

            for (int col = 0; col < nrhs; ++col) {
                RHS.data[MatrixIndex(&RHS, i, col)] -= factor * RHS.data[MatrixIndex(&RHS, k, col)];
            }
        }
    }

    error = BackSubstitutionMultiple(&U, &RHS, X, tol);

    MatrixFree(&U);
    MatrixFree(&RHS);
    return error;
}

/*
 * Solve L U x = b from no-pivot LU.
 */
MatrixError LUSolveVector(const Matrix *L, const Matrix *U, const Matrix *b, Matrix *x, REAL tol)
{
    MatrixError error = CheckLUSystem(L, U, b, x);
    if (error != MATRIX_SUCCESS) {
        return error;
    }

    Matrix y;
    MatrixInit(&y);

    error = MatrixCreate(&y, b->row, 1);
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

/*
 * Solve L U X = B from no-pivot LU.
 */
MatrixError LUSolveMatrix(const Matrix *L, const Matrix *U, const Matrix *B, Matrix *X, REAL tol)
{
    MatrixError error = CheckLUSystemMultiple(L, U, B, X);
    if (error != MATRIX_SUCCESS) {
        return error;
    }

    Matrix Y;
    MatrixInit(&Y);

    error = MatrixCreate(&Y, B->row, B->column);
    if (error != MATRIX_SUCCESS) {
        return error;
    }

    error = ForwardSubstitutionMultiple(L, B, &Y, tol);
    if (error == MATRIX_SUCCESS) {
        error = BackSubstitutionMultiple(U, &Y, X, tol);
    }

    MatrixFree(&Y);
    return error;
}

/*
 * Solve L U x = P b from pivoted LU.
 */
MatrixError LUPivotSolveVector(const Matrix *L, const Matrix *U, const int *pivot, const Matrix *b, Matrix *x, REAL tol)
{
    MatrixError error = CheckLUSystem(L, U, b, x);
    if (error != MATRIX_SUCCESS) {
        return error;
    }

    int n = b->row;

    error = CheckPivotArray(pivot, n);
    if (error != MATRIX_SUCCESS) {
        return error;
    }

    Matrix Pb, y;
    MatrixInit(&Pb);
    MatrixInit(&y);

    error = MatrixCreate(&Pb, n, 1);
    if (error != MATRIX_SUCCESS) {
        return error;
    }

    error = MatrixCreate(&y, n, 1);
    if (error != MATRIX_SUCCESS) {
        MatrixFree(&Pb);
        return error;
    }

    for (int i = 0; i < n; ++i) {
        Pb.data[i] = b->data[pivot[i]];
    }

    error = ForwardSubstitution(L, &Pb, &y, tol);
    if (error == MATRIX_SUCCESS) {
        error = BackSubstitution(U, &y, x, tol);
    }

    MatrixFree(&Pb);
    MatrixFree(&y);
    return error;
}

/*
 * Solve L U X = P B from pivoted LU.
 */
MatrixError LUPivotSolveMatrix(const Matrix *L, const Matrix *U, const int *pivot, const Matrix *B, Matrix *X, REAL tol)
{
    MatrixError error = CheckLUSystemMultiple(L, U, B, X);
    if (error != MATRIX_SUCCESS) {
        return error;
    }

    int n = B->row;
    int nrhs = B->column;

    error = CheckPivotArray(pivot, n);
    if (error != MATRIX_SUCCESS) {
        return error;
    }

    Matrix PB, Y;
    MatrixInit(&PB);
    MatrixInit(&Y);

    error = MatrixCreate(&PB, n, nrhs);
    if (error != MATRIX_SUCCESS) {
        return error;
    }

    error = MatrixCreate(&Y, n, nrhs);
    if (error != MATRIX_SUCCESS) {
        MatrixFree(&PB);
        return error;
    }

    for (int i = 0; i < n; ++i) {
        for (int col = 0; col < nrhs; ++col) {
            PB.data[MatrixIndex(&PB, i, col)] = B->data[MatrixIndex(B, pivot[i], col)];
        }
    }

    error = ForwardSubstitutionMultiple(L, &PB, &Y, tol);
    if (error == MATRIX_SUCCESS) {
        error = BackSubstitutionMultiple(U, &Y, X, tol);
    }

    MatrixFree(&PB);
    MatrixFree(&Y);
    return error;
}

/*
 * Solve A X = B by first computing no-pivot LU.
 */
MatrixError LUFactorSolveMatrix(const Matrix *A, const Matrix *B, Matrix *X, REAL tol)
{
    MatrixError error = CheckLinearSystemMultiple(A, B, X);
    if (error != MATRIX_SUCCESS) {
        return error;
    }

    int n = A->row;

    Matrix L, U;
    MatrixInit(&L);
    MatrixInit(&U);

    error = MatrixCreate(&L, n, n);
    if (error != MATRIX_SUCCESS) {
        return error;
    }

    error = MatrixCreate(&U, n, n);
    if (error != MATRIX_SUCCESS) {
        MatrixFree(&L);
        return error;
    }

    error = LUDecomposeNoPivot(A, &L, &U, tol);
    if (error == MATRIX_SUCCESS) {
        error = LUSolveMatrix(&L, &U, B, X, tol);
    }

    MatrixFree(&L);
    MatrixFree(&U);
    return error;
}

MatrixError LUFactorPivotSolveMatrix(const Matrix *A, const Matrix *B, Matrix *X, REAL tol)
{
    MatrixError error = CheckLinearSystemMultiple(A, B, X);
    if (error != MATRIX_SUCCESS) {
        return error;
    }

    int n = A->row;

    Matrix L, U;
    MatrixInit(&L);
    MatrixInit(&U);

    int *pivot = (int *)malloc(n * sizeof(int));
    if (pivot == NULL) {
        return MATRIX_ERROR_ALLOC_FAILED;
    }

    error = MatrixCreate(&L, n, n);
    if (error != MATRIX_SUCCESS) {
        free(pivot);
        return error;
    }

    error = MatrixCreate(&U, n, n);
    if (error != MATRIX_SUCCESS) {
        free(pivot);
        MatrixFree(&L);
        return error;
    }

    int swap_count;
    error = LUDecomposePartialPivot(A, &L, &U, pivot, &swap_count, tol);
    if (error == MATRIX_SUCCESS) {
        error = LUPivotSolveMatrix(&L, &U, pivot, B, X, tol);
    }

    free(pivot);
    MatrixFree(&L);
    MatrixFree(&U);
    return error;
}