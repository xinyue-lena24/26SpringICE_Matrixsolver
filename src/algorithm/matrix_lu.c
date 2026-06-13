#include "matrix_lu.h"

#include <math.h>

/*
 * Check whether a matrix is valid and square.
 */
static MatrixError CheckSquareMatrix(const Matrix *A)
{
    if (!MatrixIsValid(A)) {
        return MATRIX_ERROR_NULL_POINTER;
    }

    if (A->row != A->column) {
        return MATRIX_ERROR_NOT_SQUARE;
    }

    return MATRIX_SUCCESS;
}

/*
 * Check the input dimensions for LU decomposition.
 */
static MatrixError CheckLUDecomposeInput(const Matrix *A, Matrix *L, Matrix *U)
{
    MatrixError error = CheckSquareMatrix(A);
    if (error != MATRIX_SUCCESS) {
        return error;
    }

    if (!MatrixIsValid(L) || !MatrixIsValid(U)) {
        return MATRIX_ERROR_NULL_POINTER;
    }

    int n = A->row;
    if (!MatrixHasShape(L, n, n) || !MatrixHasShape(U, n, n)) {
        return MATRIX_ERROR_SIZE_MISMATCH;
    }

    return MATRIX_SUCCESS;
}

/*
 * Check the input dimensions for a triangular system with one right-hand side:
 *
 *   T x = b.
 */
static MatrixError CheckTriangularSolveInput(const Matrix *T, const Matrix *b, Matrix *x)
{
    MatrixError error = CheckSquareMatrix(T);
    if (error != MATRIX_SUCCESS) {
        return error;
    }

    if (!MatrixIsValid(b) || !MatrixIsValid(x)) {
        return MATRIX_ERROR_NULL_POINTER;
    }

    int n = T->row;
    if (!MatrixHasShape(b, n, 1) || !MatrixHasShape(x, n, 1)) {
        return MATRIX_ERROR_SIZE_MISMATCH;
    }

    return MATRIX_SUCCESS;
}

/*
 * Check the input dimensions for a triangular system with multiple right-hand sides:
 *
 *   T X = B.
 */
static MatrixError CheckTriangularSolveMultipleInput(const Matrix *T, const Matrix *B, Matrix *X)
{
    MatrixError error = CheckSquareMatrix(T);
    if (error != MATRIX_SUCCESS) {
        return error;
    }

    if (!MatrixIsValid(B) || !MatrixIsValid(X)) {
        return MATRIX_ERROR_NULL_POINTER;
    }

    int n = T->row;
    if (B->row != n || X->row != n || X->column != B->column) {
        return MATRIX_ERROR_SIZE_MISMATCH;
    }

    return MATRIX_SUCCESS;
}

/*
 * Local row-swap helper for the LU module.
 *
 * This avoids making matrix_lu.c depend on matrix_solve.c.
 */
static MatrixError MatrixSwapRowsLocal(Matrix *A, int r1, int r2)
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
 * Doolittle LU decomposition without pivoting:
 *
 *   A = L U.
 *
 * L is unit lower triangular and U is upper triangular.
 */
MatrixError LUDecomposeNoPivotSimple(const Matrix *A, Matrix *L, Matrix *U, REAL tol)
{
    MatrixError error = CheckLUDecomposeInput(A, L, U);
    if (error != MATRIX_SUCCESS) {
        return error;
    }

    int n = A->row;

    error = MatrixFillZero(L);
    if (error != MATRIX_SUCCESS) {
        return error;
    }

    error = MatrixFillZero(U);
    if (error != MATRIX_SUCCESS) {
        return error;
    }

    for (int i = 0; i < n; ++i) {
        L->data[MatrixIndex(L, i, i)] = 1.0;
    }

    for (int k = 0; k < n; ++k) {
        for (int j = k; j < n; ++j) {
            REAL sum = 0.0;

            for (int s = 0; s < k; ++s) {
                sum += L->data[MatrixIndex(L, k, s)] * U->data[MatrixIndex(U, s, j)];
            }

            U->data[MatrixIndex(U, k, j)] = A->data[MatrixIndex(A, k, j)] - sum;
        }

        REAL pivot = U->data[MatrixIndex(U, k, k)];
        if (fabs(pivot) < tol) {
            return MATRIX_ERROR_SINGULAR;
        }

        for (int i = k + 1; i < n; ++i) {
            REAL sum = 0.0;

            for (int s = 0; s < k; ++s) {
                sum += L->data[MatrixIndex(L, i, s)] * U->data[MatrixIndex(U, s, k)];
            }

            L->data[MatrixIndex(L, i, k)] = (A->data[MatrixIndex(A, i, k)] - sum) / pivot;
        }
    }

    return MATRIX_SUCCESS;
}

/*
 * No-pivot LU by recording Gaussian elimination multipliers:
 *
 *   A = L U.
 *
 * This version makes explicit that L stores the elimination factors.
 */
MatrixError LUDecomposeByEliminationNoPivot(const Matrix *A, Matrix *L, Matrix *U, REAL tol)
{
    MatrixError error = CheckLUDecomposeInput(A, L, U);
    if (error != MATRIX_SUCCESS) {
        return error;
    }

    int n = A->row;

    error = MatrixFillIdentity(L);
    if (error != MATRIX_SUCCESS) {
        return error;
    }

    error = MatrixCopy(A, U);
    if (error != MATRIX_SUCCESS) {
        return error;
    }

    for (int k = 0; k < n; ++k) {
        REAL pivot_value = U->data[MatrixIndex(U, k, k)];
        if (fabs(pivot_value) < tol) {
            return MATRIX_ERROR_SINGULAR;
        }

        for (int i = k + 1; i < n; ++i) {
            REAL factor = U->data[MatrixIndex(U, i, k)] / pivot_value;

            L->data[MatrixIndex(L, i, k)] = factor;
            U->data[MatrixIndex(U, i, k)] = 0.0;

            int idx_i = MatrixIndex(U, i, 0);
            int idx_k = MatrixIndex(U, k, 0);
            for (int j = k + 1; j < n; ++j) {
                U->data[idx_i + j] -= factor * U->data[idx_k + j];
            }
        }
    }

    return MATRIX_SUCCESS;
}

/*
 * Public interface for no-pivot LU:
 *
 *   A = L U.
 */
MatrixError LUDecomposeNoPivot(const Matrix *A, Matrix *L, Matrix *U, REAL tol)
{
    return LUDecomposeByEliminationNoPivot(A, L, U, tol);
}

/*
 * LU with partial pivoting by recording Gaussian elimination multipliers:
 *
 *   P A = L U.
 *
 * pivot[i] stores the original row index currently placed at row i.
 * swap_count stores the number of row swaps.
 */
MatrixError LUDecomposeByEliminationPartialPivot(const Matrix *A, Matrix *L, Matrix *U, int *pivot, int *swap_count, REAL tol)
{
    MatrixError error = CheckLUDecomposeInput(A, L, U);
    if (error != MATRIX_SUCCESS) {
        return error;
    }

    if (pivot == NULL || swap_count == NULL) {
        return MATRIX_ERROR_NULL_POINTER;
    }

    int n = A->row;

    for (int i = 0; i < n; ++i) {
        pivot[i] = i;
    }
    *swap_count = 0;

    error = MatrixFillIdentity(L);
    if (error != MATRIX_SUCCESS) {
        return error;
    }

    error = MatrixCopy(A, U);
    if (error != MATRIX_SUCCESS) {
        return error;
    }

    for (int k = 0; k < n; ++k) {
        int pivot_row = k;
        REAL pivot_abs = fabs(U->data[MatrixIndex(U, k, k)]);

        for (int i = k + 1; i < n; ++i) {
            REAL value_abs = fabs(U->data[MatrixIndex(U, i, k)]);

            if (value_abs > pivot_abs) {
                pivot_abs = value_abs;
                pivot_row = i;
            }
        }

        if (pivot_abs < tol) {
            return MATRIX_ERROR_SINGULAR;
        }

        if (pivot_row != k) {
            error = MatrixSwapRowsLocal(U, k, pivot_row);
            if (error != MATRIX_SUCCESS) {
                return error;
            }

            for (int j = 0; j < k; ++j) {
                REAL tmp = L->data[MatrixIndex(L, k, j)];
                L->data[MatrixIndex(L, k, j)] = L->data[MatrixIndex(L, pivot_row, j)];
                L->data[MatrixIndex(L, pivot_row, j)] = tmp;
            }

            int tmp_pivot = pivot[k];
            pivot[k] = pivot[pivot_row];
            pivot[pivot_row] = tmp_pivot;

            ++(*swap_count);
        }

        REAL pivot_value = U->data[MatrixIndex(U, k, k)];

        for (int i = k + 1; i < n; ++i) {
            REAL factor = U->data[MatrixIndex(U, i, k)] / pivot_value;

            L->data[MatrixIndex(L, i, k)] = factor;
            U->data[MatrixIndex(U, i, k)] = 0.0;

            int idx_i = MatrixIndex(U, i, 0);
            int idx_k = MatrixIndex(U, k, 0);
            #pragma omp simd
            for (int j = k + 1; j < n; ++j) {
                U->data[idx_i + j] -= factor * U->data[idx_k + j];
            }
        }
    }

    return MATRIX_SUCCESS;
}

/*
 * Public interface for LU with partial pivoting:
 *
 *   P A = L U.
 */
MatrixError LUDecomposePartialPivot(const Matrix *A, Matrix *L, Matrix *U, int *pivot, int *swap_count, REAL tol)
{
    return LUDecomposeByEliminationPartialPivot(A, L, U, pivot, swap_count, tol);
}

/*
 * Forward substitution:
 *
 *   L y = b.
 */
MatrixError ForwardSubstitution(const Matrix *L, const Matrix *b, Matrix *y, REAL tol)
{
    MatrixError error = CheckTriangularSolveInput(L, b, y);
    if (error != MATRIX_SUCCESS) {
        return error;
    }

    int n = L->row;

    for (int i = 0; i < n; ++i) {
        REAL sum = b->data[i];

        for (int j = 0; j < i; ++j) {
            sum -= L->data[MatrixIndex(L, i, j)] * y->data[j];
        }

        REAL diag = L->data[MatrixIndex(L, i, i)];
        if (fabs(diag) < tol) {
            return MATRIX_ERROR_SINGULAR;
        }

        y->data[i] = sum / diag;
    }

    return MATRIX_SUCCESS;
}

/*
 * Forward substitution with multiple right-hand sides:
 *
 *   L Y = B.
 */
MatrixError ForwardSubstitutionMultiple(const Matrix *L, const Matrix *B, Matrix *Y, REAL tol)
{
    MatrixError error = CheckTriangularSolveMultipleInput(L, B, Y);
    if (error != MATRIX_SUCCESS) {
        return error;
    }

    int n = L->row;
    int nrhs = B->column;

    for (int i = 0; i < n; ++i) {
        for (int col = 0; col < nrhs; ++col) {
            Y->data[MatrixIndex(Y, i, col)] = B->data[MatrixIndex(B, i, col)];
        }

        for (int j = 0; j < i; ++j) {
            for (int col = 0; col < nrhs; ++col) {
                Y->data[MatrixIndex(Y, i, col)] -= L->data[MatrixIndex(L, i, j)] * Y->data[MatrixIndex(Y, j, col)];
            }
        }

        REAL diag = L->data[MatrixIndex(L, i, i)];
        if (fabs(diag) < tol) {
            return MATRIX_ERROR_SINGULAR;
        }

        for (int col = 0; col < nrhs; ++col) {
            Y->data[MatrixIndex(Y, i, col)] /= diag;
        }
    }

    return MATRIX_SUCCESS;
}

/*
 * Back substitution:
 *
 *   U x = y.
 */
MatrixError BackSubstitution(const Matrix *U, const Matrix *y, Matrix *x, REAL tol)
{
    MatrixError error = CheckTriangularSolveInput(U, y, x);
    if (error != MATRIX_SUCCESS) {
        return error;
    }

    int n = U->row;

    for (int i = n - 1; i >= 0; --i) {
        REAL sum = y->data[i];

        for (int j = i + 1; j < n; ++j) {
            sum -= U->data[MatrixIndex(U, i, j)] * x->data[j];
        }

        REAL diag = U->data[MatrixIndex(U, i, i)];
        if (fabs(diag) < tol) {
            return MATRIX_ERROR_SINGULAR;
        }

        x->data[i] = sum / diag;
    }

    return MATRIX_SUCCESS;
}

/*
 * Back substitution with multiple right-hand sides:
 *
 *   U X = Y.
 */
MatrixError BackSubstitutionMultiple(const Matrix *U, const Matrix *Y, Matrix *X, REAL tol)
{
    MatrixError error = CheckTriangularSolveMultipleInput(U, Y, X);
    if (error != MATRIX_SUCCESS) {
        return error;
    }

    int n = U->row;
    int nrhs = Y->column;

    for (int i = n - 1; i >= 0; --i) {
        for (int col = 0; col < nrhs; ++col) {
            X->data[MatrixIndex(X, i, col)] = Y->data[MatrixIndex(Y, i, col)];
        }

        for (int j = i + 1; j < n; ++j) {
            for (int col = 0; col < nrhs; ++col) {
                X->data[MatrixIndex(X, i, col)] -= U->data[MatrixIndex(U, i, j)] * X->data[MatrixIndex(X, j, col)];
            }
        }

        REAL diag = U->data[MatrixIndex(U, i, i)];
        if (fabs(diag) < tol) {
            return MATRIX_ERROR_SINGULAR;
        }

        for (int col = 0; col < nrhs; ++col) {
            X->data[MatrixIndex(X, i, col)] /= diag;
        }
    }

    return MATRIX_SUCCESS;
}

/*
 * Determinant from U in no-pivot LU:
 *
 *   det(A) = prod_i U(i,i).
 */
MatrixError LUDeterminant(const Matrix *U, REAL *det)
{
    if (det == NULL) {
        return MATRIX_ERROR_NULL_POINTER;
    }

    MatrixError error = CheckSquareMatrix(U);
    if (error != MATRIX_SUCCESS) {
        return error;
    }

    REAL value = 1.0;

    for (int i = 0; i < U->row; ++i) {
        value *= U->data[MatrixIndex(U, i, i)];
    }

    *det = value;
    return MATRIX_SUCCESS;
}
