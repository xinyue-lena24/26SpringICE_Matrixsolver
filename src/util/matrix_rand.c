#include "matrix_rand.h"

#include <stdlib.h>

/*
 * Fill a valid matrix with random values uniformly distributed in [min, max].
 * The function does not call srand(); call srand(seed) once in main if needed.
 */
MatrixError MatrixFillRandom(Matrix *A, REAL min, REAL max)
{
    if (!MatrixIsValid(A)) {
        return MATRIX_ERROR_NULL_POINTER;
    }
    if (min > max) {
        return MATRIX_ERROR_INVALID_SIZE;
    }

    int total = A->row * A->column;
    for (int k = 0; k < total; ++k) {
        REAL r = (REAL)rand() / (REAL)RAND_MAX;
        A->data[k] = min + r * (max - min);
    }
    return MATRIX_SUCCESS;
}

/*
 * Create a matrix and fill it with random values in [min, max].
 */
MatrixError MatrixGenerateRandom(Matrix *A, int row, int column, REAL min, REAL max)
{
    MatrixError error = MatrixCreate(A, row, column);
    if (error != MATRIX_SUCCESS) {
        return error;
    }
    return MatrixFillRandom(A, min, max);
}
