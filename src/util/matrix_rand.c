#include "matrix_rand.h"
#include <stdlib.h>

MatrixError MatrixFillRandom(Matrix *A, REAL min, REAL max)
{
    if (!MatrixIsValid(A)) {
        return MATRIX_ERROR_NULL_POINTER;
    }
    for (int k = 0; k < A->row * A->column; ++k) {
        REAL r = (REAL)rand() / (REAL)RAND_MAX;
        A->data[k] = min + r * (max - min);
    }
    return MATRIX_SUCCESS;
}

MatrixError MatrixGenerateRandom(Matrix *A, int row, int column, REAL min, REAL max)
{
    MatrixError error = MatrixCreate(A, row, column);
    if (error != MATRIX_SUCCESS) {
        return error;
    }
    return MatrixFillRandom(A, min, max);
}