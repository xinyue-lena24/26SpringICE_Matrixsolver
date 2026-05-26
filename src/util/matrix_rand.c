#include "matrix_rand.h"
#include <stdlib.h>

/**
 * @brief Fill all elements of a matrix A with random values uniformly distributed between min and max.
 *       This function checks if A is a valid matrix before filling it with random values.
 * @param A Pointer to the matrix to be filled with random values.
 * @param min The minimum value for the random numbers (inclusive).
 * @param max The maximum value for the random numbers (inclusive).
 * @return MatrixError code indicating success or the type of error encountered.
 */
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

/**
 * @brief Create a matrix A with the specified number of rows and columns, and fill it with random values uniformly distributed between min and max.
 *       This function first creates a matrix A with the specified dimensions and then fills it with random values. It checks for errors during matrix creation and filling, returning the appropriate MatrixError code if any issues are encountered.
 * @param A Pointer to the matrix to be created and filled with random values.
 * @param row The number of rows for the matrix A.
 * @param column The number of columns for the matrix A.
 * @param min The minimum value for the random numbers (inclusive).
 * @param max The maximum value for the random numbers (inclusive).
 * @return MatrixError code indicating success or the type of error encountered.
 */
MatrixError MatrixGenerateRandom(Matrix *A, int row, int column, REAL min, REAL max)
{
    MatrixError error = MatrixCreate(A, row, column);
    if (error != MATRIX_SUCCESS) {
        return error;
    }
    return MatrixFillRandom(A, min, max);
}