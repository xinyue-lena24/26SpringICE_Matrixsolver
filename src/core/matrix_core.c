#include "matrix_core.h"

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

/**
 * @brief Initialize a matrix A.
 * @param A Pointer to the matrix to be initialized.
 */
inline void MatrixInit(Matrix *A)
{
    if (A == NULL) {
        return;
    }
    A->row = 0;
    A->column = 0;
    A->data = NULL;
}

/**
 * @brief Check if a matrix A is valid (non-NULL, has positive dimensions, and has allocated data).
 * @param A Pointer to the matrix to be checked.
 * @return Non-zero if the matrix is valid, zero otherwise.
 */
inline int MatrixIsValid(const Matrix *A)
{
    return (A != NULL && A->row > 0 && A->column > 0 && A->data != NULL);
}

/**
 * @brief Check if a matrix A has the specified shape (row x column).
 * @param A Pointer to the matrix to be checked.
 * @param row The expected number of rows.
 * @param column The expected number of columns.
 * @return Non-zero if the matrix has the specified shape, zero otherwise.
 */
inline int MatrixHasShape(const Matrix *A, int row, int column)
{
    return MatrixIsValid(A) && A->row == row && A->column == column;
}

/**
 * @brief Calculate the index in the data array for the element at row i and column j of matrix A.
 * @param A Pointer to the matrix.
 * @param i The row index (0-based).
 * @param j The column index (0-based).
 * @return The index in the data array corresponding to the element at (i, j).
 */
inline int MatrixIndex(const Matrix *A, int i, int j)
{
    return i * A->column + j;
}

/**
 * @brief Create a matrix A with the specified number of rows and columns, allocating memory for the data.
 *        This function checks if A is a valid pointer, if A already has allocated data (to prevent memory leaks), if the specified dimensions are positive, and if the total size does not cause overflow before allocating memory for the matrix data.
 * @param A Pointer to the matrix to be created.
 * @param row The number of rows for the matrix.
 * @param column The number of columns for the matrix.
 * @return MatrixError code indicating success or the type of error encountered.
 */
MatrixError MatrixCreate(Matrix *A, int row, int column)
{
    if (A == NULL) {
        return MATRIX_ERROR_NULL_POINTER;
    }
    if (A->data != NULL) {
        return MATRIX_ERROR_ALREADY_ALLOCATED;
    }
    if (row <= 0 || column <= 0) {
        return MATRIX_ERROR_INVALID_SIZE;
    }

    size_t nrow = (size_t)row;
    size_t ncol = (size_t)column;
    if (nrow > SIZE_MAX / ncol || nrow * ncol > SIZE_MAX / sizeof(REAL)) {
        return MATRIX_ERROR_SIZE_OVERFLOW;
    }

    A->data = (REAL *)malloc(nrow * ncol * sizeof(REAL));
    if (A->data == NULL) {
        MatrixInit(A);
        return MATRIX_ERROR_ALLOC_FAILED;
    }
    A->row = row;
    A->column = column;
    return MATRIX_SUCCESS;
}

/**
 * @brief Free the memory allocated for a matrix A.
 * @param A Pointer to the matrix whose memory is to be freed.
 */
void MatrixFree(Matrix *A)
{
    if (A == NULL) {
        return;
    }
    free(A->data);
    A->data = NULL;
    A->row = 0;
    A->column = 0;
}

/**
 * @brief Set the value of the element at row i and column j of matrix A to the specified value.
 *        This function checks if A is a valid matrix and if the specified indices are within the bounds of the matrix dimensions before setting the value.
 * @param A Pointer to the matrix.
 * @param i The row index (0-based).
 * @param j The column index (0-based).
 * @param value The value to be set at the specified position.
 * @return MatrixError code indicating success or the type of error encountered.
 */
MatrixError MatrixSet(Matrix *A, int i, int j, REAL value)
{
    if (!MatrixIsValid(A)) {
        return MATRIX_ERROR_NULL_POINTER;
    }
    if (i < 0 || i >= A->row || j < 0 || j >= A->column) {
        return MATRIX_ERROR_INDEX_OUT_OF_RANGE;
    }
    A->data[MatrixIndex(A, i, j)] = value;
    return MATRIX_SUCCESS;
}

/**
 * @brief Get the value of the element at row i and column j of matrix A, storing it in the variable pointed to by value.
 *        This function checks if A is a valid matrix, if value is a valid pointer, and if the specified indices are within the bounds of the matrix dimensions before retrieving the value.
 * @param A Pointer to the matrix.
 * @param i The row index (0-based).
 * @param j The column index (0-based).
 * @param value Pointer to a REAL variable where the retrieved value will be stored.
 * @return MatrixError code indicating success or the type of error encountered.
 */
MatrixError MatrixGet(const Matrix *A, int i, int j, REAL *value)
{
    if (!MatrixIsValid(A) || value == NULL) {
        return MATRIX_ERROR_NULL_POINTER;
    }
    if (i < 0 || i >= A->row || j < 0 || j >= A->column) {
        return MATRIX_ERROR_INDEX_OUT_OF_RANGE;
    }
    *value = A->data[MatrixIndex(A, i, j)];
    return MATRIX_SUCCESS;
}

/**
 * @brief Fill all elements of a matrix A with zeros.
 *        This function checks if A is a valid matrix before filling it with zeros.
 * @param A Pointer to the matrix to be filled with zeros.
 * @return MatrixError code indicating success or the type of error encountered.
 */
MatrixError MatrixFillZero(Matrix *A)
{
    if (!MatrixIsValid(A)) {
        return MATRIX_ERROR_NULL_POINTER;
    }
    int total = A->row * A->column;
    for (int k = 0; k < total; ++k) {
        A->data[k] = 0.0;
    }
    return MATRIX_SUCCESS;
}

/**
 * @brief Fill all elements of a matrix A with a specified value.
 *        This function checks if A is a valid matrix before filling it with the specified value.
 * @param A Pointer to the matrix to be filled.
 * @param value The value to fill the matrix with.
 * @return MatrixError code indicating success or the type of error encountered.
 */
MatrixError MatrixFillSequence(Matrix *A, REAL start, REAL step)
{
    if (!MatrixIsValid(A)) {
        return MATRIX_ERROR_NULL_POINTER;
    }
    int total = A->row * A->column;
    for (int k = 0; k < total; ++k) {
        A->data[k] = start + step * (REAL)k;
    }
    return MATRIX_SUCCESS;
}

/**
 * @brief Copy the contents of a source matrix src to a destination matrix dst.
 *        This function checks if src and dst are valid matrices and if they have the same shape before performing the copy operation.
 * @param src Pointer to the source matrix to be copied.
 * @param dst Pointer to the destination matrix where the contents will be copied.
 * @return MatrixError code indicating success or the type of error encountered.
 */
MatrixError MatrixCopy(const Matrix *src, Matrix *dst)
{
    if (!MatrixIsValid(src) || !MatrixIsValid(dst)) {
        return MATRIX_ERROR_NULL_POINTER;
    }
    if (src->row != dst->row || src->column != dst->column) {
        return MATRIX_ERROR_SIZE_MISMATCH;
    }
    int total = src->row * src->column;
    for (int k = 0; k < total; ++k) {
        dst->data[k] = src->data[k];
    }
    return MATRIX_SUCCESS;
}

/**
 * @brief Print the contents of a matrix A to the console, optionally with a name.
 *        This function checks if A is a valid matrix before printing its contents. If a name is provided, it will be printed as a header for the matrix.
 * @param A Pointer to the matrix to be printed.
 * @param name Optional name to be printed as a header for the matrix. If NULL, no name will be printed.
 */
void MatrixPrint(const Matrix *A, const char *name)
{
    if (name == NULL) {
        name = "Matrix";
    }
    if (!MatrixIsValid(A)) {
        printf("%s is an empty or invalid matrix.\n", name);
        return;
    }
    printf("%s = (%d x %d)\n", name, A->row, A->column);
    for (int i = 0; i < A->row; ++i) {
        printf("  ");
        for (int j = 0; j < A->column; ++j) {
            printf("%10.4f ", A->data[MatrixIndex(A, i, j)]);
        }
        printf("\n");
    }
}

/**
 * @brief Print the index mapping of a matrix A, showing how the 2D indices (i, j) map to the 1D data array index.
 *        This function checks if A is a valid matrix before printing the index mapping. If a name is provided, it will be used in the output to identify the matrix.
 * @param A Pointer to the matrix for which the index mapping will be printed.
 * @param name Optional name to be used in the output to identify the matrix. If NULL, "Matrix" will be used as the default name.
 */
void MatrixPrintIndexMap(const Matrix *A, const char *name)
{
    if (name == NULL) {
        name = "Matrix";
    }
    if (!MatrixIsValid(A)) {
        printf("%s is an empty or invalid matrix.\n", name);
        return;
    }
    printf("Index map for %s:\n", name);
    for (int i = 0; i < A->row; ++i) {
        for (int j = 0; j < A->column; ++j) {
            printf("  %s(%d,%d) -> data[%d]\n", name, i, j, MatrixIndex(A, i, j));
        }
    }
}

/**
 * @brief Get a human-readable error message corresponding to a MatrixError code.
 * @param error The MatrixError code for which the message is to be retrieved.
 * @return A string describing the error message corresponding to the provided MatrixError code.
 */
const char *MatrixErrorMessage(MatrixError error)
{
    switch (error) {
        case MATRIX_SUCCESS: return "success";
        case MATRIX_ERROR_NULL_POINTER: return "null pointer or invalid matrix";
        case MATRIX_ERROR_INVALID_SIZE: return "invalid matrix size";
        case MATRIX_ERROR_SIZE_OVERFLOW: return "matrix size overflow";
        case MATRIX_ERROR_ALLOC_FAILED: return "memory allocation failed";
        case MATRIX_ERROR_INDEX_OUT_OF_RANGE: return "index out of range";
        case MATRIX_ERROR_SIZE_MISMATCH: return "matrix size mismatch";
        case MATRIX_ERROR_NOT_SQUARE: return "matrix must be square";
        case MATRIX_ERROR_SINGULAR: return "matrix is singular or nearly singular";
        case MATRIX_ERROR_ALREADY_ALLOCATED: return "matrix data already allocated; free it first";
        default: return "unknown error";
    }
}
