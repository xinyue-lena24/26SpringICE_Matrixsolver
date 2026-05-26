#include "matrix_ops.h"

#include <math.h>

/**
 * @brief Check if the input matrices A, B, C are valid for binary elementwise operations (like addition and subtraction).
 *        This function checks if A, B, C are valid matrices, if A and B have the same shape, and if C has the same shape as A and B.
 * @param A Pointer to the first input matrix.
 * @param B Pointer to the second input matrix.
 * @param C Pointer to the output matrix.
 * @return MatrixError code indicating success or the type of error encountered.
 */
static MatrixError CheckBinaryElementwise(const Matrix *A, const Matrix *B, Matrix *C)
{
    if (!MatrixIsValid(A) || !MatrixIsValid(B) || !MatrixIsValid(C)) {
        return MATRIX_ERROR_NULL_POINTER;
    }
    if (A->row != B->row || A->column != B->column) {
        return MATRIX_ERROR_SIZE_MISMATCH;
    }
    if (C->row != A->row || C->column != A->column) {
        return MATRIX_ERROR_SIZE_MISMATCH;
    }
    return MATRIX_SUCCESS;
}

/**
 * @brief Perform elementwise addition of two matrices A and B, storing the result in C.
 *        This function first checks if the input matrices are valid for addition using CheckBinaryElementwise and then performs the addition if the check passes.
 * @param A Pointer to the first input matrix.
 * @param B Pointer to the second input matrix.
 * @param C Pointer to the output matrix where the result will be stored.
 * @return MatrixError code indicating success or the type of error encountered.
 */
MatrixError MatrixAdd(const Matrix *A, const Matrix *B, Matrix *C)
{
    MatrixError error = CheckBinaryElementwise(A, B, C);
    if (error != MATRIX_SUCCESS) {
        return error;
    }
    int total = A->row * A->column;
    for (int k = 0; k < total; ++k) {
        C->data[k] = A->data[k] + B->data[k];
    }
    return MATRIX_SUCCESS;
}

/**
 * @brief Perform elementwise subtraction of two matrices A and B, storing the result in C.
 *        This function first checks if the input matrices are valid for subtraction using CheckBinaryElementwise and then performs the subtraction if the check passes.
 * @param A Pointer to the first input matrix.
 * @param B Pointer to the second input matrix.
 * @param C Pointer to the output matrix where the result will be stored.
 * @return MatrixError code indicating success or the type of error encountered.
 */
MatrixError MatrixSub(const Matrix *A, const Matrix *B, Matrix *C)
{
    MatrixError error = CheckBinaryElementwise(A, B, C);
    if (error != MATRIX_SUCCESS) {
        return error;
    }
    int total = A->row * A->column;
    for (int k = 0; k < total; ++k) {
        C->data[k] = A->data[k] - B->data[k];
    }
    return MATRIX_SUCCESS;
}

/**
 * @brief Scale a matrix A by a scalar alpha, storing the result in B.
 *       This function checks if A and B are valid matrices and if they have the same shape before performing the scaling operation.
 * @param alpha The scalar value to scale the matrix A.
 * @param A Pointer to the input matrix to be scaled.
 * @param B Pointer to the output matrix where the scaled result will be stored.
 * @return MatrixError code indicating success or the type of error encountered.
 */
MatrixError MatrixScale(REAL alpha, const Matrix *A, Matrix *B)
{
    if (!MatrixIsValid(A) || !MatrixIsValid(B)) {
        return MATRIX_ERROR_NULL_POINTER;
    }
    if (A->row != B->row || A->column != B->column) {
        return MATRIX_ERROR_SIZE_MISMATCH;
    }
    int total = A->row * A->column;
    for (int k = 0; k < total; ++k) {
        B->data[k] = alpha * A->data[k];
    }
    return MATRIX_SUCCESS;
}

/**
 * @brief Transpose a matrix A, storing the result in AT.
 *      This function checks if A and AT are valid matrices and if AT has the correct shape (columns of A become rows of AT and vice versa) before performing the transposition.
 * @param A Pointer to the input matrix to be transposed.
 * @param AT Pointer to the output matrix where the transposed result will be stored.
 * @return MatrixError code indicating success or the type of error encountered.
 */
MatrixError MatrixTranspose(const Matrix *A, Matrix *AT)
{
    if (!MatrixIsValid(A) || !MatrixIsValid(AT)) {
        return MATRIX_ERROR_NULL_POINTER;
    }
    if (AT->row != A->column || AT->column != A->row) {
        return MATRIX_ERROR_SIZE_MISMATCH;
    }
    for (int i = 0; i < A->row; ++i) {
        REAL* a_row = A->data + MatrixIndex(A, i, 0);
        REAL* at_row = AT->data + MatrixIndex(AT, 0, i);
        for (int j = 0; j < A->column; ++j) {
            *at_row = a_row[j];
            at_row += AT->row;
        }
    }
    return MATRIX_SUCCESS;
}

/**
 * @brief Compute the Frobenius norm of a matrix A, storing the result in norm_value.
 *     This function checks if A is a valid matrix and if norm_value is a valid pointer before performing the computation of the Frobenius norm, which is the square root of the sum of the squares of all elements in A.
 * @param A Pointer to the input matrix for which the Frobenius norm will be computed.
 * @param norm_value Pointer to a REAL variable where the computed Frobenius norm will be stored.
 * @return MatrixError code indicating success or the type of error encountered.
 */
MatrixError MatrixNormFrobenius(const Matrix *A, REAL *norm_value)
{
    if (!MatrixIsValid(A) || norm_value == NULL) {
        return MATRIX_ERROR_NULL_POINTER;
    }
    REAL sum = 0.0;
    int total = A->row * A->column;
    for (int k = 0; k < total; ++k) {
        sum += A->data[k] * A->data[k];
    }
    *norm_value = sqrt(sum);
    return MATRIX_SUCCESS;
}

/**
 * @brief Perform matrix multiplication of A and B, storing the result in C.
 *    This function checks if A, B, and C are valid matrices and if their dimensions are compatible for multiplication (columns of A must equal rows of B, and C must have the same number of rows as A and the same number of columns as B) before performing the multiplication using a standard triple-nested loop algorithm.
 * @param A Pointer to the first input matrix.
 * @param B Pointer to the second input matrix.
 * @param C Pointer to the output matrix where the result of A*B will be stored.
 * @return MatrixError code indicating success or the type of error encountered.
 */
MatrixError MatrixMultiply(const Matrix *A, const Matrix *B, Matrix *C)
{
    if (!MatrixIsValid(A) || !MatrixIsValid(B) || !MatrixIsValid(C)) {
        return MATRIX_ERROR_NULL_POINTER;
    }
    if (A->column != B->row || C->row != A->row || C->column != B->column) {
        return MATRIX_ERROR_SIZE_MISMATCH;
    }
    MatrixFillZero(C);
    
    Matrix BT;
    MatrixInit(&BT);
    MatrixCreate(&BT, B->column, B->row);
    MatrixTranspose(B, &BT);
    for (int i = 0; i < C->row; ++i) {
        for (int j = 0; j < C->column; ++j) {
            REAL sum = 0.0;
            REAL* a_idx = A->data + MatrixIndex(A, i, 0);
            REAL* b_idx = BT.data + MatrixIndex(&BT, j, 0);
            for (int k = 0; k < A->column; ++k) {
                sum += a_idx[k] * b_idx[k];
            }
            C->data[MatrixIndex(C, i, j)] = sum;
        }
    }
    MatrixFree(&BT);
    return MATRIX_SUCCESS;
}