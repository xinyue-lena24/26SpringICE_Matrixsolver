#include "matrix_ops.h"

#include <math.h>

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

MatrixError MatrixTranspose(const Matrix *A, Matrix *AT)
{
    if (!MatrixIsValid(A) || !MatrixIsValid(AT)) {
        return MATRIX_ERROR_NULL_POINTER;
    }
    if (AT->row != A->column || AT->column != A->row) {
        return MATRIX_ERROR_SIZE_MISMATCH;
    }
    for (int i = 0; i < A->row; ++i) {
        for (int j = 0; j < A->column; ++j) {
            AT->data[MatrixIndex(AT, j, i)] = A->data[MatrixIndex(A, i, j)];
        }
    }
    return MATRIX_SUCCESS;
}

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
            for (int k = 0; k < A->column; ++k) {
                sum += A->data[MatrixIndex(A, i, k)] * BT.data[MatrixIndex(&BT, j, k)];
            }
            C->data[MatrixIndex(C, i, j)] = sum;
        }
    }
    MatrixFree(&BT);
    return MATRIX_SUCCESS;
}