/******************************************************************************
 * File        : matrix_ops.c
 * Course      : C Program and Algorithm Analysis Lab
 * Lab         : Lab 4, Matrix Multiplication
 * Author      : Gong Helin
 * Year        : 2026
 * Description : Implementation of elementary matrix operations and two matrix
 *               multiplication variants for teaching and comparison.
 ******************************************************************************/

#include "matrix_ops.h"

#include <math.h>

/*
 * 检查加法、减法等逐元素二元运算所需的尺寸条件。
 * A、B、C 三个矩阵必须均有效，且尺寸完全一致。
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

/*
 * 计算矩阵 1-范数。
 *
 * 数学定义：
 *   ||A||_1 = max_j sum_i |a_ij|
 *
 * 即：每一列取绝对值求和，然后取最大列和。
 */
MatrixError MatrixNormOne(const Matrix *A, REAL *norm_value)
{
    if (!MatrixIsValid(A) || norm_value == NULL) {
        return MATRIX_ERROR_NULL_POINTER;
    }

    REAL max_col_sum = 0.0;

    for (int j = 0; j < A->column; ++j) {
        REAL col_sum = 0.0;

        for (int i = 0; i < A->row; ++i) {
            col_sum += fabs(A->data[MatrixIndex(A, i, j)]);
        }

        if (col_sum > max_col_sum) {
            max_col_sum = col_sum;
        }
    }

    *norm_value = max_col_sum;
    return MATRIX_SUCCESS;
}

/*
 * 计算矩阵无穷范数。
 *
 * 数学定义：
 *   ||A||_infinity = max_i sum_j |a_ij|
 *
 * 即：每一行取绝对值求和，然后取最大行和。
 */
MatrixError MatrixNormInfinity(const Matrix *A, REAL *norm_value)
{
    if (!MatrixIsValid(A) || norm_value == NULL) {
        return MATRIX_ERROR_NULL_POINTER;
    }

    REAL max_row_sum = 0.0;

    for (int i = 0; i < A->row; ++i) {
        REAL row_sum = 0.0;

        for (int j = 0; j < A->column; ++j) {
            row_sum += fabs(A->data[MatrixIndex(A, i, j)]);
        }

        if (row_sum > max_row_sum) {
            max_row_sum = row_sum;
        }
    }

    *norm_value = max_row_sum;
    return MATRIX_SUCCESS;
}

/*
 * 计算两个同尺寸矩阵的最大元素绝对差。
 *
 * 数学定义：
 *   diff = max_ij |A(i,j) - B(i,j)|
 *
 * 该函数常用于比较两个算法输出是否一致。
 */
MatrixError MatrixMaxAbsDiff(const Matrix *A, const Matrix *B, REAL *diff)
{
    if (!MatrixIsValid(A) || !MatrixIsValid(B) || diff == NULL) {
        return MATRIX_ERROR_NULL_POINTER;
    }

    if (A->row != B->row || A->column != B->column) {
        return MATRIX_ERROR_SIZE_MISMATCH;
    }

    REAL max_diff = 0.0;
    int total = A->row * A->column;

    for (int k = 0; k < total; ++k) {
        REAL current_diff = fabs(A->data[k] - B->data[k]);

        if (current_diff > max_diff) {
            max_diff = current_diff;
        }
    }

    *diff = max_diff;
    return MATRIX_SUCCESS;
}

/*
 * 计算两个同尺寸矩阵之间的 Frobenius 相对误差。
 *
 * 数学定义：
 *   rel_error = ||A - B||_F / ||B||_F
 *
 * 这里 B 被视为参考矩阵。
 *
 * 特殊情况：
 *   如果 ||B||_F = 0 且 ||A-B||_F = 0，则相对误差定义为 0；
 *   如果 ||B||_F = 0 但 ||A-B||_F != 0，则返回 MATRIX_ERROR_SINGULAR。
 */
MatrixError MatrixRelativeErrorFrobenius(const Matrix *A,
                                         const Matrix *B,
                                         REAL *rel_error)
{
    if (!MatrixIsValid(A) || !MatrixIsValid(B) || rel_error == NULL) {
        return MATRIX_ERROR_NULL_POINTER;
    }

    if (A->row != B->row || A->column != B->column) {
        return MATRIX_ERROR_SIZE_MISMATCH;
    }

    REAL diff_sum = 0.0;
    REAL ref_sum = 0.0;
    int total = A->row * A->column;

    for (int k = 0; k < total; ++k) {
        REAL d = A->data[k] - B->data[k];
        diff_sum += d * d;
        ref_sum += B->data[k] * B->data[k];
    }

    REAL diff_norm = sqrt(diff_sum);
    REAL ref_norm = sqrt(ref_sum);
    if (ref_norm == 0.0) {
        if (diff_norm == 0.0) {
            *rel_error = 0.0;
            return MATRIX_SUCCESS;
        }
        return MATRIX_ERROR_SINGULAR;
    }

    *rel_error = diff_norm / ref_norm;
    return MATRIX_SUCCESS;
}

/*
 * 朴素 i-j-k 矩阵乘法。
 * 外两层循环遍历结果矩阵 C 的行和列，内层循环完成一次内积计算。
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
    for (int i = 0; i < A->row; ++i) {
        for (int j = 0; j < B->column; ++j) {
            REAL sum = 0.0;
            for (int k = 0; k < A->column; ++k) {
                sum += A->data[MatrixIndex(A, i, k)] *
                       B->data[MatrixIndex(B, k, j)];
            }
            C->data[MatrixIndex(C, i, j)] = sum;
        }
    }
    return MATRIX_SUCCESS;
}

/*
 * i-k-j 循环顺序矩阵乘法。
 * 该版本用于课堂比较：对固定的 i 和 k，连续遍历 j，C 的访问更连续。
 */
MatrixError MatrixMultiplyIKJ(const Matrix *A, const Matrix *B, Matrix *C)
{
    if (!MatrixIsValid(A) || !MatrixIsValid(B) || !MatrixIsValid(C)) {
        return MATRIX_ERROR_NULL_POINTER;
    }
    if (A->column != B->row || C->row != A->row || C->column != B->column) {
        return MATRIX_ERROR_SIZE_MISMATCH;
    }

    MatrixFillZero(C);
    for (int i = 0; i < A->row; ++i) {
        for (int k = 0; k < A->column; ++k) {
            REAL aik = A->data[MatrixIndex(A, i, k)];
            for (int j = 0; j < B->column; ++j) {
                C->data[MatrixIndex(C, i, j)] +=
                    aik * B->data[MatrixIndex(B, k, j)];
            }
        }
    }
    return MATRIX_SUCCESS;
}
