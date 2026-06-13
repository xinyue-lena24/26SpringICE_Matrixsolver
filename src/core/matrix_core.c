/******************************************************************************
 * File        : matrix_core.c
 * Course      : C Program and Algorithm Analysis Lab
 * Lab         : Lab 4, Matrix Multiplication
 * Author      : Gong Helin
 * Year        : 2026
 * Description : Implementation of matrix creation, memory management,
 *               element access, initialization and printing utilities.
 ******************************************************************************/

#include "matrix_core.h"

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

/*
 * 将矩阵对象初始化为空矩阵。所有 Matrix 变量在使用前均应先调用该函数。
 * 这样可以避免未初始化指针导致的释放错误或重复分配错误。
 */
void MatrixInit(Matrix *A)
{
    if (A == NULL) {
        return;
    }
    A->row = 0;
    A->column = 0;
    A->data = NULL;
}

/* 判断矩阵对象是否处于可用状态。 */
int MatrixIsValid(const Matrix *A)
{
    return (A != NULL && A->row > 0 && A->column > 0 && A->data != NULL);
}

/* 判断矩阵是否具有指定尺寸。 */
int MatrixHasShape(const Matrix *A, int row, int column)
{
    return MatrixIsValid(A) && A->row == row && A->column == column;
}

/*
 * row-major 存储下的二维下标到一维下标转换。
 * 调用本函数前应保证 i 和 j 已经通过边界检查。
 */
inline int MatrixIndex(const Matrix *A, int i, int j)
{
    return i * A->column + j;
}

/*
 * 创建 row x column 矩阵，并分配连续内存。
 * 本函数要求 A 已经调用 MatrixInit，且当前没有持有已分配内存。
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

/* 释放矩阵持有的动态内存，并将矩阵恢复为空状态。 */
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

/* 设置矩阵第 (i,j) 个元素。下标采用 C 语言 0-based 规则。 */
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

/* 读取矩阵第 (i,j) 个元素。结果通过 value 指针返回。 */
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

/* 将矩阵所有元素置零。 */
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

/*
 * 将矩阵所有元素填充为同一个常数 value。
 *
 * 该函数是 MatrixFillZero 的一般化版本。
 * 当 value = 0.0 时，其作用与 MatrixFillZero 相同。
 */
MatrixError MatrixFillConstant(Matrix *A, REAL value)
{
    if (!MatrixIsValid(A)) {
        return MATRIX_ERROR_NULL_POINTER;
    }

    int total = A->row * A->column;
    for (int k = 0; k < total; ++k) {
        A->data[k] = value;
    }

    return MATRIX_SUCCESS;
}

/*
 * 将矩阵填充为单位矩阵。
 *
 * 要求：
 *   1. A 必须是有效矩阵；
 *   2. A 必须是方阵。
 *
 * 数学形式：
 *   A(i,j) = 1, if i == j;
 *   A(i,j) = 0, otherwise.
 */
MatrixError MatrixFillIdentity(Matrix *A)
{
    if (!MatrixIsValid(A)) {
        return MATRIX_ERROR_NULL_POINTER;
    }

    if (A->row != A->column) {
        return MATRIX_ERROR_NOT_SQUARE;
    }

    MatrixFillZero(A);

    for (int i = 0; i < A->row; ++i) {
        A->data[MatrixIndex(A, i, i)] = 1.0;
    }

    return MATRIX_SUCCESS;
}

/*
 * 用等差序列填充矩阵，便于课堂构造测试数据。
 * data[k] = start + step * k。
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

/* 将 src 的全部元素拷贝到同尺寸矩阵 dst。 */
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

/* 打印矩阵，主要用于小规模矩阵的课堂验证。 */
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

/* 打印二维下标到一维数组下标的映射关系。 */
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

/* 将错误码转换为便于输出的说明文字。 */
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
