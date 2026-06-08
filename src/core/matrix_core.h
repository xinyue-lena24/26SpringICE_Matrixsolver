/******************************************************************************
 * File        : matrix_core.h
 * Course      : C Program and Algorithm Analysis Lab
 * Lab         : Lab 4, Matrix Multiplication
 * Author      : Gong Helin
 * Year        : 2026
 * Description : Basic matrix data structure and core utility interfaces.
 ******************************************************************************/

#ifndef MATRIX_CORE_H
#define MATRIX_CORE_H

#include <stddef.h>

/* 本课程统一使用 REAL 表示矩阵元素类型。后续可改为 float 或其他数值类型。 */
typedef double REAL;

/*
 * MatrixError 用于统一表示矩阵库函数的返回状态。
 * 约定：MATRIX_SUCCESS 表示成功，其他值表示具体错误类型。
 */
typedef enum {
    MATRIX_SUCCESS = 0,
    MATRIX_ERROR_NULL_POINTER = 1,
    MATRIX_ERROR_INVALID_SIZE = 2,
    MATRIX_ERROR_SIZE_OVERFLOW = 3,
    MATRIX_ERROR_ALLOC_FAILED = 4,
    MATRIX_ERROR_INDEX_OUT_OF_RANGE = 5,
    MATRIX_ERROR_SIZE_MISMATCH = 6,
    MATRIX_ERROR_NOT_SQUARE = 7,
    MATRIX_ERROR_SINGULAR = 8,
    MATRIX_ERROR_ALREADY_ALLOCATED = 9
} MatrixError;

/*
 * Matrix 结构体由两部分组成：
 * 1. row 和 column 记录矩阵尺寸；
 * 2. data 指向一段连续的一维数组，按 row-major 方式存储矩阵元素。
 */
typedef struct {
    int row;
    int column;
    REAL *data;
} Matrix;

void MatrixInit(Matrix *A);
MatrixError MatrixCreate(Matrix *A, int row, int column);
void MatrixFree(Matrix *A);

int MatrixIsValid(const Matrix *A);
int MatrixHasShape(const Matrix *A, int row, int column);
int MatrixIndex(const Matrix *A, int i, int j);

MatrixError MatrixSet(Matrix *A, int i, int j, REAL value);
MatrixError MatrixGet(const Matrix *A, int i, int j, REAL *value);
MatrixError MatrixFillZero(Matrix *A);
MatrixError MatrixFillConstant(Matrix *A, REAL value);              //新增
MatrixError MatrixFillSequence(Matrix *A, REAL start, REAL step);
MatrixError MatrixFillIdentity(Matrix *A);                          //新增
MatrixError MatrixCopy(const Matrix *src, Matrix *dst);

void MatrixPrint(const Matrix *A, const char *name);
void MatrixPrintIndexMap(const Matrix *A, const char *name);
const char *MatrixErrorMessage(MatrixError error);

#endif
