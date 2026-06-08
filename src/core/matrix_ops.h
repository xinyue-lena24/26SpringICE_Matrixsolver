/******************************************************************************
 * File        : matrix_ops.h
 * Course      : C Program and Algorithm Analysis Lab
 * Lab         : Lab 4, Matrix Multiplication
 * Author      : Gong Helin
 * Year        : 2026
 * Description : Interfaces of elementary matrix operations used in Lab 4.
 ******************************************************************************/

#ifndef MATRIX_OPS_H
#define MATRIX_OPS_H

#include "matrix_core.h"

// 基础运算
MatrixError MatrixAdd(const Matrix *A, const Matrix *B, Matrix *C);
MatrixError MatrixSub(const Matrix *A, const Matrix *B, Matrix *C);
MatrixError MatrixScale(REAL alpha, const Matrix *A, Matrix *B);
MatrixError MatrixTranspose(const Matrix *A, Matrix *AT);

// 范数计算
MatrixError MatrixNormFrobenius(const Matrix *A, REAL *norm_value);
MatrixError MatrixNormOne(const Matrix *A, REAL *norm_value);                                   //新增
MatrixError MatrixNormInfinity(const Matrix *A, REAL *norm_value);                              //新增
MatrixError MatrixMaxAbsDiff(const Matrix *A, const Matrix *B, REAL *diff);                     //新增
MatrixError MatrixRelativeErrorFrobenius(const Matrix *A, const Matrix *B, REAL *rel_error);    //新增

// 乘法
MatrixError MatrixMultiply(const Matrix *A, const Matrix *B, Matrix *C);
MatrixError MatrixMultiplyIKJ(const Matrix *A, const Matrix *B, Matrix *C);

#endif
