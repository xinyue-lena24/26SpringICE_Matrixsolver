#ifndef MATRIX_OPS_H
#define MATRIX_OPS_H

#include "matrix_core.h"

MatrixError MatrixAdd(const Matrix *A, const Matrix *B, Matrix *C);
MatrixError MatrixSub(const Matrix *A, const Matrix *B, Matrix *C);
MatrixError MatrixScale(REAL alpha, const Matrix *A, Matrix *B);
MatrixError MatrixTranspose(const Matrix *A, Matrix *AT);
MatrixError MatrixNormFrobenius(const Matrix *A, REAL *norm_value);
MatrixError MatrixMultiply(const Matrix *A, const Matrix *B, Matrix *C);

#endif
