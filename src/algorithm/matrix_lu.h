#ifndef MATRIX_LU_H
#define MATRIX_LU_H

#include "matrix_core.h"

MatrixError LUDecomposeNoPivot(const Matrix *A, Matrix *L, Matrix *U, REAL tol);
MatrixError ForwardSubstitution(const Matrix *L, const Matrix *b, Matrix *y, REAL tol);
MatrixError ForwardSubstitutionMultiple(const Matrix *L, const Matrix *B, Matrix *Y, REAL tol);
MatrixError BackSubstitution(const Matrix *U, const Matrix *y, Matrix *x, REAL tol);
MatrixError BackSubstitutionMultiple(const Matrix *U, const Matrix *Y, Matrix *X, REAL tol);
MatrixError LUDeterminant(const Matrix *U, REAL *det);

#endif
