#ifndef MATRIX_RAND_H
#define MATRIX_RAND_H

#include "matrix_core.h"

MatrixError MatrixFillRandom(Matrix *A, REAL min, REAL max);
MatrixError MatrixGenerateRandom(Matrix *A, int row, int column, REAL min, REAL max);

#endif
