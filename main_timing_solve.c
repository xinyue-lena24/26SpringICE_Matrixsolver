#include "timer.h"
#include "matrix_core.h"
#include "matrix_ops.h"
#include "matrix_lu.h"
#include "matrix_solve.h"
#include "matrix_rand.h"

#include <stdio.h>
#include <stdlib.h>
// #include <cblas.h>

/*
 * Check and print matrix library errors.
 */
static int CheckError(MatrixError error, const char *func_name)
{
    if (error != MATRIX_SUCCESS) {
        printf("[ERROR] %s failed: %s\n", func_name, MatrixErrorMessage(error));
        return 0;
    }
    return 1;
}

typedef struct {
    int size_A;
    int col_B;
} Scale;
Scale scales[105];
int total_scales = 0;

/*
Create data scale for timing tests.
*/
void CreateScales() {
    for (int n = 10; n <= 1000; n*=10) {
        for (int m = 1; m <= 10000; m*=10) {
            scales[total_scales].size_A = n;
            scales[total_scales].col_B = m;
            total_scales++;
        }
    }
}

/*
@param solving_function: A pointer to the matrix solving function to be tested.
@param fun_name: A string representing the name of the solving function (for printing purposes).
@param scale: A pointer to the scale parameters for the test.
*/
void TestSolveFunction(MatrixError (*solving_function)(const Matrix*, const Matrix*, Matrix*, REAL), const char *fun_name, Scale* scale) {
    srand(0); // Set a fixed seed for reproducibility

    printf("\n========== Timing %s with A: %dx%d, B: %dx%d ==========\n", fun_name, scale->size_A, scale->size_A, scale->size_A, scale->col_B);

    Matrix A, B, X;
    MatrixInit(&A);
    MatrixInit(&B);
    MatrixInit(&X);

    if (!CheckError(MatrixCreate(&A, scale->size_A, scale->size_A), "MatrixCreate(A)") ||
        !CheckError(MatrixCreate(&B, scale->size_A, scale->col_B), "MatrixCreate(B)") ||
        !CheckError(MatrixCreate(&X, scale->size_A, scale->col_B), "MatrixCreate(X)")) {
        MatrixFree(&A);
        MatrixFree(&B);
        MatrixFree(&X);
        return;
    }

    if (!CheckError(MatrixFillRandom(&A, -10.0, 10.0), "MatrixFillRandom(A)") ||
        !CheckError(MatrixFillRandom(&B, -10.0, 10.0), "MatrixFillRandom(B)")) {
        MatrixFree(&A);
        MatrixFree(&B);
        MatrixFree(&X);
        return;
    }

    Timer timer;
    timer_start(&timer);
    
    int test_num = 100; // Number of times to repeat the test for averaging
    for (int i = 0; i < test_num; i++) {
        if (!CheckError(solving_function(&A, &B, &X, 1e-10), fun_name)) {
            MatrixFree(&A);
            MatrixFree(&B);
            MatrixFree(&X);
            return;
        }
    }

    double elapsed = timer_elapsed_ms(&timer);
    printf("Time taken: %.6f seconds\n", elapsed);

    MatrixFree(&A);
    MatrixFree(&B);
    MatrixFree(&X);
}

int main() {
    CreateScales();
    for (int i = 0; i < total_scales; i++) {
        TestSolveFunction(LUDecomposeSolveMultiple, "LUDecomposeSolveMultiple", &scales[i]);
    }
    return 0;
}