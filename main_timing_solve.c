#include "timer.h"
#include "matrix_core.h"
#include "matrix_ops.h"
#include "matrix_lu.h"
#include "matrix_solve.h"
#include "matrix_rand.h"

#include <stdio.h>
#include <stdlib.h>

// #define COMPARE_TO_BLAS

#ifdef COMPARE_TO_BLAS
#include <cblas.h>
#include <lapacke.h>
#endif

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

inline static long long total_ops_gauss(int n, int m) {
    return (long long)n * n * n * m; // Approximate total operations for solving
}
inline static long long total_ops_lu(int n, int m) {
    return (long long)n * n * n + 2 * (long long)n * n * m; // Approximate total operations for LU decomposition and solve
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
        for (int m = 1; m <= 1000; m*=10) {
            long long total_ops = total_ops_gauss(n, m);
            if (total_ops > 1e9) { // Skip scales that are too large
                continue;
            }
            scales[total_scales].size_A = n;
            scales[total_scales].col_B = m;
            total_scales++;
        }
    }
}

/*
* Tests a matrix solving function with specific scale parameters and logs the result to a CSV file.
* @param solving_function: A pointer to the matrix solving function to be tested.
* @param fun_name: A string representing the name of the solving function (for printing purposes).
* @param scale: A pointer to the scale parameters for the test.
* @param csv: A pointer to the CSV file for logging results.
*/
void TestSolveFunctionForScale(MatrixError (*solving_function)(const Matrix*, const Matrix*, Matrix*, REAL), const char *fun_name, Scale* scale, FILE *csv) {
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
    
    // Adjust the number of tests based on the size of the problem to keep total time reasonable
    int test_num = 1e8 / total_ops_gauss(scale->size_A, scale->col_B); // Aim for around 100 million operations total
    test_num = test_num < 5 ? 5 : test_num; // Ensure at least 5 test is run

    printf("Running %d tests...\n", test_num);

    for (int i = 0; i < test_num; i++) {
        if (!CheckError(solving_function(&A, &B, &X, 1e-10), fun_name)) {
            MatrixFree(&A);
            MatrixFree(&B);
            MatrixFree(&X);
            return;
        }
    }

    double elapsed_ms = timer_elapsed_ms(&timer);
    double avg_ms = elapsed_ms / (double)test_num;

    printf("Total time: %.6f ms\n", elapsed_ms);
    printf("Average time per solve: %.6f ms\n", avg_ms);
    if (csv != NULL) {
        fprintf(csv, "%d,%d,%.6f,%.6f,%d\n", scale->size_A, scale->col_B, elapsed_ms, avg_ms, test_num);
    }

    MatrixFree(&A);
    MatrixFree(&B);
    MatrixFree(&X);
}

/*
 * Tests a matrix solving function with different scale parameters and logs results to a CSV file.
 * @param solving_function: A pointer to the matrix solving function to be tested.
 * @param fun_name: A string representing the name of the solving function (for printing purposes).
 * @param csv_filename: The name of the CSV file to which results will be written.
 */
void TestSolveFunction(MatrixError (*solving_function)(const Matrix*, const Matrix*, Matrix*, REAL), const char *fun_name, const char *csv_filename) {
    FILE *csv = fopen(csv_filename, "w");
    if (csv == NULL) {
        fprintf(stderr, "Error: Could not open CSV file for writing.\n");
        return;
    }
    fprintf(csv, "Size_A,Col_B,Total_Time_ms,Avg_Time_ms,Test_Num\n");
    for (int i = 0; i < total_scales; i++) {
        TestSolveFunctionForScale(solving_function, fun_name, &scales[i], csv);
    }
    fclose(csv);
}

#ifdef COMPARE_TO_BLAS
/* Adjust blas for tests */
MatrixError BlasSolve(const Matrix *A, const Matrix *B, Matrix *X, REAL tol)
{
    if (!MatrixIsValid(A) || !MatrixIsValid(B) || !MatrixIsValid(X)) {
        return MATRIX_ERROR_NULL_POINTER;
    }
    if (A->row != A->column) {
        return MATRIX_ERROR_NOT_SQUARE;
    }
    if (B->row != A->row) {
        return MATRIX_ERROR_SIZE_MISMATCH;
    }
    if (X->row != A->row || X->column != B->column) {
        return MATRIX_ERROR_SIZE_MISMATCH;
    }

    int n = A->row;
    int nrhs = B->column;

    Matrix Acopy;
    MatrixInit(&Acopy);

    MatrixError err = MatrixCreate(&Acopy, n, n);
    if (err != MATRIX_SUCCESS) {
        return err;
    }

    MatrixCopy(A, &Acopy);
    MatrixCopy(B, X);

    lapack_int *ipiv =
        (lapack_int *)malloc(n * sizeof(lapack_int));

    if (ipiv == NULL) {
        MatrixFree(&Acopy);
        return MATRIX_ERROR_ALLOC_FAILED;
    }

    lapack_int info =
        LAPACKE_dgesv(
            LAPACK_ROW_MAJOR,
            n,
            nrhs,
            Acopy.data,
            n,
            ipiv,
            X->data,
            nrhs
        );

    free(ipiv);
    MatrixFree(&Acopy);
    
    if (info > 0) {
        return MATRIX_ERROR_SINGULAR;
    }

    return MATRIX_SUCCESS;
}
#endif

int main() {
    CreateScales();
    
    TestSolveFunction(LUFactorSolveMatrix, "LUFactorSolveMatrix", "results/lu_solve_matrix.csv");
    TestSolveFunction(GaussianSolveMatrix, "GaussianSolveMatrix", "results/gaussian_solve_matrix.csv");
    TestSolveFunction(LUFactorPivotSolveMatrix, "LUFactorPivotSolveMatrix", "results/lu_pivot_solve_matrix.csv");
    TestSolveFunction(GaussianSolveMatrixBatch, "GaussianSolveMatrixBatch", "results/gaussian_batch_solve_matrix.csv");
    
#ifdef COMPARE_TO_BLAS
    // test for blas
    TestSolveFunction(BlasSolve, "BlasSolve", "results/blas_solve.csv");
#endif
    return 0;
}