/******************************************************************************
 * File        : main_timing_ops.c
 * Project     : Small C Matrix Computation Library
 * Description : Timing comparison for four basic matrix operations.
 *
 * This program measures the average running time of:
 *
 *   1. Matrix addition
 *   2. Matrix scalar multiplication
 *   3. Matrix transpose
 *   4. Matrix multiplication
 *
 * Matrix sizes:
 *
 *   100 x 100
 *   1000 x 1000
 *   2000 x 2000
 *
 * Each operation is repeated several times, and the average running time is
 * reported. Matrix creation, memory allocation, initialization and printing
 * are not included in the timing region.
 ******************************************************************************/

#include "matrix_core.h"
#include "matrix_ops.h"

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

/*
 * g_sink is used to prevent the compiler from completely removing
 * computations during optimization. 
 * (As this sript is only used for test, the results matrices will not be used in the following procedure)
 *
 * The variable is declared volatile, so the compiler must assume that its
 * value may be used outside the current optimization context.
 */
static volatile REAL g_sink = 0.0;

/*
 * Convert the difference of two clock_t values into seconds.
 *
 * clock() measures CPU time used by the current process. For a simple course
 * assignment, this is acceptable. For more serious benchmarking, one may use
 * clock_gettime() on Linux or other high-resolution timers.
 */
static double ElapsedSeconds(clock_t start, clock_t end)
{
    return (double)(end - start) / (double)CLOCKS_PER_SEC;
}

/*
 * Print error information if a matrix operation fails.
 */
static int CheckError(MatrixError error, const char *func_name)
{
    if (error != MATRIX_SUCCESS) {
        printf("[ERROR] %s failed: %s\n", func_name, MatrixErrorMessage(error));
        return 0;
    }
    return 1;
}

/*
 * Create four n x n matrices:
 *
 *   A and B are input matrices.
 *   C is used as the output matrix for addition, scalar multiplication
 *     and matrix multiplication.
 *   T is used as the output matrix for transpose.
 *
 * For square matrices, the transpose of an n x n matrix is still n x n.
 */
static int CreateSquareMatrices(Matrix *A, Matrix *B, Matrix *C, Matrix *T, int n)
{
    MatrixInit(A);
    MatrixInit(B);
    MatrixInit(C);
    MatrixInit(T);

    if (!CheckError(MatrixCreate(A, n, n), "MatrixCreate(A)")) {
        return 0;
    }

    if (!CheckError(MatrixCreate(B, n, n), "MatrixCreate(B)")) {
        MatrixFree(A);
        return 0;
    }

    if (!CheckError(MatrixCreate(C, n, n), "MatrixCreate(C)")) {
        MatrixFree(A);
        MatrixFree(B);
        return 0;
    }

    if (!CheckError(MatrixCreate(T, n, n), "MatrixCreate(T)")) {
        MatrixFree(A);
        MatrixFree(B);
        MatrixFree(C);
        return 0;
    }

    return 1;
}

/*
 * Release all matrices created by CreateSquareMatrices.
 */
static void FreeSquareMatrices(Matrix *A, Matrix *B, Matrix *C, Matrix *T)
{
    MatrixFree(A);
    MatrixFree(B);
    MatrixFree(C);
    MatrixFree(T);
}

/*
 * TimeAdd
 * -------
 * Measure the average running time of matrix addition:
 *
 *   C = A + B
 *
 * The operation is repeated repeat times. The returned value is the average
 * time per call.
 */
static double TimeAdd(const Matrix *A, const Matrix *B, Matrix *C, int repeat)
{
    clock_t start = clock();

    for (int r = 0; r < repeat; ++r) {
        MatrixError error = MatrixAdd(A, B, C);
        if (error != MATRIX_SUCCESS) {
            printf("[ERROR] MatrixAdd failed: %s\n", MatrixErrorMessage(error));
            return -1.0;
        }
    }

    clock_t end = clock();

    /*
     * Read one element from the result so that the compiler cannot assume
     * the computed result is completely unused.
     */
    g_sink += C->data[0];

    return ElapsedSeconds(start, end) / (double)repeat;
}

/*
 * TimeScale
 * ---------
 * Measure the average running time of scalar multiplication:
 *
 *   B = 2.0 * A
 */
static double TimeScale(const Matrix *A, Matrix *B, int repeat)
{
    clock_t start = clock();

    for (int r = 0; r < repeat; ++r) {
        MatrixError error = MatrixScale(2.0, A, B);
        if (error != MATRIX_SUCCESS) {
            printf("[ERROR] MatrixScale failed: %s\n", MatrixErrorMessage(error));
            return -1.0;
        }
    }

    clock_t end = clock();

    g_sink += B->data[0];

    return ElapsedSeconds(start, end) / (double)repeat;
}

/*
 * TimeTranspose
 * -------------
 * Measure the average running time of matrix transpose:
 *
 *   AT = transpose(A)
 */
static double TimeTranspose(const Matrix *A, Matrix *AT, int repeat)
{
    clock_t start = clock();

    for (int r = 0; r < repeat; ++r) {
        MatrixError error = MatrixTranspose(A, AT);
        if (error != MATRIX_SUCCESS) {
            printf("[ERROR] MatrixTranspose failed: %s\n", MatrixErrorMessage(error));
            return -1.0;
        }
    }

    clock_t end = clock();

    g_sink += AT->data[0];

    return ElapsedSeconds(start, end) / (double)repeat;
}

/*
 * TimeMultiply
 * ------------
 * Measure the average running time of matrix multiplication:
 *
 *   C = A * B
 *
 * For n x n square matrices, the computational complexity is O(n^3).
 * Therefore this operation is expected to be much more expensive than
 * addition, scalar multiplication and transpose.
 */
static double TimeMultiply(const Matrix *A, const Matrix *B, Matrix *C, int repeat)
{
    clock_t start = clock();

    for (int r = 0; r < repeat; ++r) {
        MatrixError error = MatrixMultiply(A, B, C);
        if (error != MATRIX_SUCCESS) {
            printf("[ERROR] MatrixMultiply failed: %s\n", MatrixErrorMessage(error));
            return -1.0;
        }
    }

    clock_t end = clock();

    g_sink += C->data[0];

    return ElapsedSeconds(start, end) / (double)repeat;
}

int main(int argc, char *argv[])
{
    /*
     * Default repeat count.
     *
     * For the final report, repeat = 5 is recommended.
     * If 2000x2000 multiplication is too slow on your machine, you may use
     * repeat = 3, but the report should honestly state the actual repeat count.
     */
    int repeat = 5;

    /*
     * The repeat count can be changed from the command line.
     *
     * Example:
     *   ./main_timing_ops 3
     *   ./main_timing_ops 5
     */
    if (argc >= 2) {
        repeat = atoi(argv[1]);
        if (repeat <= 0) {
            repeat = 5;
        }
    }

    /*
     * Matrix sizes required by the assignment.
     */
    int sizes[] = {100, 500, 1000, 2000};
    int nsizes = (int)(sizeof(sizes) / sizeof(sizes[0]));

    printf("============================================================\n");
    printf("Timing Test for Four Basic Matrix Operations\n");
    printf("============================================================\n");
    printf("Repeat count: %d\n", repeat);
    printf("Timing excludes matrix creation, initialization and printing.\n\n");

    printf("%-12s %-14s %-14s %-14s %-14s\n",
           "Size", "Add(s)", "Scale(s)", "Transpose(s)", "Multiply(s)");
    printf("--------------------------------------------------------------------------\n");

    FILE *csv = fopen("results/ops_results.csv", "w");
    if (csv != NULL) {
        fprintf(csv, "size,repeat,add_s,scale_s,transpose_s,multiply_s\n");
    } else {
        printf("[WARNING] Could not open results/ops_results.csv for writing.\n");
    }

    for (int idx = 0; idx < nsizes; ++idx) {
        int n = sizes[idx];

        Matrix A, B, C, T;

        if (!CreateSquareMatrices(&A, &B, &C, &T, n)) {
            printf("[ERROR] Failed to create matrices of size %d x %d.\n", n, n);
            continue;
        }

        /*
         * Initialize matrices before timing.
         *
         * The initialization is outside the timing region because the report
         * focuses on the cost of the four matrix operations themselves.
         */
        MatrixFillSequence(&A, 1.0, 0.001);
        MatrixFillSequence(&B, 2.0, 0.001);
        MatrixFillZero(&C);
        MatrixFillZero(&T);

        double t_add = TimeAdd(&A, &B, &C, repeat);
        double t_scale = TimeScale(&A, &C, repeat);
        double t_transpose = TimeTranspose(&A, &T, repeat);
        double t_multiply = TimeMultiply(&A, &B, &C, repeat);

        printf("%-12d %-14.6f %-14.6f %-14.6f %-14.6f\n",
               n, t_add, t_scale, t_transpose, t_multiply);

        if (csv != NULL) {
            fprintf(csv, "%d,%d,%.10f,%.10f,%.10f,%.10f\n",
                    n, repeat, t_add, t_scale, t_transpose, t_multiply);
        }

        FreeSquareMatrices(&A, &B, &C, &T);
    }

    if (csv != NULL) {
        fclose(csv);
        printf("Results written to results/ops_results.csv\n");
    }

    printf("\ng_sink = %.6f\n", g_sink);
    printf("Timing test finished.\n");

    return 0;
}
