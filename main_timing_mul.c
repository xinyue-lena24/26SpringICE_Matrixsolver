/******************************************************************************
 * File        : main_timing_mul.c
 * Project     : Small C Matrix Computation Library
 * Description : Timing comparison between two matrix multiplication versions.
 *
 * This program compares:
 *
 *   1. MatrixMultiply
 *      The basic i-j-k implementation in the source code.
 *
 *   2. MatrixMultiplyIKJ
 *      The i-k-j implementation. For row-major storage, this version often
 *      has better cache locality because the innermost loop traverses j,
 *      which corresponds to continuous memory access in a matrix row.
 *
 * The purpose of this program is to support the optional improvement section
 * of the report.
 ******************************************************************************/

#include "matrix_core.h"
#include "matrix_ops.h"

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

/*
 * Used to prevent the compiler from removing computations whose results
 * appear to be unused.
 */
static volatile REAL g_sink = 0.0;

/*
 * Convert CPU clock ticks to seconds.
 */
static double ElapsedSeconds(clock_t start, clock_t end)
{
    return (double)(end - start) / (double)CLOCKS_PER_SEC;
}

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

/*
 * Create four n x n matrices:
 *
 *   A and B are input matrices.
 *   C1 stores the result of MatrixMultiply.
 *   C2 stores the result of MatrixMultiplyIKJ.
 *
 * C1 and C2 are separated so that the two multiplication versions do not
 * overwrite each other's result. This also allows us to compare their results.
 */
static int CreateMatrices(Matrix *A, Matrix *B, Matrix *C1, Matrix *C2, int n)
{
    MatrixInit(A);
    MatrixInit(B);
    MatrixInit(C1);
    MatrixInit(C2);

    if (!CheckError(MatrixCreate(A, n, n), "MatrixCreate(A)")) {
        return 0;
    }

    if (!CheckError(MatrixCreate(B, n, n), "MatrixCreate(B)")) {
        MatrixFree(A);
        return 0;
    }

    if (!CheckError(MatrixCreate(C1, n, n), "MatrixCreate(C1)")) {
        MatrixFree(A);
        MatrixFree(B);
        return 0;
    }

    if (!CheckError(MatrixCreate(C2, n, n), "MatrixCreate(C2)")) {
        MatrixFree(A);
        MatrixFree(B);
        MatrixFree(C1);
        return 0;
    }

    return 1;
}

/*
 * Release matrices created by CreateMatrices.
 */
static void FreeMatrices(Matrix *A, Matrix *B, Matrix *C1, Matrix *C2)
{
    MatrixFree(A);
    MatrixFree(B);
    MatrixFree(C1);
    MatrixFree(C2);
}

/*
 * TimeMultiplyVersion
 * -------------------
 * Measure the average running time of a given matrix multiplication function.
 *
 * The first argument is a function pointer. This allows the same timing
 * function to test both MatrixMultiply and MatrixMultiplyIKJ.
 *
 * The expected function signature is:
 *
 *   MatrixError multiply(const Matrix *A, const Matrix *B, Matrix *C);
 */
static double TimeMultiplyVersion(
    MatrixError (*multiply)(const Matrix *, const Matrix *, Matrix *),
    const Matrix *A,
    const Matrix *B,
    Matrix *C,
    int repeat)
{
    clock_t start = clock();

    for (int r = 0; r < repeat; ++r) {
        MatrixError error = multiply(A, B, C);
        if (error != MATRIX_SUCCESS) {
            printf("[ERROR] matrix multiplication failed: %s\n",
                   MatrixErrorMessage(error));
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
     */
    int repeat = 5;

    /*
     * The repeat count can be specified from the command line.
     *
     * Example:
     *   ./main_timing_mul 3
     *   ./main_timing_mul 5
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
    printf("Timing Comparison of Matrix Multiplication Implementations\n");
    printf("============================================================\n");
    printf("Repeat count: %d\n", repeat);
    printf("Timing excludes matrix creation, initialization and printing.\n\n");

    printf("%-12s %-14s %-14s %-14s %-14s\n",
           "Size", "IJK(s)", "IKJ(s)", "Speedup", "MaxDiff");
    printf("--------------------------------------------------------------------------\n");

    FILE *csv = fopen("results/mul_results.csv", "w");
    if (csv != NULL) {
        fprintf(csv, "size,repeat,ijk_s,ikj_s,speedup,max_abs_diff\n");
    } else {
        printf("[WARNING] Could not open results/mul_results.csv for writing.\n");
    }

    for (int idx = 0; idx < nsizes; ++idx) {
        int n = sizes[idx];

        Matrix A, B, C_ijk, C_ikj;

        if (!CreateMatrices(&A, &B, &C_ijk, &C_ikj, n)) {
            printf("[ERROR] Failed to create matrices of size %d x %d.\n", n, n);
            continue;
        }

        /*
         * Fill input matrices outside the timing region.
         *
         * The same A and B are used for both multiplication implementations
         * to ensure a fair comparison.
         */
        MatrixFillSequence(&A, 1.0, 0.001);
        MatrixFillSequence(&B, 2.0, 0.001);
        MatrixFillZero(&C_ijk);
        MatrixFillZero(&C_ikj);

        /*
         * Measure the two multiplication functions.
         *
         * According to the current source code:
         *   MatrixMultiply    is the basic i-j-k version.
         *   MatrixMultiplyIKJ is the improved i-k-j version.
         */
        double t_ijk = TimeMultiplyVersion(MatrixMultiply, &A, &B, &C_ijk, repeat);
        double t_ikj = TimeMultiplyVersion(MatrixMultiplyIKJ, &A, &B, &C_ikj, repeat);

        /*
         * Check that both implementations produce the same result.
         */
        REAL max_diff = 0.0;
        MatrixMaxAbsDiff(&C_ijk, &C_ikj, &max_diff);

        /*
         * Speedup is defined as:
         *
         *   speedup = time of i-j-k version / time of i-k-j version
         *
         * If speedup > 1, then the i-k-j version is faster.
         */
        double speedup = 0.0;
        if (t_ikj > 0.0) {
            speedup = t_ijk / t_ikj;
        }

        printf("%-12d %-14.6f %-14.6f %-14.2f %-14.6e\n",
               n, t_ijk, t_ikj, speedup, max_diff);

        if (csv != NULL) {
            fprintf(csv, "%d,%d,%.10f,%.10f,%.10f,%.10e\n",
                    n, repeat, t_ijk, t_ikj, speedup, max_diff);
        }

        FreeMatrices(&A, &B, &C_ijk, &C_ikj);
    }

    if (csv != NULL) {
        fclose(csv);
        printf("Results written to results/mul_results.csv\n");
    }

    printf("\ng_sink = %.6f\n", g_sink);
    printf("Multiplication comparison finished.\n");

    return 0;
}
