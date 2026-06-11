/******************************************************************************
 * File        : main_basic.c
 * Project     : Small C Matrix Computation Library
 * Description : Basic correctness test for the matrix library.
 *
 * This program tests the four required operations in the order specified
 * by the assignment:
 *
 *   1. Matrix addition
 *   2. Matrix scalar multiplication
 *   3. Matrix transpose
 *   4. Matrix multiplication
 *
 * The test matrices are intentionally small, so that the expected results
 * can be checked by hand. This program is mainly used to demonstrate that
 * the matrix functions are mathematically correct and that the basic workflow
 * of the library is clear.
 ******************************************************************************/

#include "matrix_core.h"
#include "matrix_ops.h"
#include "matrix_lu.h"
#include "matrix_solve.h"

#include <stdio.h>

/*
 * CheckError
 * ----------
 * This helper function checks whether a matrix library function has returned
 * MATRIX_SUCCESS.
 *
 * If an error occurs, the function prints the name of the failed operation
 * together with the corresponding error message.
 *
 * Return value:
 *   1  if the operation is successful;
 *   0  otherwise.
 *
 * `static` makes a function or variable visible only inside the current `.c` file.
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
 * TestMatrixAdd
 * -------------
 * Test matrix addition with a small 2x2 example.
 *
 * A = [1  2]
 *     [3  4]
 *
 * B = [5  6]
 *     [7  8]
 *
 * Expected result:
 *
 * C = A + B = [ 6   8]
 *             [10  12]
 */
static int TestMatrixAdd(void)
{
    printf("\n========== 1. Matrix Addition ==========\n");

    Matrix A, B, C;
    MatrixInit(&A);
    MatrixInit(&B);
    MatrixInit(&C);

    if (!CheckError(MatrixCreate(&A, 2, 2), "MatrixCreate(A)") ||
        !CheckError(MatrixCreate(&B, 2, 2), "MatrixCreate(B)") ||
        !CheckError(MatrixCreate(&C, 2, 2), "MatrixCreate(C)")) {
        MatrixFree(&A);
        MatrixFree(&B);
        MatrixFree(&C);
        return 0;
    }

    /*
     * Fill A and B by row-major order.
     *
     * MatrixFillSequence(&A, 1.0, 1.0) gives:
     * A->data = {1, 2, 3, 4}
     *
     * Since the matrix is stored in row-major order, this corresponds to:
     * A = [1 2; 3 4].
     */
    MatrixFillSequence(&A, 1.0, 1.0);
    MatrixFillSequence(&B, 5.0, 1.0);

    if (!CheckError(MatrixAdd(&A, &B, &C), "MatrixAdd")) {
        MatrixFree(&A);
        MatrixFree(&B);
        MatrixFree(&C);
        return 0;
    }

    MatrixPrint(&A, "A");
    MatrixPrint(&B, "B");
    MatrixPrint(&C, "C = A + B");
    printf("Expected: C = [6 8; 10 12]\n");

    MatrixFree(&A);
    MatrixFree(&B);
    MatrixFree(&C);
    return 1;
}

/*
 * TestMatrixScale
 * ---------------
 * Test matrix scalar multiplication with a 2x3 example.
 *
 * A = [1  2  3]
 *     [4  5  6]
 *
 * alpha = 2
 *
 * Expected result:
 *
 * B = alpha * A = [2   4   6]
 *                 [8  10  12]
 */
static int TestMatrixScale(void)
{
    printf("\n========== 2. Matrix Scalar Multiplication ==========\n");

    Matrix A, B;
    MatrixInit(&A);
    MatrixInit(&B);

    if (!CheckError(MatrixCreate(&A, 2, 3), "MatrixCreate(A)") ||
        !CheckError(MatrixCreate(&B, 2, 3), "MatrixCreate(B)")) {
        MatrixFree(&A);
        MatrixFree(&B);
        return 0;
    }

    MatrixFillSequence(&A, 1.0, 1.0);

    REAL alpha = 2.0;

    if (!CheckError(MatrixScale(alpha, &A, &B), "MatrixScale")) {
        MatrixFree(&A);
        MatrixFree(&B);
        return 0;
    }

    MatrixPrint(&A, "A");
    printf("alpha = %.4f\n", alpha);
    MatrixPrint(&B, "B = alpha * A");
    printf("Expected: B = [2 4 6; 8 10 12]\n");

    MatrixFree(&A);
    MatrixFree(&B);
    return 1;
}

/*
 * TestMatrixTranspose
 * -------------------
 * Test matrix transpose with a 2x3 matrix.
 *
 * A = [1  2  3]
 *     [4  5  6]
 *
 * Expected result:
 *
 * A^T = [1  4]
 *       [2  5]
 *       [3  6]
 *
 * Notice that if A is 2x3, then A^T must be 3x2.
 */
static int TestMatrixTranspose(void)
{
    printf("\n========== 3. Matrix Transpose ==========\n");

    Matrix A, AT;
    MatrixInit(&A);
    MatrixInit(&AT);

    if (!CheckError(MatrixCreate(&A, 2, 3), "MatrixCreate(A)") ||
        !CheckError(MatrixCreate(&AT, 3, 2), "MatrixCreate(AT)")) {
        MatrixFree(&A);
        MatrixFree(&AT);
        return 0;
    }

    MatrixFillSequence(&A, 1.0, 1.0);

    if (!CheckError(MatrixTranspose(&A, &AT), "MatrixTranspose")) {
        MatrixFree(&A);
        MatrixFree(&AT);
        return 0;
    }

    MatrixPrint(&A, "A");
    MatrixPrint(&AT, "AT = transpose(A)");
    printf("Expected: AT = [1 4; 2 5; 3 6]\n");

    MatrixFree(&A);
    MatrixFree(&AT);
    return 1;
}

/*
 * TestMatrixMultiply
 * ------------------
 * Test matrix multiplication with a 2x3 matrix and a 3x2 matrix.
 *
 * A = [1  2  3]
 *     [4  5  6]
 *
 * B = [ 7   8]
 *     [ 9  10]
 *     [11  12]
 *
 * Since A is 2x3 and B is 3x2, the product C = AB is 2x2.
 *
 * Expected result:
 *
 * C(0,0) = 1*7 + 2*9 + 3*11 = 58
 * C(0,1) = 1*8 + 2*10 + 3*12 = 64
 * C(1,0) = 4*7 + 5*9 + 6*11 = 139
 * C(1,1) = 4*8 + 5*10 + 6*12 = 154
 *
 * Therefore:
 *
 * C = [ 58   64]
 *     [139  154]
 */
static int TestMatrixMultiply(void)
{
    printf("\n========== 4. Matrix Multiplication ==========\n");

    Matrix A, B, C;
    MatrixInit(&A);
    MatrixInit(&B);
    MatrixInit(&C);

    if (!CheckError(MatrixCreate(&A, 2, 3), "MatrixCreate(A)") ||
        !CheckError(MatrixCreate(&B, 3, 2), "MatrixCreate(B)") ||
        !CheckError(MatrixCreate(&C, 2, 2), "MatrixCreate(C)")) {
        MatrixFree(&A);
        MatrixFree(&B);
        MatrixFree(&C);
        return 0;
    }

    MatrixFillSequence(&A, 1.0, 1.0);

    /*
     * Manually fill B so that the example is easy to check by hand.
     */
    MatrixSet(&B, 0, 0, 7.0);
    MatrixSet(&B, 0, 1, 8.0);
    MatrixSet(&B, 1, 0, 9.0);
    MatrixSet(&B, 1, 1, 10.0);
    MatrixSet(&B, 2, 0, 11.0);
    MatrixSet(&B, 2, 1, 12.0);

    if (!CheckError(MatrixMultiply(&A, &B, &C), "MatrixMultiply")) {
        MatrixFree(&A);
        MatrixFree(&B);
        MatrixFree(&C);
        return 0;
    }

    MatrixPrint(&A, "A");
    MatrixPrint(&B, "B");
    MatrixPrint(&C, "C = A * B");
    printf("Expected: C = [58 64; 139 154]\n");

    MatrixFree(&A);
    MatrixFree(&B);
    MatrixFree(&C);
    return 1;
}

/*
 * Test Solving Linear Systems
 * ------------------
 * This part is not required by the assignment, but it is a good demonstration
 * of how to use the matrix library for solving linear systems.
 *
 * The test case is:
 *
 * A = [4  -2  1]
 *     [3   6  2]
 *     [1   0  5]
 *
 * b = [11; 2; 7]
 *
 * The expected solution is:
 *
 * x = [2; -1; 1]
*/

static int TestSolveLinearSystem(void)
{
    printf("\n========== 5. Solving Linear Systems Example ==========\n");

    Matrix A, b, x;
    MatrixInit(&A);
    MatrixInit(&b);
    MatrixInit(&x);

    if (!CheckError(MatrixCreate(&A, 3, 3), "MatrixCreate(A)") ||
        !CheckError(MatrixCreate(&b, 3, 1), "MatrixCreate(b)") ||
        !CheckError(MatrixCreate(&x, 3, 1), "MatrixCreate(x)")) {
        MatrixFree(&A);
        MatrixFree(&b);
        MatrixFree(&x);
        return 0;
    }

    // Fill A
    MatrixSet(&A, 0, 0, 4.0); MatrixSet(&A, 0, 1, -2.0); MatrixSet(&A, 0, 2, 1.0);
    MatrixSet(&A, 1, 0, 3.0); MatrixSet(&A, 1, 1, 6.0); MatrixSet(&A, 1, 2, 2.0);
    MatrixSet(&A, 2, 0, 1.0); MatrixSet(&A, 2, 1, 0.0); MatrixSet(&A, 2, 2, 5.0);

    // Fill b
    MatrixSet(&b, 0, 0, 11.0);
    MatrixSet(&b, 1, 0, 2.0);
    MatrixSet(&b, 2, 0, 7.0);

    // Gaussian elimination
    if (!CheckError(GaussianSolveVector(&A, &b, &x, 1e-9), "GaussianSolveVector")) {
        MatrixFree(&A);
        MatrixFree(&b);
        MatrixFree(&x);
        return 0;
    }

    MatrixPrint(&A, "Coefficient matrix A");
    MatrixPrint(&b, "Right-hand side vector b");

    puts("Gaussian elimination with partial pivoting:");
    MatrixPrint(&x, "Solution vector x");
    printf("Expected: x = [2; -1; 1]\n");

    // LU decomposition and solve
    if (!CheckError(LUFactorSolvMatrix(&A, &b, &x, 1e-9), "LUFactorSolvMatrix")) {
        MatrixFree(&A);
        MatrixFree(&b);
        MatrixFree(&x);
        return 0;
    }

    MatrixPrint(&A, "Coefficient matrix A");
    MatrixPrint(&b, "Right-hand side vector b");

    puts("LU decomposition and solve:");
    MatrixPrint(&x, "Solution vector x");
    printf("Expected: x = [2; -1; 1]\n");

    MatrixFree(&b);
    MatrixFree(&x);

    // Test multiple right-hand sides
    Matrix B, X;
    MatrixInit(&B);
    MatrixInit(&X);
    if (!CheckError(MatrixCreate(&B, 3, 2), "MatrixCreate(B)") ||
        !CheckError(MatrixCreate(&X, 3, 2), "MatrixCreate(X)")) {
        MatrixFree(&A);
        MatrixFree(&B);
        MatrixFree(&X);
        return 0;
    }
    // Fill B with two right-hand side vectors
    MatrixSet(&B, 0, 0, 11.0); MatrixSet(&B, 0, 1, 11.0);
    MatrixSet(&B, 1, 0, 2.0); MatrixSet(&B, 1, 1, 2.0);
    MatrixSet(&B, 2, 0, 7.0); MatrixSet(&B, 2, 1, 7.0);
    MatrixPrint(&A, "Coefficient matrix A");
    MatrixPrint(&B, "Right-hand side matrix B");

    puts("Gaussian elimination with partial pivoting for multiple right-hand sides:");
    if (!CheckError(GaussianSolveMatrix(&A, &B, &X, 1e-9), "GaussianSolveMatrix with multiple RHS")) {
        MatrixFree(&A);
        MatrixFree(&B);
        MatrixFree(&X);
        return 0;
    }
    MatrixPrint(&X, "Solution matrix X");

    puts("LU decomposition and solve with multiple right-hand sides:");
    if (!CheckError(LUFactorSolvMatrix(&A, &B, &X, 1e-9), "LUFactorSolvMatrix with multiple RHS")) {
        MatrixFree(&A);
        MatrixFree(&B);
        MatrixFree(&X);
        return 0;
    }
    MatrixPrint(&X, "Solution matrix X");
    printf("Expected: X = [[2 2]; [-1 -1]; [1 1]]\n");
    return 1;
}

/*
 * TestErrorHandling
 * -----------------
 * This function demonstrates that the library can detect obvious dimension
 * errors instead of producing meaningless results.
 *
 * This part is not the main correctness test, but it is useful for the report
 * because the assignment requires reasonable handling of invalid input cases.
 */
static void TestErrorHandling(void)
{
    printf("\n========== 5. Error Handling Example ==========\n");

    Matrix A, B, C;
    MatrixInit(&A);
    MatrixInit(&B);
    MatrixInit(&C);

    MatrixCreate(&A, 2, 3);
    MatrixCreate(&B, 2, 3);
    MatrixCreate(&C, 3, 2);

    /*
     * A and B have the same size, but C is 3x2.
     * For matrix addition, the output C should also be 2x3.
     * Therefore this call should return MATRIX_ERROR_SIZE_MISMATCH.
     */
    MatrixError error = MatrixAdd(&A, &B, &C);
    printf("MatrixAdd with wrong output size: %s\n",
           MatrixErrorMessage(error));

    /*
     * For matrix multiplication A * B, the number of columns of A must be
     * equal to the number of rows of B.
     *
     * Here A is 2x3 and B is 2x3, so A.column = 3 but B.row = 2.
     * Therefore this call should also return MATRIX_ERROR_SIZE_MISMATCH.
     */
    error = MatrixMultiply(&A, &B, &C);
    printf("MatrixMultiply with A.column != B.row: %s\n",
           MatrixErrorMessage(error));

    MatrixFree(&A);
    MatrixFree(&B);
    MatrixFree(&C);
}

int main(void)
{
    printf("============================================================\n");
    printf("Basic Test for Small C Matrix Computation Library\n");
    printf("============================================================\n");

    if (!TestMatrixAdd()) return 1;
    if (!TestMatrixScale()) return 1;
    if (!TestMatrixTranspose()) return 1;
    if (!TestMatrixMultiply()) return 1;
    if (!TestSolveLinearSystem()) return 1;

    TestErrorHandling();

    printf("\nAll basic tests finished.\n");
    return 0;
}
