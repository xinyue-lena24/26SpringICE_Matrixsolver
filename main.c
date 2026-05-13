#include "matrix_core.h"
#include "matrix_ops.h"
#include "matrix_rand.h"
#include "timer.h"

#include <stdio.h>

static void Check(MatrixError error, const char *name)
{
    printf("%-28s : %s\n", name, MatrixErrorMessage(error));
}

int main(void)
{
    struct Timer timer;
    Matrix A, B, C;

    // test
    MatrixInit(&A); MatrixInit(&B); MatrixInit(&C);

    srand(0); // 固定随机数种子，确保每次运行结果一致

    MatrixGenerateRandom(&A, 100, 100, 0.0, 1.0);
    MatrixGenerateRandom(&B, 100, 100000, 0.0, 1.0);
    MatrixError error = MatrixCreate(&C, 100, 100000);
    Check(error, "MatrixCreate C");

    printf("%d\n", MatrixIsValid(&C));

    MatrixFillSequence(&A, 1.0, 1.0);
    MatrixFillSequence(&B, 10.0, 5.0);

    Check(MatrixMultiply(&A, &B, &C), "MatrixMultiply");

    timer_start(&timer);
    int n = 1;
    for (int i = 0; i < n; ++i) {
        MatrixMultiply(&A, &B, &C);
    }
    double elapsed_ms = timer_elapsed_ms(&timer);
    printf("Time for %d multiplications: %.3f ms\n", n, elapsed_ms);

    MatrixFree(&A); MatrixFree(&B); MatrixFree(&C);
    return 0;
}
