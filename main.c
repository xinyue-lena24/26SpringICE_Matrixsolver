#include "matrix_core.h"
#include "matrix_ops.h"
#include "matrix_rand.h"
#include "timer.h"

#include <stdio.h>

static void Check(MatrixError error, const char *name)
{
    printf("%-28s : %s\n", name, MatrixErrorMessage(error));
}

typedef struct testScale {
    int n;
    int m;
    int p;
} testScale;

int test_num;
testScale test_scales[31*31*31];

void init_test_scales() {
    srand(0); // Reset random seed for reproducibility

    // test scales required
    test_scales[0] = (testScale){100, 100, 100};
    test_scales[1] = (testScale){1000, 1000, 1000};
    test_scales[2] = (testScale){2000, 2000, 2000};
    test_num = 3;

    // Generate test scales with n, m, p
    for (int n = 32; n <= 2048; n *= 4) {
        for (int m = 32; m <= 2048; m *= 4) {
            for (int p = 32; p <= 2048; p *= 4) {
                test_scales[test_num].n = n;
                test_scales[test_num].m = m;
                test_scales[test_num].p = p;
                test_num++;
            }
        }
    }
    // for (int i = 0; i < test_num; ++i) {
    //     test_scales[i].n = rand() % 1000 + 1;
    //     test_scales[i].m = rand() % 1000 + 1;
    //     test_scales[i].p = rand() % 1000 + 1;
    // }
}

void test_multiply()
{
    struct Timer timer;
    Matrix A, B, C;
    FILE *resultFile = fopen("multiply_results.csv", "w+");
    fprintf(resultFile, "n,m,p,elapsed_ms,test_num\n");

    srand(0); // Reset random seed for reproducibility
    
    for (int i = 0; i < test_num; ++i) {
        int n = test_scales[i].n;
        int m = test_scales[i].m;
        int p = test_scales[i].p;

        printf("Test %d: A(%d x %d) * B(%d x %d) = C(%d x %d)\n", i + 1, n, p, p, m, n, m);

        MatrixInit(&A); MatrixInit(&B); MatrixInit(&C);
        MatrixGenerateRandom(&A, n, p, 0.0, 1.0);
        MatrixGenerateRandom(&B, p, m, 0.0, 1.0);
        MatrixError error = MatrixCreate(&C, n, m);
        // Check(error, "MatrixCreate C");
        // printf("%d\n", MatrixIsValid(&C));

        error = MatrixMultiply(&A, &B, &C);
        // Check(error, "MatrixMultiply");

        int scale = n * m * p;
        printf("Scale (n*m*p): %d\n", scale);
        int test_num = (1e8 / scale) + 10; // Adjust test_num based on the scale
        printf("Running multiplication test %d times...\n", test_num);

        timer_start(&timer);
        for(int i = 0; i < test_num; ++i) {
            MatrixMultiply(&A, &B, &C);
        }

        double elapsed_ms = timer_elapsed_ms(&timer) / test_num;
        printf("Time for multiplication (average for %d times): %.3f ms\n", test_num, elapsed_ms);
        fprintf(resultFile, "%d,%d,%d,%.3f,%d\n", n, m, p, elapsed_ms, test_num);

        MatrixFree(&A); MatrixFree(&B); MatrixFree(&C);
    }
    fclose(resultFile);
}

int main(void)
{
    init_test_scales();

    // test multiply
    test_multiply();

    return 0;
}
