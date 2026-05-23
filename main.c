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

int test_scale_num;
testScale test_scales[31*31*31];

void init_test_scales() {
    srand(0); // Reset random seed for reproducibility

    // test scales required
    test_scales[0] = (testScale){100, 100, 100};
    test_scales[1] = (testScale){1000, 1000, 1000};
    test_scales[2] = (testScale){2000, 2000, 2000};
    test_scale_num = 3;
    // test_scale_num = 0;

    // Generate test scales with n, m, p
    // for (int n = 16; n <= 1024; n *= 4) {
    //     for (int m = 16; m <= 1024; m *= 4) {
    //         for (int p = 16; p <= 1024; p *= 4) {
    //             test_scales[test_scale_num].n = n;
    //             test_scales[test_scale_num].m = m;
    //             test_scales[test_scale_num].p = p;
    //             test_scale_num++;
    //         }
    //     }
    // }

    // m small
    // for (int n = 32; n <= 2048; n *= 4) {
    //     for (int m = 1; m <= 512; m *= 2) {
    //         for (int p = 32; p <= 2048; p *= 4) {
    //             test_scales[test_scale_num].n = n;
    //             test_scales[test_scale_num].m = m;
    //             test_scales[test_scale_num].p = p;
    //             test_scale_num++;
    //         }
    //     }
    // }

    // random scales
    // for (int i = 0; i < test_scale_num; ++i) {
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
    
    for (int i = 0; i < test_scale_num; ++i) {
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

        long long scale = (long long)n * m * p;
        printf("Scale (n*m*p): %lld\n", scale);
        int test_num = (1e8 / scale) + 5; // Adjust test_num based on the scale
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

void test_all_four_operations()
{
    struct Timer timer;
    Matrix A, B, C;

    srand(0); // Reset random seed for reproducibility
    
    // addition test
    FILE *resultFile = fopen("addition_results.csv", "w+");
    fprintf(resultFile, "n,m,elapsed_ms,test_num\n");
    for (int i = 0; i < test_scale_num; ++i) {
        int n = test_scales[i].n;
        int m = test_scales[i].m;
        MatrixInit(&A); MatrixInit(&B); MatrixInit(&C);
        MatrixGenerateRandom(&A, n, m, 0.0, 1.0);
        MatrixGenerateRandom(&B, n, m, 0.0, 1.0);
        MatrixCreate(&C, n, m);

        int test_num = 2000 * 2000 / (n * m) * 100; 
        
        printf("Test %d: A(%d x %d) + B(%d x %d) = C(%d x %d)\n", i + 1, n, m, n, m, n, m);
        timer_start(&timer);
        for(int i = 0; i < test_num; ++i) {
            MatrixAdd(&A, &B, &C);
        }
        double elapsed_ms = timer_elapsed_ms(&timer) / test_num;
        printf("Time for addition (average for %d times): %.3f ms\n", test_num, elapsed_ms);
        MatrixFree(&A); MatrixFree(&B); MatrixFree(&C);
        fprintf(resultFile, "%d,%d,%.3f,%d\n", n, m, elapsed_ms, test_num);
    }
    fclose(resultFile);

    // scaling test
    resultFile = fopen("scaling_results.csv", "w+");
    fprintf(resultFile, "n,m,elapsed_ms,test_num\n");
    for (int i = 0; i < test_scale_num; ++i) {
        int n = test_scales[i].n;
        int m = test_scales[i].m;
        MatrixInit(&A); MatrixInit(&C);
        MatrixGenerateRandom(&A, n, m, 0.0, 1.0);
        MatrixCreate(&C, n, m);
        
        int test_num = 2000 * 2000 / (n * m) * 100; 
        
        printf("Test %d: B = alpha * A, A(%d x %d) -> B(%d x %d)\n", i + 1, n, m, n, m);
        REAL alpha = rand() / RAND_MAX * 2.0 - 1.0; // Random alpha between -1 and 1
        timer_start(&timer);
        for(int i = 0; i < test_num; ++i) {
            MatrixScale(alpha, &A, &C);
        }
        double elapsed_ms = timer_elapsed_ms(&timer) / test_num;
        printf("Time for scaling (average for %d times): %.3f ms\n", test_num, elapsed_ms);
        fprintf(resultFile, "%d,%d,%.3f,%d\n", n, m, elapsed_ms, test_num);
        MatrixFree(&A); MatrixFree(&C);
    }
    fclose(resultFile);

    // transpose test
    resultFile = fopen("transpose_results.csv", "w+");
    fprintf(resultFile, "n,m,elapsed_ms,test_num\n");
    for (int i = 0; i < test_scale_num; ++i) {
        int n = test_scales[i].n;
        int m = test_scales[i].m;
        MatrixInit(&A); MatrixInit(&C);
        MatrixGenerateRandom(&A, n, m, 0.0, 1.0);
        MatrixCreate(&C, m, n);
        
        int test_num = 2000 * 2000 / (n * m) * 100; 
        
        printf("Test %d: AT = A^T, A(%d x %d) -> AT(%d x %d)\n", i + 1, n, m, m, n);
        timer_start(&timer);
        for(int i = 0; i < test_num; ++i) {
            MatrixTranspose(&A, &C);
        }
        double elapsed_ms = timer_elapsed_ms(&timer) / test_num;
        printf("Time for transpose (average for %d times): %.3f ms\n", test_num, elapsed_ms);
        fprintf(resultFile, "%d,%d,%.3f,%d\n", n, m, elapsed_ms, test_num);
        MatrixFree(&A); MatrixFree(&C);
    }
    fclose(resultFile);

    // multiplication test
    test_multiply();
}

int main(void)
{
    init_test_scales();

    // test all four functions
    test_all_four_operations();

    // test multiply
    // test_multiply();

    return 0;
}
