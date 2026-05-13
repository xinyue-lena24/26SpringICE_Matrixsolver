#include "matrix_core.h"
#include <stdio.h>

int main(void)
{
    Matrix I;
    MatrixInit(&I);

    MatrixError error = MatrixCreate(&I, 3, 3);
    if (error != MATRIX_SUCCESS) {
        printf("Create failed: %s\n", MatrixErrorMessage(error));
        return 1;
    }

    MatrixFillZero(&I);
    for (int i = 0; i < I.row; ++i) {
        MatrixSet(&I, i, i, 1.0);
    }

    MatrixPrint(&I, "I");
    printf("\n");

    REAL value = 0.0;
    error = MatrixGet(&I, 3, 0, &value);
    printf("MatrixGet(&I, 3, 0): %s\n", MatrixErrorMessage(error));

    MatrixFree(&I);
    MatrixFree(&I);
    printf("Finished Lab 2 practice.\n");
    return 0;
}
