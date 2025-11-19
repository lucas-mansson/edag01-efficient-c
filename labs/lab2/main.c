#include "simplex.h"
#include <stdio.h>
#include <stdlib.h>

double** alloc_matrix(int rows, int cols)
{
        double** a;

        a = calloc(rows, sizeof(double*));
        for (int i = 0; i < rows; i++) {
                a[i] = calloc(cols, sizeof(double));
        }
        return a;
}

int main(int argc, char** argv)
{
        // First two numbers are m and n
        int m, n;

        scanf("%d %d", &m, &n);

        // then we have a vector with n number of c-coefficients
        double* c = calloc(n, sizeof(double));
        for (int i = 0; i < n; i++) {
                scanf("%lf", &c[i]);
        }
        printf("\n");

        // then we have a matrix with m rows and n cols
        double** matrix = alloc_matrix(m, n);
        for (int i = 0; i < m; i++) {
                for (int j = 0; j < n; j++) {
                        scanf("%lf", &matrix[i][j]);
                }
        }

        // last is a vector with m number of b-values
        double* b = calloc(m, sizeof(double));
        for (int i = 0; i < m; i++) {
                scanf("%lf", &b[i]);
        }

        double* x = calloc(n, sizeof(double));
        double y = 0;

        double res = simplex(m, n, matrix, b, c, x, y);

        printf("Result: %lf \n", res);

        free(c);
        c = NULL; // good practice to set pointers to null after free
        for (int i = 0; i < m; i++) {
                free(matrix[i]);
                matrix[i] = NULL;
        }
        free(matrix);
        c = NULL;

        free(b);
        b = NULL;

        return 0;
}
