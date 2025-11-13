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
        printf("Rows: %d Cols: %d \n", m, n);

        printf("max z = ");

        // then we have a vector with n number of c-coefficients
        double* c = calloc(n, sizeof(double));
        for (int i = 0; i < n; i++) {
                if (i > 0) {
                        printf(" + ");
                }
                scanf("%lf", &c[i]);
                printf("%2.1lfx_%d", c[i], i);
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
        double* b = calloc(n, sizeof(double));
        for (int i = 0; i < n; i++) {
                scanf("%lf", &b[i]);
        }

        for (int i = 0; i < m; i++) {
                for (int j = 0; j < n; j++) {
                        if (j > 0) {
                                printf(" + ");
                        }
                        printf("%2.1lfx_%d", matrix[i][j], j);
                }
                printf(" <= %2.1lf", b[i]);
                printf("\n");
        }

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
