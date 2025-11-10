#include <stdio.h>
#include <stdlib.h>

double** make_matrix(int rows, int cols)
{
        double** a;
        int i;

        a = calloc(rows, sizeof(double*));
        for (i = 0; i < rows; i++) {
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

        // then we have a vector with n number of c-coefficients

        // then we have a matrix with m rows and n cols
        double** matrix = make_matrix(m, n);
        for (int i = 0; i < m; i++) {
                for (int j = 0; i < n; i++) {
                        double a;
                        scanf("%lf", &a);
                        printf("%lf ", a);
                }
        }

        // last is a vector with m number of b-values

        return 0;
}
