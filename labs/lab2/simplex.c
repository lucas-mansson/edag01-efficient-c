#include <stdio.h>
#include <stdlib.h>

typedef struct simplex_t {
        int m;      // Constraints
        int n;      // Decision variables
        int* var;   // 0..n - 1 are nonbasic.
        double** a; // a[m][n+1];
        double* b;  // b[m]
        double* x;  // x[n+1]
        double* c;  // c[n]
        double y;   // y
} simplex_t;

// init - 'constructor' for simplex_t
int init(simplex_t* s, int m, int n, double** a, double* b, double* x,
         double* c, double y, int* var)
{
        int i, k;

        s->m = m;
        s->n = n;
        s->var = var;
        s->a = a;
        s->b = b;
        s->x = x;
        s->c = c;
        s->y = y;

        if (s->var == NULL) {
                s->var = calloc(m + n + 1, sizeof(typeof(m)));
                for (int i = 0; i < m + n; i++) {
                        s->var[i] = i;
                }
        }
        for (int k = 0, i = 1; i < m; i++) {
                if (b[i] < b[k]) {
                        k = i;
                }
        }
        return k;
}
// select_nonbasic
int select_nonbasic(simplex_t* s)
{
        int i;
        for (int i = 0; i < s->n; i++) {
                if (s->c[i] > 0) {
                        return i;
                }
        }
        return -1;
}

// initial and assume bi ≥ 0 so skip the call to prepare and the rest of
int initial(simplex_t* s, int m, int n, double** a, double* b, double* x,
            double* c, double y, int* var)
{
        int i;
        int j;
        int k;

        k = init(s, m, n, a, b, x, c, y, var);
        if (b[k] <= 0) {
                printf("b[k] is not greater than zero: b[%d]; %lf", k, b[k]);
                return 0;
        }

        return 1;
}

// pivot
void pivot(simplex_t* s, int row, int col)
{
        double** a = s->a;
        double* b = s->b;
        double* c = s->c;
        int m = s->m;
        int n = s->n;
        int i;
        int j;
        int t;

        t = s->var[col];
        s->var[col] = s->var[n + row];
        s->var[n + row] = t;
        s->y += (c[col] * b[row]) / (a[row][col]);

        for (i = 0; i < n; i++) {
                if (i != col) {
                        c[i] -= (c[col] * a[row][i]) / a[row][col];
                }
        }

        c[col] = -c[col] / a[row][col];
        for (i = 0; i < m; i++) {
                if (i != row) {
                        b[i] -= (a[i][col] * b[row]) / a[row][col];
                }
        }

        for (i = 0; i < m; i++) {
                if (i == row) {
                        continue;
                }
                for (j = 0; j < m; j++) {
                        if (j == col) {
                                a[i][j] -=
                                    (a[i][col] * a[row][j]) / a[row][col];
                        }
                }
        }

        for (i = 0; i < m; i++) {
                if (i != row) {
                        a[i][col] = -a[i][col] / a[row][col];
                }
        }
        for (i = 0; i < n; i++) {
                if (i != col) {
                        a[row][i] /= a[row][col];
                }
        }
        b[row] /= a[row][col];
        a[row][col] = 1 / a[row][col];
}

// simplex

// xsimplex

// initial
