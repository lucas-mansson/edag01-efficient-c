#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#define epsilon 1e-6

typedef struct simplex {
        int m;      // Constraints
        int n;      // Decision variables
        int* var;   // 0..n - 1 are nonbasic.
        double** a; // a[m][n+1];
        double* b;  // b[m]
        double* c;  // c[n]
        double* x;  // x[n+1]
        double y;   // y
} simplex_t;

void pivot(simplex_t* s, int row, int col);
double xsimplex(int m, int n, double** a, double* b, double* c, double* x,
                double y, int* var, int h);

int init(simplex_t* s, int m, int n, double** a, double* b, double* c,
         double* x, double y, int* var)
{
        int i;
        int k;

        s->m = m;
        s->n = n;
        s->var = var;
        s->a = a;
        s->b = b;
        s->c = c;
        s->x = x;
        s->y = y;

        if (s->var == NULL) {
                s->var = calloc(m + n + 1, sizeof(int));
                for (int i = 0; i < m + n; i++) {
                        s->var[i] = i;
                }
        }
        for (k = 0, i = 1; i < m; i++) {
                if (s->b[i] < s->b[k]) {
                        k = i;
                }
        }
        return k;
}

int select_nonbasic(simplex_t* s)
{
        for (int i = 0; i < s->n; i++) {
                if (s->c[i] > epsilon) {
                        return i;
                }
        }
        return -1;
}

void prepare(simplex_t* s, int k)
{
        int m = s->m;
        int n = s->n;
        int i;

        // make room for xm+n at s.var[n] by moving s.var[n..n+m-1] one step to
        // the right.
        for (i = m + n; i > n; i--) {
                s->var[i] = s->var[i - 1];
        }
        s->var[n] = m + n;

        // add xm+n to each constraint
        n++;
        for (i = 0; i < m; i++) {
                // s.a[i][n-1] <- -1 ??
                s->a[i][n - 1] = -1; // double check this
        }

        s->x = calloc(m + n, sizeof(double));
        s->c = calloc(m, sizeof(double));
        s->c[n - 1] = -1;
        s->n = n;

        pivot(s, k, n - 1);
}

int initial(simplex_t* s, int m, int n, double** a, double* b, double* c,
            double* x, double y, int* var)
{
        int i;
        int j;
        int k;

        // double w;

        k = init(s, m, n, a, b, c, x, y, var);

        if (s->b[k] >= 0) {
                return 1;
        }

        prepare(s, k);

        n = s->n;
        s->y = xsimplex(m, n, s->a, s->b, s->c, s->x, 0, s->var, 1);

        for (i = 0; i < m + n; i++) {
                if (s->var[i] == m + n - 1) {
                        if (abs((int)s->x[i]) > epsilon) {
                                // delete s.x
                                // delete s.c
                                free(s->x);
                                s->x = NULL;
                                free(s->c);
                                s->c = NULL;
                                return 0;
                        } else {
                                break;
                        }
                }
        }

        if (i >= n) {
                // x_n+m is basic. find good nonbasic.
                int j;
                for (j = k = 0; k < n; k++) {
                        if (fabs(s->a[i - n][k]) > fabs(s->a[i - n][j])) {
                                j = k;
                        }
                }
                pivot(s, i - n, j);
                i = j;
        }
        if (i < n - 1) {
                // xn+m is nonbasic and not last. swap columns i and n-1
                k = s->var[i];
                s->var[i] = s->var[n - 1];
                s->var[n - 1] = k;
                for (k = 0; k < m; k++) {
                        double w = s->a[k][n - 1];

                        s->var[n - 1] = s->a[k][i];
                        s->a[k][i] = w;
                }
        } else {
                // xn+m is nonbasic and last. forget it.
        }
        // delete s.c
        free(s->c);
        s->c = NULL;
        s->c = c;
        s->y = y;

        for (k = n - 1; k < n + m - 1; k++) {
                s->var[k] = s->var[k + 1];
        }

        n = s->n = s->n - 1;
        double* t = calloc(n, sizeof(double));

        int next_k = 0;
        for (k = 0; k < n; k++) {
                for (j = 0; j < n; j++) {
                        if (k == s->var[j]) {
                                // x_k is nonbasic. add ck
                                t[j] = t[j] + s->c[k];
                                next_k = 1; // goto next_k?
                                break;
                        }
                }
                if (next_k == 1) {
                        next_k = 0;
                        continue; // skip rest of this k-loop
                }
                // xk is basic.
                for (j = 0; j < m; j++) {
                        if (s->var[n + j] == k) {
                                // x_k is at row j
                                break;
                        }
                }
                s->y = s->y + (s->c[k] * s->b[j]);
                for (i = 0; i < n; i++) {
                        t[i] = t[i] - (s->c[k] * s->a[j][i]);
                }
        }
        for (i = 0; i < n; i++) {
                s->c[i] = t[i];
        }
        // delete t and s.x ????
        free(t);
        t = NULL;
        free(s->x);
        s->x = NULL;

        return 1;
}

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
                        b[i] = b[i] - (a[i][col] * b[row]) / a[row][col];
                }
        }

        for (i = 0; i < m; i++) {
                if (i == row) {
                        continue;
                }
                for (j = 0; j < n; j++) {
                        if (j != col) {
                                a[i][j] = a[i][j] - ((a[i][col] * a[row][j]) /
                                                     a[row][col]);
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
        b[row] = b[row] / a[row][col];
        a[row][col] = 1 / a[row][col];
}

double xsimplex(int m, int n, double** a, double* b, double* c, double* x,
                double y, int* var, int h)
{
        simplex_t* s = calloc(1, sizeof(simplex_t));
        int i;
        int row;
        int col;

        if (!initial(s, m, n, a, b, c, x, y, var)) {
                free(s->var);
                s->var = NULL;
                free(s);
                return -1;
        }

        while ((col = select_nonbasic(s)) >= 0) {
                row = -1;
                for (i = 0; i < m; i++) {
                        if (a[i][col] > epsilon &&
                            (row < 0 ||
                             ((b[i] / a[i][col]) < (b[row] / a[row][col])))) {
                                row = i;
                        }
                }
                if (row < 0) {
                        free(s->var);
                        s->var = NULL;
                        free(s);
                        return INFINITY;
                }
                pivot(s, row, col);
        }

        if (h == 0) {
                for (i = 0; i < n; i++) {
                        if (s->var[i] < n) {
                                x[s->var[i]] = 0;
                        }
                }
                for (i = 0; i < n; i++) {
                        if (s->var[n + i] < n) {
                                x[s->var[n + i]] = s->b[i];
                        }
                }
                free(s->var);
                s->var = NULL;
        } else {
                for (i = 0; i < n; i++) {
                        x[i] = 0;
                }
                for (i = n; i < n + m; i++) {
                        x[i] = s->b[i - n];
                }
        }
        double res = s->y;
        free(s);
        s = NULL;
        return res;
}

double simplex(int m, int n, double** a, double* b, double* c, double* x,
               double y)
{
        return xsimplex(m, n, a, b, c, x, y, NULL, 0);
}

int main()
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
        double** a = calloc(m + n, sizeof(double*));
        for (int i = 0; i < m; i++) {
                a[i] = calloc(n + 1, sizeof(double));
                for (int j = 0; j < n; j++) {
                        scanf("%lf", &a[i][j]);
                }
        }

        // last is a vector with m number of b-values
        double* b = calloc(m, sizeof(double));
        for (int i = 0; i < m; i++) {
                b[i] = 0;
                scanf("%lf", &b[i]);
        }

        for (int i = 0; i < m; i++) {
                for (int j = 0; j < n; j++) {
                        if (j > 0) {
                                printf(" + ");
                        }
                        printf("%2.1lfx_%d", a[i][j], j);
                }
                printf(" <= %2.1lf", b[i]);
                printf("\n");
        }

        double* x = calloc(m + n, sizeof(double));
        double y = 0.0;

        double res = simplex(m, n, a, b, c, x, y);

        printf("Result: %lf \n", res);

        free(c);
        c = NULL;
        for (int i = 0; i < m; i++) {
                free(a[i]);
                a[i] = NULL;
        }
        free(a);
        a = NULL;
        free(b);
        b = NULL;
        free(x);
        x = NULL;

        return 0;
}
