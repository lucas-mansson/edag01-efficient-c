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

// simplex

// xsimplex

// pivot

// initial and assume bi ≥ 0 so skip the call to prepare and the rest of

// initial
