struct simplex_t {
        int m;      // Constraints
        int n;      // Decision variables
        int* var;   // 0..n - 1 are nonbasic.
        double** a; // a[m][n+1];
        double* b;  // b[m]
        double* x;  // x[n+1]
        double* c;  // c[n]
        double y;   // y
};
