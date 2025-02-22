#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// Function declarations
double **Mateigval(double **a);
double **Matquad(double a, double b, double c);
double **Mateye(int m);
double **createMat(int m, int n);
void printMat(double **p, int m, int n);
double **loadMat(char *str, int m, int n);
double Matnorm(double **a, int m);
double Matdot(double **a, double **b, int m);
double **Matsub(double **a, double **b, int m, int n);
double **Matadd(double **a, double **b, int m, int n);
double **Matscale(double **a, int m, int n, double k);
double **Matinv(double **mat, int m);
double **Matmul(double **a, double **b, int m, int n, int p);
double **Mathstack(double **a, double **b, int m, int n, int p);
double **transposeMat(double **a, int m, int n);
double **rotMat(double theta);
double **normVec(double **a);
double **circulantMat(double **a, int m);
double **Matsec(double **a, double **b, int m, double k);
double Mattrace(double **a, int m);
double Matdet(double **a);
double **Matcol(double **a, int m, int n);
double **Matrow(double **a, int m, int n);
double **Matunit(double **a, int m);
double **Mateigvec(double **a, double **lam);
double **Matquad(double a, double b, double c);


double **Matsec(double **a, double **b, int m, double k) {
    double **temp = Matscale(Matadd(a, Matscale(b, m, 1, k), m, 1), m, 1, 1.0 / (k + 1));
    return temp;
}

double **Matadd(double **a, double **b, int m, int n) {
    int i, j;
    double **c;
    c = createMat(m, n);

    for (i = 0; i < m; i++) {
        for (j = 0; j < n; j++) {
            c[i][j] = a[i][j] + b[i][j];
        }
    }
    return c;
}

double **Matscale(double **a, int m, int n, double k) {
    int i, j;
    double **c;
    c = createMat(m, n);

    for (i = 0; i < m; i++) {
        for (j = 0; j < n; j++) {
            c[i][j] = k * a[i][j];
        }
    }
    return c;
}

double **normVec(double **m) {
    double **temp = Matmul(rotMat(M_PI / 2), m, 2, 2, 1);
    return temp;
}

double **createMat(int m, int n) {
    int i;
    double **a;

    // Allocate memory to the pointer
    a = (double **)malloc(m * sizeof(*a));
    for (i = 0; i < m; i++)
        a[i] = (double *)malloc(n * sizeof(*a[i]));

    return a;
}

double **Matcol(double **a, int m, int n) {
    int i;
    double **b = createMat(m, 1);

    // Extract column vector
    for (i = 0; i < m; i++) {
        b[i][0] = a[i][n];
    }
    return b;
}

double **Matrow(double **a, int m, int n) {
    int i;
    double **b = createMat(1, n);

    // Extract row vector
    for (i = 0; i < n; i++) {
        b[0][i] = a[m][i];
    }
    return b;
}

double **loadMat(char *str, int m, int n) {
    FILE *fp;
    double **a;
    int i, j;

    a = createMat(m, n);
    fp = fopen(str, "r");

    for (i = 0; i < m; i++) {
        for (j = 0; j < n; j++) {
            fscanf(fp, "%lf", &a[i][j]);
        }
    }

    fclose(fp);
    return a;
}

void printMat(double **p, int m, int n) {
    int i, j;

    for (i = 0; i < m; i++) {
        for (j = 0; j < n; j++)
            printf("%lf ", p[i][j]);
        printf("\n");
    }
}

double **rotMat(double theta) {
    double **temp = createMat(2, 2);
    double c = cos(theta), s = sin(theta);
    temp[0][0] = c;
    temp[0][1] = -s;
    temp[1][0] = s;
    temp[1][1] = c;

    return temp;
}

double Matdot(double **a, double **b, int m) {
    double **temp = Matmul(transposeMat(a, m, 1), b, 1, m, 1);
    return temp[0][0];
}

double Matnorm(double **a, int m) {
    return sqrt(Matdot(a, a, m));
}

double **Matsub(double **a, double **b, int m, int n) {
    int i, j;
    double **c;
    c = createMat(m, n);

    for (i = 0; i < m; i++) {
        for (j = 0; j < n; j++) {
            c[i][j] = a[i][j] - b[i][j];
        }
    }
    return c;
}

double **Matinv(double **a, int m) {
    double **c, det = 0;
    int i, j;
    c = createMat(m, m);

    if (m == 2) {
        det = a[0][0] * a[1][1] - a[0][1] * a[1][0];

        c[0][0] = a[1][1] / det;
        c[0][1] = -a[1][0] / det;
        c[1][0] = -a[0][1] / det;
        c[1][1] = a[0][0] / det;
    } else if (m == 3) {
        for (i = 0; i < m; i++)
            det += a[0][i] * (a[1][(i + 1) % 3] * a[2][(i + 2) % 3] - a[1][(i + 2) % 3] * a[2][(i + 1) % 3]);

        for (i = 0; i < m; i++) {
            for (j = 0; j < m; j++)
                c[i][j] = ((a[(i + 1) % 3][(j + 1) % 3] * a[(i + 2) % 3][(j + 2) % 3]) - (a[(i + 1) % 3][(j + 2) % 3] * a[(i + 2) % 3][(j + 1) % 3])) / det;
        }
    } else {
        printf("Invalid input \n");
        exit(0);
    }

    return c;
}

double **Matmul(double **a, double **b, int m, int n, int p) {
    int i, j, k;
    double **c, temp = 0;
    c = createMat(m, p);

    for (i = 0; i < m; i++) {
        for (k = 0; k < p; k++) {
            for (j = 0; j < n; j++) {
                temp = temp + a[i][j] * b[j][k];
            }
            c[i][k] = temp;
            temp = 0;
        }
    }
    return c;
}

double **transposeMat(double **a, int m, int n) {
    int i, j;
    double **c;
    c = createMat(n, m);

    for (i = 0; i < n; i++) {
        for (j = 0; j < m; j++) {
            c[i][j] = a[j][i];
        }
    }
    return c;
}

double Mattrace(double **a, int m) {
    double c = 0;
    for (int i = 0; i < m; i++) {
        c += a[i][i];
    }
    return c;
}

double Matdet(double **a) {
    return a[0][0] * a[1][1] - a[0][1] * a[1][0];
}

double **Mateigval(double **a) {
    double b = -Mattrace(a, 2);
    double c = Matdet(a);
    return Matquad(1, b, c);
}

double **Mateye(int m) {
    int i, j; //dummy integers
    double **I = createMat(m, m);

    for (i = 0; i < m; i++) {
        for (j = 0; j < m; j++) {
            if (i == j)
                I[i][j] = 1;
            else
                I[i][j] = 0;
        }
    }
    return I;
}

double **Matunit(double **a, int m) {
    double **c = createMat(2, 1);
    double temp = Matnorm(a, m);
    for (int i = 0; i < m; i++) {
        c[i][0] = a[i][0] / temp;
    }
    return c;
}

double **Matquad(double a, double b, double c) {
    double **lam = createMat(2, 1);
    double D = sqrt(pow(b, 2.0) - 4 * a * c);
    double den = 2.0 * a;
    lam[0][0] = (-b + D) / den;
    lam[1][0] = (-b - D) / den;
    return lam;
}

double **Mateigvec(double **a, double **lam) {
    double **b1, **b2;
    double **p1, **p2;
    double **temp1, **temp2;
    double **omat = rotMat(M_PI / 2);
    // A-lambda I
    b1 = Matadd(a, Matscale(Mateye(2), 2, 2, -lam[0][0]), 2, 2);
    b2 = Matadd(a, Matscale(Mateye(2), 2, 2, -lam[1][0]), 2, 2);
    // Extract 1st row
    temp1 = Matrow(b1, 0, 2);
    temp2 = Matrow(b2, 0, 2);
    // free the matrices
    free(b1);
    free(b2);
    // Generate unit vector
    b1 = Matunit(temp1, 2);
    b2 = Matunit(temp2, 2);
    // free temp vectors
    free(temp1);
    free(temp2);
    // Find eigen vector
    p1 = Matmul(omat, b1, 2, 2, 1);
    p2 = Matmul(omat, b2, 2, 2, 1);
    // free vectors
    free(b1);
    free(b2);
    return Mathstack(p1, p2, 2, 1, 1);
}

double **Mathstack(double **a, double **b, int m, int n, int p) {
    double **c = createMat(m, n + p);
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            c[i][j] = a[i][j];
        }
        for (int k = 0; k < p; k++) {
            c[i][n + k] = b[i][k];
        }
    }
    return c;
}

double **circulantMat(double **a, int m) {
    int i, j, k;
    printMat(a, 4, 1);
    double **c = createMat(m, m);
    for (i = 0; i < m; i++) {
        for (j = 0; j < m; j++) {
            if (i >= j)
                c[i][j] = a[i - j][0];
            else
                c[i][j] = a[m - j + i][0];
        }
    }
    return c;
}

