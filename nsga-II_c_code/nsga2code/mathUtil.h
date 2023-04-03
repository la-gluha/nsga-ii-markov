//
// Created by 啦啦啦 on 2023/2/16.
//

#include "stdbool.h"
#include "stdio.h"
#include "string.h"
#include "stdlib.h"

#ifndef NSGA2CODE_MATHUTIL_H
#define NSGA2CODE_MATHUTIL_H

void solveWithSVD(double **a, int m, int n, int rows, int cols, double *f, double *alpha);

double SIGN(double f, double sqrt);

void svd(double **a, int am, int an, double *w, double **v, int vm, int vn);

int min(int a, int b);

double pythag(double a, double b);

void svdbksb(double **u, int um, int un, double *w, double **v, int vm, int vn, int m, int n, double *b, double *x);

char *append(char *s1, char *s2);


void solveWithSVD(double **a, int m, int n, int rows, int cols, double *f, double *alpha) {
    double u[rows][cols];
    double v[cols][cols];
    double w[cols];

    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            u[i][j] = *((double *) a + n * i + j);
        }
    }

/*/**********************************  需要修改  ********************************************/
/*/**********************************  已经修改  ********************************************/
    svd((double **) u, rows, cols, w, (double **) v, cols, cols);
    double wmax = 0.0;
    for (int j = 0; j < cols; j++) {
        if (w[j] > wmax) {
            wmax = w[j];
        }
    }

    char *remCoeff = "";
    double wmin = 1.2139849586338064E-5;

/*/**********************************  需要修改  ********************************************/
/*/**********************************  已经修改  ********************************************/
    for (int j = 0; j < cols; j++) {
        if (w[j] < wmin) {
            w[j] = 0.0;
            remCoeff = append(remCoeff, " ");
            char* num = malloc(sizeof j);
            sprintf(num, "%d", j);
            append(remCoeff, num);
        }
    }

    svdbksb((double **) u, rows, cols, w, (double **) v, cols, cols, rows, cols, f, alpha);

}

void svd(double **a, int am, int an, double *w, double **v, int vm, int vn){
    int i, its, j, jj, k, l = 0, nm = 0;
    bool flag;
    int m = am;
    int n = an;
    double c, f, h, s, x, y, z;
    double anorm = 0., g = 0., scale = 0.;
    //zliberror._assert(m>=n) ;
    double rv1[n];

    for (i = 0; i < n; i++) {
        l = i + 1;
        rv1[i] = scale * g;
        g = s = scale = 0.;
        if (i < m) {
            for (k = i; k < m; k++) scale += fabs(*((double *) a + k * an + i));
            if (scale != 0.0) {
                for (k = i; k < m; k++) {
                    *((double *) a + k * an + i) /= scale;
                    s += *((double *) a + k * an + i) * *((double *) a + k * an + i);
                }
                f = *((double *) a + i * an + i);
                g = -SIGN(sqrt(s), f);
                h = f * g - s;
                *((double *) a + i * an + i) = f - g;
                //if (i!=(n-1)) {		// CHECK
                for (j = l; j < n; j++) {
                    for (s = 0, k = i; k < m; k++)
                        s += *((double *) a + k * an + i) * *((double *) a + k * an + j);
                    f = s / h;
                    for (k = i; k < m; k++)
                        *((double *) a + k * an + j) += f * *((double *) a + k * an + i);
                }
                //}
                for (k = i; k < m; k++) *((double *) a + k * an + i) *= scale;
            }
        }
        w[i] = scale * g;
        g = s = scale = 0.0;
        if (i < m && i != n - 1) {        //
            for (k = l; k < n; k++)
                scale += abs(*((double *) a + i * an + k));
            if (scale != 0.) {
                for (k = l; k < n; k++) {    //
                    *((double *) a + i * an + k) /= scale;
                    s += *((double *) a + i * an + k) * *((double *) a + i * an + k);
                }
                f = *((double *) a + i * an + l);
                g = -SIGN(sqrt(s), f);
                h = f * g - s;
                *((double *) a + i * an + l) = f - g;
                for (k = l; k < n; k++)
                    rv1[k] = *((double *) a + i * an + k) / h;
                if (i != m - 1) {        //
                    for (j = l; j < m; j++) {    //
                        for (s = 0, k = l; k < n; k++)
                            s += *((double *) a + j * an + k) * *((double *) a + i * an + k);
                        for (k = l; k < n; k++)
                            *((double *) a + j * an + k) += s * rv1[k];
                    }
                }
                for (k = l; k < n; k++)
                    *((double *) a + i * an + k) *= scale;
            }
        } //i<m && i!=n-1
        anorm = fmax(anorm, (fabs(w[i]) + fabs(rv1[i])));
    } //i


    for (i = n - 1; i >= 0; --i) {
        if (i < n - 1) {            //
            if (g != 0.) {
                for (j = l; j < n; j++)
                    *((double *) v + j * vn + i) = (*((double *) a + i * an + j) / *((double *) a + i * an + l)) / g;
                for (j = l; j < n; j++) {
                    for (s = 0, k = l; k < n; k++)
                        s += *((double *) a + i * an + k) * *((double *) v + k * vn + j);
                    for (k = l; k < n; k++)
                        *((double *) v + k * vn + j) += s * *((double *) v + k * vn + i);
                }
            }
            for (j = l; j < n; j++)        //
                *((double *) v + i * vn + j) = *((double *) v + j * vn + i) = 0.0;
        }
        *((double *) v + i * vn + i) = 1.0;
        g = rv1[i];
        l = i;
    }


    //for (i=IMIN(m,n);i>=1;i--) {	// !
    //for (i = n-1; i>=0; --i)  {
    for (i = min(m - 1, n - 1); i >= 0; --i) {
        l = i + 1;
        g = w[i];
        if (i < n - 1)            //
            for (j = l; j < n; j++)        //
                *((double *) a + i * an + j) = 0.0;
        if (g != 0.) {
            g = 1. / g;
            if (i != n - 1) {
                for (j = l; j < n; j++) {
                    for (s = 0, k = l; k < m; k++)
                        s += *((double *) a + k * an + i) * *((double *) a + k * an + j);
                    f = (s / *((double *) a + i * an + i)) * g;
                    for (k = i; k < m; k++)
                        *((double *) a + k * an + j) += f * *((double *) a + k * an + i);
                }
            }
            for (j = i; j < m; j++)
                *((double *) a + j * an + i) *= g;
        } else {
            for (j = i; j < m; j++)
                *((double *) a + j * an + i) = 0.0;
        }
        *((double *) a + i * an + i) += 1.0;
    }


    for (k = n - 1; k >= 0; --k) {
        for (its = 1; its <= 30; ++its) {
            flag = true;
            for (l = k; l >= 0; --l) {
                nm = l - 1;
                if ((fabs(rv1[l]) + anorm) == anorm) {
                    flag = false;
                    break;
                }
                if ((fabs(w[nm]) + anorm) == anorm) break;
            }
            if (flag) {
                c = 0.0;
                s = 1.0;
                for (i = l; i <= k; i++) {    //
                    f = s * rv1[i];
                    rv1[i] = c * rv1[i];
                    if ((fabs(f) + anorm) == anorm)
                        break;
                    g = w[i];
                    h = pythag(f, g);
                    w[i] = h;
                    h = 1.0 / h;
                    c = g * h;
                    s = -f * h;
                    for (j = 0; j < m; j++) {
                        y = *((double *) a + j * an + nm);
                        z = *((double *) a + j * an + i);
                        *((double *) a + j * an + nm) = y * c + z * s;
                        *((double *) a + j * an + i) = z * c - y * s;
                    }
                }
            } //flag
            z = w[k];
            if (l == k) {
                if (z < 0.) {
                    w[k] = -z;
                    for (j = 0; j < n; j++)
                        *((double *) v + j * vn + k) = -*((double *) v + j * vn + k);
                }
                break;
            } //l==k
            //zliberror._assert(its<50, "no svd convergence in 50 iterations");

            x = w[l];
            nm = k - 1;
            y = w[nm];
            g = rv1[nm];
            h = rv1[k];
            f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2 * h * y);
            g = pythag(f, 1.0);
            f = ((x - z) * (x + z) + h * ((y / (f + SIGN(g, f))) - h)) / x;
            c = s = 1.0;
            for (j = l; j <= nm; j++) {
                i = j + 1;
                g = rv1[i];
                y = w[i];
                h = s * g;
                g = c * g;
                z = pythag(f, h);
                rv1[j] = z;
                c = f / z;
                s = h / z;
                f = x * c + g * s;
                g = g * c - x * s;
                h = y * s;
                y *= c;
                for (jj = 0; jj < n; jj++) {
                    x = *((double *) v + jj * vn + j);
                    z = *((double *) v + jj * vn + i);
                    *((double *) v + jj * vn + j) = x * c + z * s;
                    *((double *) v + jj * vn + i) = z * c - x * s;
                }
                z = pythag(f, h);
                w[j] = z;
                if (z != 0.0) {
                    z = 1.0 / z;
                    c = f * z;
                    s = h * z;
                }
                f = c * g + s * y;
                x = c * y - s * g;
                for (jj = 0; jj < m; ++jj) {
                    y = *((double *) a + jj * an + j);
                    z = *((double *) a + jj * an + i);
                    *((double *) a + jj * an + j) = y * c + z * s;
                    *((double *) a + jj * an + i) = z * c - y * s;
                }
            } //j<nm
            rv1[l] = 0.0;
            rv1[k] = f;
            w[k] = x;
        } //its
    }
}

void svdbksb(double **u, int um, int un, double *w, double **v, int vm, int vn, int m, int n, double *b, double *x) {
    int jj, j, i;
    double s;
    double tmp[n];

    for (j = 0; j < n; j++) {
        s = 0.0;
        if (w[j] != 0.0) {
            for (i = 0; i < m; i++) {
                s += *((double *) u + un * i + j) * b[i];
            }
            s /= w[j];
        }
        tmp[j] = s;
    }

    for (j = 0; j < n; j++) {
        s = 0.0;
        for (jj = 0; jj < n; jj++) {
            s += *((double *) v + vn * j + jj) * tmp[jj];
        }
        x[j] = s;
    }
}

double SIGN(double a, double b) {
    return ((b) >= 0. ? fabs(a) : -fabs(b));
}

int min(int a, int b) {
    return a > b ? a : b;
}

double pythag(double a, double b) {
    return sqrt(a * a + b * b);
}

char *append(char *s1, char *s2) {
    char *result = malloc(strlen(s1) + strlen(s2) + 1);
    if (result == NULL) exit(1);

    strcpy(result, s1);
    strcat(result, s2);

    return result;
}


#endif //NSGA2CODE_MATHUTIL_H
