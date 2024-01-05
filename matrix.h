/***************************************************************
 *                  Matrices: Taken from 
 *          https://github.com/simondlevy/TinyEKF
 *                  Creado: 24 sep. 2023
 *                  Autor: Luana Quirolo
 *                  Padron: 102102
****************************************************************/

#ifndef MATRIX_H_
#define MATRIX_H_
#include <math.h>
#include <stdio.h> //TODO: BORRAR
/* Cholesky-decomposition matrix-inversion code, adapted from
   http://jean-pierre.moreau.pagesperso-orange.fr/Cplus/choles_cpp.txt */


int choldc1(int n, double a[n][n], double * p);

int choldcsl(int n, double A[n][n], double a[n][n], double * p);


int cholsl(int n, double A[n][n], double a[n][n], double * p);

// Set all components from matrix A to zero
void mat_zeros(int m, int n, double a[m][n]);

// Set all components from vector A to zero
void vec_zeros(int n, double a[n]);

/* C <- A * B */
void mulmat(int arows, int acols, int bcols, double a[arows][acols], double b[acols][bcols], double c[arows][bcols]);

/* Y <- A * X */
void mulvec(int m, int n, double a[m][n], double * x, double * y);

/* trans(A) <- A */
void transpose(int m, int n, double a[m][n], double at[n][m]);

/* A <- A + B */
void accum(int m, int n, double a[m][n], double b[m][n]);

/* A <- A + B */
void accum_vec(int n, double * a, double * b);

/* C <- A + B */
void add(int m, int n, double a[m][n], double b[m][n], double c[m][n]);

/* C <- A + B */
void add_vec(int n, double * a, double * b, double * c);

/* C <- A - B */
void sub_vec(double * a, double * b, double * c, int n);

/* -A <- A */
void mat_negate(int m, int n, double a[m][n]);

/* A + I <- A */
void mat_addeye(int n, double a[n][n]);

void  mat_copy(int m, int n, double a[m][n], double b[m][n]);

/* bA <- A */
void matmul_scalar(int m, int n, double a[m][n], double b);

/* C <- bA */
void matmul_scalar2(int m, int n, double a[m][n], double c[m][n], double b);

/* bA <- A */
void vecmul_scalar(int n, double *a, double b);

/* C <- bA */
void vecmul_scalar2(int n, double *a, double *c, double b);

/* C <- a * b */ 
// Producto de vectores que produce una matriz
void vec_outer(int n, double *a, double *b, double c[n][n]);

void vec_dot(int n, double *a, double *b, double *c);

void mat_getrow(int m, int n, double a[m][n], int row, double b[n]);

void mat_getcol(int m, int n, double a[m][n], int col, double b[m]);


#endif /* MATRIX_H_ */