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


int choldc1(int n, float a[n][n], float * p);

int choldcsl(int n, float A[n][n], float a[n][n], float * p);


int cholsl(int n, float A[n][n], float a[n][n], float * p);

// Set all components from matrix A to zero
void mat_zeros(int m, int n, float a[m][n]);

// Set all components from vector A to zero
void vec_zeros(int n, float a[n]);

/* C <- A * B */
void mulmat(int arows, int acols, int bcols, float a[arows][acols], float b[acols][bcols], float c[arows][bcols]);

/* Y <- A * X */
void mulvec(int m, int n, float a[m][n], float * x, float * y);

/* trans(A) <- A */
void transpose(int m, int n, float a[m][n], float at[n][m]);

/* A <- A + B */
void accum(int m, int n, float a[m][n], float b[m][n]);

/* A <- A + B */
void accum_vec(int n, float * a, float * b);

/* C <- A + B */
void add(int m, int n, float a[m][n], float b[m][n], float c[m][n]);

/* C <- A + B */
void add_vec(int n, float * a, float * b, float * c);

/* C <- A - B */
void sub_vec(float * a, float * b, float * c, int n);

/* -A <- A */
void mat_negate(int m, int n, float a[m][n]);

/* A + I <- A */
void mat_addeye(int n, float a[n][n]);

void  mat_copy(int m, int n, float a[m][n], float b[m][n]);

/* bA <- A */
void matmul_scalar(int m, int n, float a[m][n], float b);

/* C <- bA */
void matmul_scalar2(int m, int n, float a[m][n], float c[m][n], float b);

/* bA <- A */
void vecmul_scalar(int n, float *a, float b);

/* C <- bA */
void vecmul_scalar2(int n, float *a, float *c, float b);

/* C <- a * b */ 
// Producto de vectores que produce una matriz
void vec_outer(int n, float *a, float *b, float c[n][n]);

void vec_dot(int n, float *a, float *b, float *c);

void mat_getrow(int m, int n, float a[m][n], int row, float b[n]);

void mat_getcol(int m, int n, float a[m][n], int col, float b[m]);


#endif /* MATRIX_H_ */