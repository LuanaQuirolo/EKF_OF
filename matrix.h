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

/* Cholesky-decomposition matrix-inversion code, adapted from
   http://jean-pierre.moreau.pagesperso-orange.fr/Cplus/choles_cpp.txt */


int choldc1(double * a, double * p, int n);

int choldcsl(double * A, double * a, double * p, int n);


int cholsl(double * A, double * a, double * p, int n);

// Set all components from matrix A to zero
void mat_zeros(double * a, int m, int n);

/* C <- A * B */
void mulmat(double * a, double * b, double * c, int arows, int acols, int bcols);

/* Y <- A * X */
void mulvec(double * a, double * x, double * y, int m, int n);

/* trans(A) <- A */
void transpose(double * a, double * at, int m, int n);

/* A <- A + B */
void accum(double * a, double * b, int m, int n);

/* C <- A + B */
void add(double * a, double * b, double * c, int m, int n);

/* C <- A - B */
void sub(double * a, double * b, double * c, int n);

/* -A <- A */
void mat_negate(double * a, int m, int n);

/* A + I <- A */
void mat_addeye(double * a, int n);

/* bA <- A */
void matmul_scalar(double * a, int m, int n, double b);

/* C <- bA */
void matmul_scalar2(double * a, double * c, int m, int n, double b);

#endif /* MATRIX_H_ */