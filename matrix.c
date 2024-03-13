/***************************************************************
 *                  Matrices: Taken from 
 *          https://github.com/simondlevy/TinyEKF
 *                  Creado: 24 sep. 2023
 *                  Autor: Luana Quirolo
 *                  Padron: 102102
****************************************************************/

#include "matrix.h"

/* Cholesky-decomposition matrix-inversion code, adapted from
   http://jean-pierre.moreau.pagesperso-orange.fr/Cplus/choles_cpp.txt */


int choldc1(int n, float a[n][n], float * p){
    int i,j,k;
    float sum;

    for (i = 0; i < n; i++) {
        for (j = i; j < n; j++) {
            sum = a[i][j];
            for (k = i - 1; k >= 0; k--) {
                sum -= (float)a[i][k] * a[j][k];
            }
            if (i == j) {
                if (sum <= 0) {
                    return 1; /* error */
                }
                p[i] = (float)sqrtf((float)sum);
            }
            else {
                a[j][i] = (float)sum / p[i];
            }
        }
    }

    return 0; /* success */
}

int choldcsl(int n, float A[n][n], float a[n][n], float * p){
    int i,j,k; float sum;
    for (i = 0; i < n; i++) 
        for (j = 0; j < n; j++) 
            a[i][j] = (float)A[i][j];
    if (choldc1(n, a, p)) return 1;
    for (i = 0; i < n; i++) {
        a[i][i] = (float)1 / p[i];
        for (j = i + 1; j < n; j++) {
            sum = 0;
            for (k = i; k < j; k++) {
                sum -= a[j][k] * a[k][i];
            }
            a[j][i] = (float)sum / p[j];
        }
    }

    return 0; /* success */
}


int cholsl(int n, float A[n][n], float a[n][n], float * p){
    int i,j,k;
    if (choldcsl(n,A,a,p)) return 1;
    for (i = 0; i < n; i++) {
        for (j = i + 1; j < n; j++) {
            a[i][j] = 0.0;
        }
    }
    for (i = 0; i < n; i++) {
        a[i][i] *= a[i][i];
        for (k = i + 1; k < n; k++) {
            a[i][i] += a[k][i] * a[k][i];
        }
        for (j = i + 1; j < n; j++) {
            for (k = j; k < n; k++) {
                a[i][j] += a[k][i] * a[k][j];
            }
        }
    }
    for (i = 0; i < n; i++) {
        for (j = 0; j < i; j++) {
            a[i][j] = a[j][i];
        }
    }

    return 0; /* success */
}

// Set all components from matrix A to zero
void mat_zeros(int m, int n, float a[m][n]){
    for(int i = 0; i < m; i++){
        for(int j = 0; j < n; j++) {
            a[i][j] = 0.0f;
        }
    }
}

// Set all components from vector A to zero
void vec_zeros(int n, float a[n]){
    for(int j = 0; j < n; j++){
        a[j] = 0.0f;
    }
}

/* C <- A * B */
void mulmat(int arows, int acols, int bcols, float a[arows][acols], float b[acols][bcols], float c[arows][bcols]){
    for(int i = 0; i < arows; i++)
        for(int j = 0; j < bcols; j++) {
            c[i][j] = 0;
            for(int l = 0; l < acols; l++){
                c[i][j] +=  (float)a[i][l] * b[l][j];
            }
        }
}

/* Y <- A * X */
void mulvec(int m, int n, float a[m][n], float * x, float * y){
    for(int i = 0;  i <m; i++) {
        y[i] = 0.0f;
        for(int j = 0; j < n; j++){
            y[i] += (float) x[j] * a[i][j];
        }
    }
}

/* trans(A) <- A */
void transpose(int m, int n, float a[m][n], float at[n][m]){
    for(int i = 0; i < m; i++){
        for(int j = 0; j < n; j++) {
            at[j][i] =  (float)a[i][j];
        }
    }
}

/* A <- A + B */
void accum(int m, int n, float a[m][n], float b[m][n]){        
    for(int i = 0; i < m; i++){
        for(int j = 0; j < n; j++){
            a[i][j] += (float)b[i][j];
        }
    }
}

/* A <- A + B */
void accum_vec(int n, float * a, float * b){
    for(int j = 0; j < n; j++){
        a[j] +=(float) b[j];
    }
}

/* C <- A + B */
void add(int m, int n, float a[m][n], float b[m][n], float c[m][n]){
    for(int i = 0; i < m; i++){
        for(int j = 0; j < n; j++){
            c[i][j] = (float)b[i][j] + a[i][j];
        }
    }
}

/* C <- A + B */
void add_vec(int n, float * a, float * b, float * c){
    for(int j = 0; j < n; j++){
        c[j] =(float) a[j] +  (float)b[j];
    }
}

/* C <- A - B */
void sub_vec(float * a, float * b, float * c, int n){
    for(int j = 0; j < n; j++){
        c[j] = (float) a[j] -  (float)b[j];
    }    
}

/* -A <- A */
void mat_negate(int m, int n, float a[m][n]){        
    for(int i = 0; i < m; i++){
        for(int j = 0; j < n; j++){
            a[i][j] = (float)-a[i][j];
        }
    }
}

/* A + I <- A */
void  mat_addeye(int n, float a[n][n]){
    for (int i = 0; i < n; i++){
        a[i][i] += (float)1.0f;
    }
}

/* B <- A */
void  mat_copy(int m, int n, float a[m][n], float b[m][n]){
    for(int i = 0; i < m; i++){
        for(int j = 0; j < n; j++){
            b[i][j] =(float) a[i][j];
        }
    }
}

/* bA <- A */
void matmul_scalar(int m, int n, float a[m][n], float b){        
    for(int i = 0; i < m; i++){
        for(int j = 0; j < n; j++){
            a[i][j] =(float) b * a[i][j];
        }
    }
}

/* C <- bA */
void matmul_scalar2(int m, int n, float a[m][n], float c[m][n], float b){        
    for(int i = 0; i < m; i++){
        for(int j = 0; j < n; j++){
            c[i][j] =(float) b * a[i][j];
        }
    }
}

/* bA <- A */
void vecmul_scalar(int n, float *a, float b){        
    for(int i = 0; i < n; i++){
        a[i] = (float)b * a[i];
    }
}

/* C <- bA */
void vecmul_scalar2(int n, float *a, float *c, float b){        
    for(int i = 0; i < n; i++){
        c[i] = (float)b * a[i];
    }
}

/* C <- a * b */ 
// Producto de vectores que produce una matriz
void vec_outer(int n, float *a, float *b, float c[n][n]){
    for(int i = 0; i < n; i++){
        for(int j = 0; j < n; j++){
            c[i][j] = (float)a[i] * b[j];
        }
    }
}

/* C <- a * b */ 
// Producto de vectores que produce un escalar
void vec_dot(int n, float *a, float *b, float *c){
    *c = 0;
    for(int i = 0; i < n; i++){
        *c += (float)a[i] * b[i];
    }
}

void mat_getrow(int m, int n, float a[m][n], int row, float b[n]){
    for (int i = 0; i < n; i++){
        b[i] = (float)a[row][i];
    }
}

void mat_getcol(int m, int n, float a[m][n], int col, float b[m]){
    for (int i = 0; i < m; i++){
        b[i] = (float)a[i][col];
    }
}