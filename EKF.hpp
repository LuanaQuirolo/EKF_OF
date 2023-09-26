/***************************************************************
 *            Algoritmo de fusion de sensores
 *                  Creado: 11 sep. 2023
 *                  Autor: Luana Quirolo
 *                  Padron: 102102
****************************************************************/

#ifndef EKF_H_
#define EKF_H_

#include <math.h>
#include <stdio.h>
#include "matrix.hpp"

#define N_STATES 19
#define N_OBS_00 6
#define N_OBS_01 8
#define N_OBS_10 7
#define N_OBS_11 9

typedef struct mediciones{
  double dt;
  double ax;
  double ay;
  double az;
  double wx;
  double wy;
  double wz;
  double ofx;
  double ofy;
  double range;
  } mediciones_t;

/***************************************************************
 * 				   DEFINICION ESTRUCTURA OFS_EKF
 ***************************************************************/
typedef struct ofs_ekf {
    uint8_t N; // Cantidad de estados
    bool beta; // Indica si hay una lectura nueva del OFS
    bool gamma; // Indica si hay una lectura nueva del sensor de distancia
    double g[3]; 
    uint8_t M00; // Cantidad de observaciones con beta=gamma=0
    uint8_t M01; // Cantidad de observaciones con beta=0, gamma=1
    uint8_t M10; // Cantidad de observaciones con beta=1, gamma=0
    uint8_t M11; // Cantidad de observaciones con beta=1, gamma=1
    double states[N_STATES]; //p, v, q, w, ba, bw
    double cov[N_STATES][N_STATES]; // Matriz de covarianza de estados
    uint8_t Npix; // Cantidad de píxeles
    float FOV_OF; // FOV del sensor de OF
    float f;  // Factor de conversión
} ofs_ekf_t; 

  void ofs_ekf_init(ofs_ekf_t* filtro);

/* Realiza el paso de prediccion. Es decir, propaga la accion de control.
* Pre: El filtro debe estar inicializado.
* Pos: El estado se encuentra actualizado segun el control aplicado. */
  void prediction_step(ofs_ekf_t* filtro, mediciones_t u, double dt);

/* Realiza el paso de correccion luego del paso de prediccion. Es decir,
* una vez propagada la accion de control, se corrige el estado en base a
* las mediciones disponibles.
* Pre: El filtro debe estar inicializado y haber llamado anteriormente a 'prediction_step'.
* Pos: El estado se encuentra corregido en base a las mediciones. */
  void correction_step(ofs_ekf_t* filtro, mediciones_t z, double dt);

#endif /* EKF_H_ */