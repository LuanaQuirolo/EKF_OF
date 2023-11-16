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
#include "matrix.h"
#include "quaternions.h"
#include <stdint.h>

#define N_OBS_00 6
#define N_OBS_01 8
#define N_OBS_10 7
#define N_OBS_11 9
#define N_P 3
#define N_V 3
#define N_Q 4
#define N_STATES 10 //N_P+N_V+N_Q -> La cantidad de estados total es la suma de sus componentes
#define N_PROC_NOISE 7 // a, w, z
#define N_CORR_NOISE 9 // a, w, flow_x, flow_y, z
#define N_IMU 6
#define N_OFS 2
#define N_TOFS 1
#define g 9.81

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
    int beta; // Indica si hay una lectura nueva del OFS
    int gamma; // Indica si hay una lectura nueva del sensor de distancia
    quaternion_t qg; 
    uint8_t M00; // Cantidad de observaciones con beta=gamma=0
    uint8_t M01; // Cantidad de observaciones con beta=0, gamma=1
    uint8_t M10; // Cantidad de observaciones con beta=1, gamma=0
    uint8_t M11; // Cantidad de observaciones con beta=1, gamma=1
    double states[N_STATES]; //p, v, q
    double cov[N_STATES][N_STATES]; // Matriz de covarianza de estados
    double F[N_STATES][N_STATES]; // Derivada de vector de estados respecto de si mismo
    double W[N_STATES][N_PROC_NOISE]; // Derivada de vector de estados respecto de ruidos
    double H[N_CORR_NOISE][N_STATES]; // Derivada de modelo de medicion respecto a los estados
    double Q[N_PROC_NOISE][N_PROC_NOISE]; // Matriz ruidos
    uint8_t Npix; // Cantidad de píxeles
    float FOV_OF; // FOV del sensor de OF
    float f;  // Factor de conversión
    int puntero;
} ofs_ekf_t; 

  void ofs_ekf_init(ofs_ekf_t* filtro);

/* Realiza el paso de prediccion. Es decir, propaga la accion de control.
* Pre: El filtro debe estar inicializado.
* Pos: El estado se encuentra actualizado segun el control aplicado. */
  void prediction_step(ofs_ekf_t* filtro, mediciones_t u);

/* Realiza el paso de correccion luego del paso de prediccion. Es decir,
* una vez propagada la accion de control, se corrige el estado en base a
* las mediciones disponibles.
* Pre: El filtro debe estar inicializado y haber llamado anteriormente a 'prediction_step'.
* Pos: El estado se encuentra corregido en base a las mediciones. */
  void correction_step(ofs_ekf_t* filtro, mediciones_t z, double dt);

/* Calcula el jacobiano del modelo de la IMU respecto al estado
* Pre: El filtro debe estar inicializado y haber llamado anteriormente a 'prediction_step'.
* Pos: El jacobiano tiene calculada la parte de la IMU posicionada en las filas segun el offset. */
  void IMU_states(ofs_ekf_t* filtro, int8_t offset);

/* Calcula el jacobiano del modelo del sensor de flujo optico respecto al estado
* Pre: El filtro debe estar inicializado y haber llamado anteriormente a 'prediction_step'.
* Pos: El jacobiano tiene calculada la parte del sensor de flujo optico posicionada en las filas segun el offset. */
  void OFS_states(ofs_ekf_t* filtro, int8_t offset);

  /* Calcula el jacobiano del modelo del sensor de diistancia respecto al estado
* Pre: El filtro debe estar inicializado y haber llamado anteriormente a 'prediction_step'.
* Pos: El jacobiano tiene calculada la parte del sensor de diistancia posicionada en las filas segun el offset. */
  void TOFS_states(ofs_ekf_t* filtro, int8_t offset);

#endif /* EKF_H_ */