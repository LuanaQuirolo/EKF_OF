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

#define N_P 3
#define N_V 3
#define N_Q 4
#define N_STATES N_P+N_V+N_Q // La cantidad de estados total es la suma de sus componentes
#define N_PROC_NOISE 7 // a, w
#define N_CORR_NOISE 9 // a, mag, of, range
#define N_OBS_00 6
#define N_OBS_01 8
#define N_OBS_10 7
#define N_OBS_11 9 // a, mag, of, range
#define U_P_A pow(0.01, 2) // Ruido proceso acelerometro
#define U_P_W pow(0.001, 2) // Ruido proceso giroscopio
#define U_A pow(0.01, 2) // Ruido medicion acelerometro
#define U_FLOW pow(0.2, 2) // Ruido medicion flujo optico
#define U_RANGE pow(0.0001, 2) // Ruido medicion distancia
#define U_MAG pow(0.0001, 2) // Ruido medicion magnetometro
#define N_IMU 3 // Para el paso de correccion, en realidad son 6 mediciones
#define N_MAG 3
#define N_OFS 2
#define N_TOFS 1
#define g 9.81
#define MIN_HEIGHT 0.04
#define MAX_HEIGHT 100

typedef struct mediciones{
  double dt;
  double tau;
  double ax;
  double ay;
  double az;
  double wx;
  double wy;
  double wz;
  double mx;
  double my;
  double mz;
  double ofx;
  double ofy;
  double range;
  } mediciones_t;

/***************************************************************
 * 				   DEFINICION ESTRUCTURA OFS_EKF
 ***************************************************************/
typedef struct ofs_ekf {
    double states[N_STATES]; //p, v, q
    double exp_meas[N_OBS_11]; //a, w, flow, z
    double cov[N_STATES][N_STATES]; // Matriz de covarianza de estados
    double F[N_STATES][N_STATES]; // Derivada de vector de estados respecto de si mismo
    double H[N_OBS_11][N_STATES]; // Derivada de modelo de medicion respecto a los estados
    double G[N_STATES][N_OBS_11]; // Ganancia de Kalman
    double WQWt[N_STATES][N_STATES]; // Matriz W * Q * Wt
    double R[N_OBS_11][N_OBS_11]; // Matriz ruidos mediciones
    double m[N_MAG]; //Campo magnético inicial
    double f;  // Factor de conversión
    int beta; // Indica si hay una lectura nueva del OFS
    int gamma; // Indica si hay una lectura nueva del sensor de distancia
    /* Variables auxiliares */
    quaternion_t q; // Cuaternion para pasar de mundo a cuerpo, es parte de states
    double p[N_P]; // Posicion, es parte de states
    double v[N_V]; // Velocidad, es parte de states
    quaternion_t qg;
    quaternion_t qm; 
    quaternion_t qa_meas;
    quaternion_t qw_meas;
    quaternion_t qm_meas;
    /* Matrices auxiliares */
    double aux4[N_PROC_NOISE][N_STATES];
    double aux5[N_STATES][N_STATES];
    double aux6[N_STATES][N_STATES];
    double aux7[N_STATES][N_STATES];
    double aux9[N_STATES][N_OBS_11];
    double aux10[N_OBS_11][N_OBS_11];
    double aux11[N_OBS_11][N_OBS_11];
    double Wt[N_PROC_NOISE][N_STATES];
    double Ft[N_STATES][N_STATES];
    double Ht[N_STATES][N_OBS_11];
    /* Vectores temporales */
    double aux2[N_V];
    double aux3[N_V];
    double aux8[N_STATES];
    double aux12[N_OBS_11];
    double meas[N_OBS_11];
    /* Cuaterniones temporales */
    quaternion_t aux;
    quaternion_t aux13;
    /* Escalares auxiliares */
    uint8_t meas_counter;
} ofs_ekf_t; 

  void ofs_ekf_init(ofs_ekf_t* filtro, double m[N_MAG], double a[N_IMU], double z);

/* Realiza el paso de prediccion. Es decir, propaga la accion de control.
* Pre: El filtro debe estar inicializado.
* Pos: El estado se encuentra actualizado segun el control aplicado. */
  void prediction_step(ofs_ekf_t* filtro, mediciones_t *u);

/* Realiza el paso de correccion luego del paso de prediccion. Es decir,
* una vez propagada la accion de control, se corrige el estado en base a
* las mediciones disponibles.
* Pre: El filtro debe estar inicializado y haber llamado anteriormente a 'prediction_step'.
* Pos: El estado se encuentra corregido en base a las mediciones. */
  void correction_step(ofs_ekf_t *filtro, mediciones_t *z);

/* Calcula el jacobiano del modelo de la IMU respecto al estado
* Pre: El filtro debe estar inicializado y haber llamado anteriormente a 'prediction_step'.
* Pos: El jacobiano tiene calculada la parte de la IMU posicionada en las filas segun el offset. */
  void IMU_states(ofs_ekf_t *filtro);

/* Calcula el jacobiano del modelo del magnetometro respecto al estado
* Pre: El filtro debe estar inicializado y haber llamado anteriormente a 'prediction_step'.
* Pos: El jacobiano tiene calculada la parte del MAG posicionada en las filas segun el offset. */
 void MAG_states(ofs_ekf_t* filtro);

/* Calcula el jacobiano del modelo del sensor de flujo optico respecto al estado
* Pre: El filtro debe estar inicializado y haber llamado anteriormente a 'prediction_step'.
* Pos: El jacobiano tiene calculada la parte del sensor de flujo optico posicionada en las filas segun el offset. */
  void OFS_states(ofs_ekf_t* filtro, mediciones_t *z);

  /* Calcula el jacobiano del modelo del sensor de distancia respecto al estado
* Pre: El filtro debe estar inicializado y haber llamado anteriormente a 'prediction_step'.
* Pos: El jacobiano tiene calculada la parte del sensor de diistancia posicionada en las filas segun el offset. */
  void TOFS_states(ofs_ekf_t* filtro, int8_t offset);

/* Calcula la traza de la covarianza
* Pre: El filtro debe estar inicializado.
* Pos: Devuelve la traza de la matriz de covarianza. */
  double calc_trace_cov(ofs_ekf_t *filtro);
#endif /* EKF_H_ */