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
#define N_OBS 9 // a, mag, of, range
#define U_P_A pow(0.157, 2) // Ruido proceso acelerometro
#define U_P_W pow(0.0167, 2) // Ruido proceso giroscopio
#define U_A pow(0.157, 2) // Ruido medicion acelerometro
#define U_FLOWX pow(0.3888, 2) // Ruido medicion flujo optico
#define U_FLOWY pow(0.6861, 2) // Ruido medicion flujo optico
#define U_RANGE pow(0.00308, 2) // Ruido medicion distancia
#define U_MAGX pow(0.01145, 2) // Ruido medicion magnetometro
#define U_MAGY pow(0.008894, 2) // Ruido medicion magnetometro
#define U_MAGZ pow(0.00585, 2) // Ruido medicion magnetometro
#define N_IMU 3 // Para el paso de correccion, en realidad son 6 mediciones
#define N_MAG 3
#define N_OFS 2
#define N_TOFS 1
#define g 9.81
#define MIN_HEIGHT 0.04
#define MAX_HEIGHT 600

typedef struct mediciones{
  double dt;
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
    double exp_meas[N_OBS]; //a, w, flow, z
    double cov[N_STATES][N_STATES]; // Matriz de covarianza de estados
    double F[N_STATES][N_STATES]; // Derivada de vector de estados respecto de si mismo
    double H[N_OBS][N_STATES]; // Derivada de modelo de medicion respecto a los estados
    double G[N_STATES][N_OBS]; // Ganancia de Kalman
    double WQWt[N_STATES][N_STATES]; // Matriz W * Q * Wt
    double R[N_OBS][N_OBS]; // Matriz ruidos mediciones
    double m[N_MAG]; //Campo magnético inicial
    double f;  // Factor de conversión
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
    double m_aux[N_PROC_NOISE][N_STATES];
    double m_aux2[N_STATES][N_STATES];
    double m_aux3[N_STATES][N_STATES];
    double m_aux4[N_STATES][N_OBS];
    double m_aux5[N_OBS][N_OBS];
    double m_aux6[N_OBS][N_OBS];
    double m_aux7[N_STATES][N_OBS];
    /* Vectores temporales */
    double v_aux[N_V];
    double v_aux2[N_V];
    double v_aux3[N_STATES];
    double v_aux4[N_OBS];
    double meas[N_OBS];
    double v_aux5[N_STATES]; //TODO BORRAR
    double v_aux6[N_STATES]; //TODO BORRAR
    double m_aux8[N_STATES][N_STATES]; //TODO BORRAR
    double S;//TODO BORRAR
    /* Cuaterniones temporales */
    quaternion_t q_aux;
    quaternion_t q_aux2;
} ofs_ekf_t; 

/* Inicializa el filtro.
* Pre: Debe haber mediciones de la IMU, el magnetometro y el sensor de distancia.
* Pos: Se crean las matrices pertinentes y los estados se inicializan con los valores dados por las mediciones. */
  void ofs_ekf_init(ofs_ekf_t* filtro, mediciones_t *meas);

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

/* Calcula la traza de la covarianza
* Pre: El filtro debe estar inicializado.
* Pos: Devuelve la traza de la matriz de covarianza. */
  double calc_trace_cov(ofs_ekf_t *filtro);
#endif /* EKF_H_ */