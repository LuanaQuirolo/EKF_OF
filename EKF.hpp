/***************************************************************
 *            Algoritmo de fusion de sensores
 *                  Creado: 11 sep. 2023
 *                  Autor: Luana Quirolo
 *                  Padron: 102102
****************************************************************/

#ifndef EKF_H_
#define EKF_H_

#include <cmath>
#define N_STATES 18
#define N_OBS_00 6
#define N_OBS_01 8
#define N_OBS_10 7
#define N_OBS_11 9

typedef struct mediciones{
  double dt;
  double ax = 0;
  double ay = 0;
  double az = 0;
  double wx = 0;
  double wy = 0;
  double wz = 0;
  double ofx = 0;
  double ofy = 0;
  double range = 0;
  } mediciones_t;

/***************************************************************
 * 				   DEFINICION CLASE OFS_EKF
 ***************************************************************/
class OFS_EKF {
private:
    uint8_t N =  N_STATES; // Cantidad de estados
    bool beta; // Indica si hay una lectura nueva del OFS
    bool gamma; // Indica si hay una lectura nueva del sensor de distancia
    Vector3d g; 
    uint8_t M00 = N_OBS_00; // Cantidad de observaciones con beta=gamma=0
    uint8_t M01 = N_OBS_01; // Cantidad de observaciones con beta=0, gamma=1
    uint8_t M10 = N_OBS_10; // Cantidad de observaciones con beta=1, gamma=0
    uint8_t M11 = N_OBS_11; // Cantidad de observaciones con beta=1, gamma=1
    Vector<double, N_STATES> states; //p, v, angles, w, ba, bw
    Matrix<double, N_STATES, N_STATES> cov; // Matriz de covarianza de estados
    uint8_t Npix = 35; // Cantidad de píxeles
    float FOV_OF = 42 * M_PI / 180; // FOV del sensor de OF
    float f = this->Npix / (2 * sin(this->FOV_OF / 2));  // Factor de conversión

public:
/* Crea una instancia de la clase.
* Pre: 
* Pos: Devuelve el objeto con sus atributos inicializados. */
    OFS_EKF(uint8_t N); 

/* Devuelve el vector de estados del filtro.
* Pre: Debe existir una instancia de la clase.*/
  Vector<double, N_STATES> get_states(); 

/* Devuelve la matriz de covarianza de estados del filtro.
* Pre: Debe existir una instancia de la clase.*/
  Matrix<double, N_STATES, N_STATES> get_cov(); 

/* Realiza el paso de prediccion. Es decir, propaga la accion de control.
* Pre: Debe existir una instancia de la clase.
* Pos: El estado se encuentra actualizado segun el control aplicado. */
    void prediction_step(mediciones_t u, double dt);

/* Realiza el paso de correccion luego del paso de prediccion. Es decir,
* una vez propagada la accion de control, se corrige el estado en base a
* las mediciones disponibles.
* Pre: Debe existir una instancia de la clase y haber llamado anteriormente a 'prediction_step'.
* Pos: El estado se encuentra corregido en base a las mediciones. */
    void correction_step(mediciones_t z, double dt);
}; 

#endif /* EKF_H_ */