#include "EKF.hpp"

#include <iostream>
void ofs_ekf_init(ofs_ekf_t* filtro){
    (*filtro).N =  N_STATES; // Cantidad de estados
    (*filtro).beta = false; // Indica si hay una lectura nueva del OFS
    (*filtro).gamma = false; // Indica si hay una lectura nueva del sensor de distancia
    (*filtro).g.q1 = 0.0;
    (*filtro).g.q2 = 0.0;
    (*filtro).g.q3 = 0.0;
    (*filtro).g.q4 = 9.81;
    (*filtro).M00 = N_OBS_00; // Cantidad de observaciones con beta=gamma=0
    (*filtro).M01 = N_OBS_01; // Cantidad de observaciones con beta=0, gamma=1
    (*filtro).M10 = N_OBS_10; // Cantidad de observaciones con beta=1, gamma=0
    (*filtro).M11 = N_OBS_11; // Cantidad de observaciones con beta=1, gamma=1
    zeros((*filtro).states, N_STATES, 1); //p, v, q, w, ba, bw
    zeros(*(*filtro).cov, N_STATES, N_STATES); // Matriz de covarianza de estados
    (*filtro).Npix = 35; // Cantidad de píxeles
    (*filtro).FOV_OF = 4.2 * M_PI / 180; // FOV del sensor de OF
    (*filtro).f = (*filtro).Npix / (2 * tan((*filtro).FOV_OF / 2));  // Factor de conversión
};
    
void prediction_step(ofs_ekf_t* filtro, mediciones_t u, double dt){

};

void correction_step(ofs_ekf_t* filtro, mediciones_t z, double dt){

};