#include "EKF.hpp"

void ofs_ekf_init(ofs_ekf_t* filtro){
    (*filtro).N =  N_STATES; // Cantidad de estados
    (*filtro).beta = false; // Indica si hay una lectura nueva del OFS
    (*filtro).gamma = false; // Indica si hay una lectura nueva del sensor de distancia
    double temp[3] = {0, 0, 9.81};
    (*filtro).g = quat_extend(temp);
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
    
    double p_old[] = {(*filtro).states[0], (*filtro).states[1], (*filtro).states[2]};
    double v_old[] = {(*filtro).states[N_P+0], (*filtro).states[N_P+1], (*filtro).states[N_P+2]};
    add(p_old, v_old, (*filtro).states, N_P); // Pk+1 = Pk + dt * Vk
};

void correction_step(ofs_ekf_t* filtro, mediciones_t z, double dt){

};