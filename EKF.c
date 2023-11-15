#include "EKF.h"

void ofs_ekf_init(ofs_ekf_t* filtro){
    (*filtro).N =  N_STATES; // Cantidad de estados
    (*filtro).beta = 0; // Indica si hay una lectura nueva del OFS
    (*filtro).gamma = 0; // Indica si hay una lectura nueva del sensor de distancia
    double temp[3] = {0, 0, 9.81};
    (*filtro).qg = vec2quat(temp);
    (*filtro).M00 = N_OBS_00; // Cantidad de observaciones con beta=gamma=0
    (*filtro).M01 = N_OBS_01; // Cantidad de observaciones con beta=0, gamma=1
    (*filtro).M10 = N_OBS_10; // Cantidad de observaciones con beta=1, gamma=0
    (*filtro).M11 = N_OBS_11; // Cantidad de observaciones con beta=1, gamma=1
    mat_zeros((*filtro).states, N_STATES, 1); //p, v, q
    (*filtro).states[N_P + N_V] = 1; //q1 = 1;
    mat_zeros(*(*filtro).cov, N_STATES, N_STATES); // Matriz de covarianza de estados
    mat_addeye(*(*filtro).cov, N_STATES);
    mat_zeros(*(*filtro).F, N_STATES, N_STATES);
    mat_zeros(*(*filtro).W, N_STATES, N_NOISE);
    /* Wk */
    // d(velocidad) / d(uax)
    (*filtro).W[N_P][0] = 1;
    // d(velocidad) / d(uay)
    (*filtro).W[N_P+1][1] = 1;
    // d(velocidad) / d(uaz)
    (*filtro).W[N_P+2][2] = 1;
    mat_zeros(*(*filtro).Q, N_NOISE, N_NOISE);
    (*filtro).Q[0][0] = 0.1; //uax
    (*filtro).Q[1][1] = 0.1; //uay
    (*filtro).Q[2][2] = 0.1; //uaz
    (*filtro).Q[3][3] = 0.1; //uwx
    (*filtro).Q[4][4] = 0.1; //uwy
    (*filtro).Q[5][5] = 0.1; //uwz
    (*filtro).Q[6][6] = 0.1; //ubsz
    (*filtro).Npix = 35; // Cantidad de píxeles
    (*filtro).FOV_OF = 4.2 * M_PI / 180; // FOV del sensor de OF
    (*filtro).f = (*filtro).Npix / (2 * atan2f((*filtro).FOV_OF, 2));  // Factor de conversión
};
    
void prediction_step(ofs_ekf_t* filtro, mediciones_t u){
    
    /* Auxiliares */
    (*filtro).puntero = 0;
    double p[] = {(*filtro).states[(*filtro).puntero], (*filtro).states[(*filtro).puntero+1], (*filtro).states[(*filtro).puntero+2]};
    (*filtro).puntero += N_P;
    double v[] = {(*filtro).states[(*filtro).puntero+0], (*filtro).states[(*filtro).puntero+1], (*filtro).states[(*filtro).puntero+2]};
    (*filtro).puntero += N_V;
    quaternion_t q = {(*filtro).states[(*filtro).puntero+0], (*filtro).states[(*filtro).puntero+1], (*filtro).states[(*filtro).puntero+2], (*filtro).states[(*filtro).puntero+3]};
    quaternion_t qa_meas = {0, u.ax, u.ay, u.az};
    quaternion_t qw_meas = {0, u.wx, u.wy, u.wz};
    quaternion_t aux;
    double aux2[N_V];
    double aux3[N_V];
    double aux4[N_STATES][N_STATES];
    double aux5[N_STATES][N_STATES];
    double aux6[N_STATES][N_STATES];
    double aux7[N_STATES][N_STATES];
    double Wt[N_NOISE][N_STATES];
    double Ft[N_STATES][N_STATES];

    /* Fk */
    mat_addeye(*(*filtro).F, N_STATES);
    //d(posicion) / d(velocidad)
    (*filtro).F[0][N_P] = u.dt;
    (*filtro).F[1][N_P + 1] = u.dt;
    (*filtro).F[2][N_P + 2] = u.dt;
    //d(velocidad) / d(quaternion)
    (*filtro).F[N_P][N_P+N_V] = u.dt * (2*(u.ax) * q.q1 - 2*(u.ay)*q.q4 + 2*(u.az)*q.q3);
    (*filtro).F[N_P+1][N_P+N_V] = u.dt * (2*(u.ax) * q.q4 + 2*(u.ay)*q.q1 - 2*(u.az)*q.q2);
    (*filtro).F[N_P+2][N_P+N_V] = u.dt * (-2*(u.ax) * q.q3 + 2*(u.ay)*q.q2 + 2*(u.az)*q.q1);
    (*filtro).F[N_P][N_P+N_V+1] = u.dt * (2*(u.ax) * q.q2 + 2*(u.ay)*q.q3 + 2*(u.az)*q.q4);
    (*filtro).F[N_P+1][N_P+N_V+1] = u.dt * (2*(u.ax) * q.q3 - 2*(u.ay)*q.q2 - 2*(u.az)*q.q1);
    (*filtro).F[N_P+2][N_P+N_V+1] = u.dt * (2*(u.ax) * q.q4 + 2*(u.ay)*q.q1 - 2*(u.az)*q.q2);
    (*filtro).F[N_P][N_P+N_V+2] = u.dt * (-2*(u.ax) * q.q3 + 2*(u.ay)*q.q2 + 2*(u.az)*q.q1);
    (*filtro).F[N_P+1][N_P+N_V+2] = u.dt * (2*(u.ax) * q.q2 + 2*(u.ay)*q.q3 + 2*(u.az)*q.q4);
    (*filtro).F[N_P+2][N_P+N_V+2] = u.dt * (-2*(u.ax) * q.q1 + 2*(u.ay)*q.q4 - 2*(u.az)*q.q3);
    (*filtro).F[N_P][N_P+N_V+3] = u.dt * (-2*(u.ax) * q.q4 - 2*(u.ay)*q.q1 + 2*(u.az)*q.q2);
    (*filtro).F[N_P+1][N_P+N_V+3] = u.dt * (2*(u.ax) * q.q1 - 2*(u.ay)*q.q4 + 2*(u.az)*q.q3);
    (*filtro).F[N_P+2][N_P+N_V+3] = u.dt * (2*(u.ax) * q.q2 + 2*(u.ay)*q.q3 + 2*(u.az)*q.q4);
    //d(q) / d(q1)
    (*filtro).F[N_P+N_V][N_P+N_V] = 1;
    (*filtro).F[N_P+N_V+1][N_P+N_V] = -u.wx * u.dt / 2;
    (*filtro).F[N_P+N_V+2][N_P+N_V] = -u.wy * u.dt / 2;
    (*filtro).F[N_P+N_V+3][N_P+N_V] = -u.wz * u.dt / 2;
    //d(q) / d(q2)
    (*filtro).F[N_P+N_V][N_P+N_V+1] = u.wx * u.dt / 2;
    (*filtro).F[N_P+N_V+1][N_P+N_V+1] = 1;
    (*filtro).F[N_P+N_V+2][N_P+N_V+1] = u.wz * u.dt / 2;
    (*filtro).F[N_P+N_V+3][N_P+N_V+1] = -u.wy * u.dt / 2;
    //d(q) / d(q3)
    (*filtro).F[N_P+N_V][N_P+N_V+2] = u.wy * u.dt / 2;
    (*filtro).F[N_P+N_V+1][N_P+N_V+2] = -u.wz * u.dt / 2;
    (*filtro).F[N_P+N_V+2][N_P+N_V+2] = 1;
    (*filtro).F[N_P+N_V+3][N_P+N_V+2] = u.wx * u.dt / 2;
    //d(q) / d(q4)
    (*filtro).F[N_P+N_V][N_P+N_V+3] = u.wz * u.dt / 2;
    (*filtro).F[N_P+N_V+1][N_P+N_V+3] = u.wy * u.dt / 2;
    (*filtro).F[N_P+N_V+2][N_P+N_V+3] = -u.wx * u.dt / 2;
    (*filtro).F[N_P+N_V+3][N_P+N_V+3] = 1;

    /* Covarianza a priori */
    transpose(*(*filtro).W, *Wt, N_STATES, N_NOISE);
    mulmat(*(*filtro).Q, *Wt, *aux4, N_NOISE, N_NOISE, N_STATES);
    mulmat(*(*filtro).W, *aux4, *aux6, N_STATES, N_STATES, N_NOISE);
    transpose(*(*filtro).F, *Ft, N_STATES, N_STATES);
    mulmat(*(*filtro).cov, *Ft, *aux5, N_STATES, N_STATES, N_STATES);
    mulmat(*(*filtro).F, *aux5, *aux7, N_STATES, N_STATES, N_STATES);
    add(*aux6, *aux7, *(*filtro).cov, N_STATES, N_STATES); 

    /* Paso de prediccion */
    // Posicion
    matmul_scalar2(v, aux2, N_V, 1, u.dt); //Vk = dt * Vk
    add(p, aux2, (*filtro).states, N_P, 1); // Pk+1 = Pk + dt * Vk
    // Velocidad
    // Queremos pasar la acel de cuerpo a mundo
    aux = quat_mult(qa_meas, quat_conjugate(q)); // (qa_meas) * -q
    aux = quat_mult(q, aux); // q * (qa_meas) * -q
    quat_sub(&aux, aux, (*filtro).qg); // q * (qa_meas) * -q - g
    quat2vec(aux, aux3);
    matmul_scalar(aux3, N_V, 1, u.dt); // dt (a_meas - g)
    add(v, aux3, (*filtro).states + N_P, N_V, 1); // Vk+1 = Vk + dt * (a - g)
    // Quaternion
    aux = quat_mult(qw_meas, q); // qw * qk
    quat_scalar(&aux, u.dt / 2); // (dt / 2) qw * qk
    quat_add(&aux, aux, q); // qk + (dt / 2) qw * qk
    quat_Normalization(&aux);
    (*filtro).states[N_P + N_V] = aux.q1;
    (*filtro).states[N_P + N_V + 1] = aux.q2;
    (*filtro).states[N_P + N_V + 2] = aux.q3;
    (*filtro).states[N_P + N_V + 3] = aux.q4;
};

void correction_step(ofs_ekf_t* filtro, mediciones_t z, double dt){


if((*filtro).beta == 1 && (*filtro).gamma == 1){

}
(*filtro).beta = 0;
(*filtro).gamma = 0;
};