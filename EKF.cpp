#include "EKF.hpp"

void ofs_ekf_init(ofs_ekf_t* filtro){
    (*filtro).N =  N_STATES; // Cantidad de estados
    (*filtro).beta = false; // Indica si hay una lectura nueva del OFS
    (*filtro).gamma = false; // Indica si hay una lectura nueva del sensor de distancia
    double temp[3] = {0, 0, 9.81};
    (*filtro).g = vec2quat(temp);
    (*filtro).M00 = N_OBS_00; // Cantidad de observaciones con beta=gamma=0
    (*filtro).M01 = N_OBS_01; // Cantidad de observaciones con beta=0, gamma=1
    (*filtro).M10 = N_OBS_10; // Cantidad de observaciones con beta=1, gamma=0
    (*filtro).M11 = N_OBS_11; // Cantidad de observaciones con beta=1, gamma=1
    zeros((*filtro).states, N_STATES, 1); //p, v, q, w, ba, bw
    (*filtro).states[N_P + N_V] = 1; //q1 = 1;
    zeros(*(*filtro).cov, N_STATES, N_STATES); // Matriz de covarianza de estados
    mat_addeye(*(*filtro).cov, N_STATES);
    zeros(*(*filtro).F, N_STATES, N_STATES);
    zeros(*(*filtro).W, N_STATES, N_NOISE);
    zeros(*(*filtro).Q, N_NOISE, N_NOISE);
    (*filtro).Q[0][0] = 0.1; //uax
    (*filtro).Q[1][1] = 0.1; //uay
    (*filtro).Q[2][2] = 0.1; //uaz
    (*filtro).Q[3][3] = 0.1; //uwx
    (*filtro).Q[4][4] = 0.1; //uwy
    (*filtro).Q[5][5] = 0.1; //uwz
    (*filtro).Q[6][6] = 0.01; //ubax
    (*filtro).Q[7][7] = 0.01; //ubay
    (*filtro).Q[8][8] = 0.01; //ubaz
    (*filtro).Q[9][9] = 0.01; //ubwx
    (*filtro).Q[10][10] = 0.01; //ubwy
    (*filtro).Q[11][11] = 0.01; //ubwz
    (*filtro).Q[12][12] = 0.1; //ubsz
    (*filtro).Npix = 35; // Cantidad de píxeles
    (*filtro).FOV_OF = 4.2 * M_PI / 180; // FOV del sensor de OF
    (*filtro).f = (*filtro).Npix / (2 * tan((*filtro).FOV_OF / 2));  // Factor de conversión
};
    
void prediction_step(ofs_ekf_t* filtro, mediciones_t u){
    
    /* Auxiliares */
    int puntero = 0;
    double p[] = {(*filtro).states[puntero], (*filtro).states[puntero+1], (*filtro).states[puntero+2]};
    puntero += N_P;
    double v[] = {(*filtro).states[puntero+0], (*filtro).states[puntero+1], (*filtro).states[puntero+2]};
    puntero += N_V;
    quaternion_t q = {(*filtro).states[puntero+0], (*filtro).states[puntero+1], (*filtro).states[puntero+2], (*filtro).states[puntero+3]};
    puntero += N_Q;
    double w[] = {(*filtro).states[puntero+0], (*filtro).states[puntero+1], (*filtro).states[puntero+2]};
    quaternion_t qw = vec2quat(w);
    puntero += N_W;
    double ba[] = {(*filtro).states[puntero+0], (*filtro).states[puntero+1], (*filtro).states[puntero+2]};
    puntero += N_BA;
    double bw[] = {(*filtro).states[puntero+0], (*filtro).states[puntero+1], (*filtro).states[puntero+2]};
    quaternion_t qa_meas = {0, u.ax, u.ay, u.az};
    quaternion_t qw_meas = {0, u.wx, u.wy, u.wz};
    quaternion_t qba = vec2quat(ba);
    quaternion_t qbw = vec2quat(bw);
    quaternion_t aux;
    double aux2[N_V];
    double aux3[N_A];
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
    (*filtro).F[N_P][N_P+N_V] = 2*(u.ax - ba[0]) * q.q1 + 2*(u.ay - ba[1])*q.q4 - 2*(u.az - ba[2])*q.q3;
    (*filtro).F[N_P+1][N_P+N_V] = -2*(u.ax - ba[0]) * q.q4 + 2*(u.ay - ba[1])*q.q1 + 2*(u.az - ba[2])*q.q2;
    (*filtro).F[N_P+2][N_P+N_V] = 2*(u.ax - ba[0]) * q.q3 - 2*(u.ay - ba[1])*q.q2 + 2*(u.az - ba[2])*q.q1;
    (*filtro).F[N_P][N_P+N_V+1] = 2*(u.ax - ba[0]) * q.q2 + 2*(u.ay - ba[1])*q.q3 + 2*(u.az - ba[2])*q.q4;
    (*filtro).F[N_P+1][N_P+N_V+1] = 2*(u.ax - ba[0]) * q.q3 - 2*(u.ay - ba[1])*q.q2 + 2*(u.az - ba[2])*q.q1;
    (*filtro).F[N_P+2][N_P+N_V+1] = 2*(u.ax - ba[0]) * q.q4 - 2*(u.ay - ba[1])*q.q1 - 2*(u.az - ba[2])*q.q2;
    (*filtro).F[N_P][N_P+N_V+2] = -2*(u.ax - ba[0]) * q.q3 + 2*(u.ay - ba[1])*q.q2 - 2*(u.az - ba[2])*q.q1;
    (*filtro).F[N_P+1][N_P+N_V+2] = 2*(u.ax - ba[0]) * q.q2 + 2*(u.ay - ba[1])*q.q3 + 2*(u.az - ba[2])*q.q4;
    (*filtro).F[N_P+2][N_P+N_V+2] = 2*(u.ax - ba[0]) * q.q1 + 2*(u.ay - ba[1])*q.q4 - 2*(u.az - ba[2])*q.q3;
    (*filtro).F[N_P][N_P+N_V+3] = -2*(u.ax - ba[0]) * q.q4 + 2*(u.ay - ba[1])*q.q1 + 2*(u.az - ba[2])*q.q2;
    (*filtro).F[N_P+1][N_P+N_V+3] = -2*(u.ax - ba[0]) * q.q1 - 2*(u.ay - ba[1])*q.q4 + 2*(u.az - ba[2])*q.q3;
    (*filtro).F[N_P+2][N_P+N_V+3] = 2*(u.ax - ba[0]) * q.q2 + 2*(u.ay - ba[1])*q.q3 + 2*(u.az - ba[2])*q.q4;
    //d(velocidad) / d(ba)
    (*filtro).F[N_P][N_P+N_V+N_Q] = -pow(q.q1, 2) - pow(q.q2, 2) + pow(q.q3, 2) + pow(q.q4, 2);
    (*filtro).F[N_P+1][N_P+N_V+N_Q] = 2 * q.q1 * q.q4 -2 * q.q2 * q.q3;
    (*filtro).F[N_P+2][N_P+N_V+N_Q] = -2 * q.q1 * q.q3 -2 * q.q2 * q.q4;
    (*filtro).F[N_P][N_P+N_V+N_Q+1] = -2 * q.q1 * q.q4 -2 * q.q2 * q.q3;
    (*filtro).F[N_P+1][N_P+N_V+N_Q+1] = -pow(q.q1, 2) + pow(q.q2, 2) - pow(q.q3, 2) + pow(q.q4, 2);
    (*filtro).F[N_P+2][N_P+N_V+N_Q+1] = 2 * q.q1 * q.q2 -2 * q.q3 * q.q4;
    (*filtro).F[N_P][N_P+N_V+N_Q+2] = 2 * q.q1 * q.q3 -2 * q.q2 * q.q4;
    (*filtro).F[N_P+1][N_P+N_V+N_Q+2] = -2 * q.q1 * q.q2 -2 * q.q3 * q.q4;
    (*filtro).F[N_P+2][N_P+N_V+N_Q+2] = -pow(q.q1, 2) + pow(q.q2, 2) + pow(q.q3, 2) - pow(q.q4, 2);
    //d(q) / d(q1)
    (*filtro).F[N_P+N_V][N_P+N_V] = 1;
    (*filtro).F[N_P+N_V+1][N_P+N_V] = w[0] * u.dt / 2;
    (*filtro).F[N_P+N_V+2][N_P+N_V] = w[1] * u.dt / 2;
    (*filtro).F[N_P+N_V+3][N_P+N_V] = w[2] * u.dt / 2;
    //d(q) / d(q2)
    (*filtro).F[N_P+N_V][N_P+N_V+1] = -w[0] * u.dt / 2;
    (*filtro).F[N_P+N_V+1][N_P+N_V+1] = 1;
    (*filtro).F[N_P+N_V+2][N_P+N_V+1] = -w[2] * u.dt / 2;
    (*filtro).F[N_P+N_V+3][N_P+N_V+1] = w[1] * u.dt / 2;
    //d(q) / d(q3)
    (*filtro).F[N_P+N_V][N_P+N_V+2] = -w[1] * u.dt / 2;
    (*filtro).F[N_P+N_V+1][N_P+N_V+2] = w[2] * u.dt / 2;
    (*filtro).F[N_P+N_V+2][N_P+N_V+2] = 1;
    (*filtro).F[N_P+N_V+3][N_P+N_V+2] = -w[0] * u.dt / 2;
    //d(q) / d(q4)
    (*filtro).F[N_P+N_V][N_P+N_V+3] = -w[2] * u.dt / 2;
    (*filtro).F[N_P+N_V+1][N_P+N_V+3] = -w[1] * u.dt / 2;
    (*filtro).F[N_P+N_V+2][N_P+N_V+3] = w[0] * u.dt / 2;
    (*filtro).F[N_P+N_V+3][N_P+N_V+3] = 1;
    //d(q) / d(wx)
    (*filtro).F[N_P+N_V][N_P+N_V+N_Q] = -q.q2 * u.dt / 2;
    (*filtro).F[N_P+N_V+1][N_P+N_V+N_Q] = q.q1 * u.dt / 2;
    (*filtro).F[N_P+N_V+2][N_P+N_V+N_Q] = q.q4 * u.dt / 2;
    (*filtro).F[N_P+N_V+3][N_P+N_V+N_Q] = -q.q2 * u.dt / 2;
    //d(q) / d(wy)
    (*filtro).F[N_P+N_V][N_P+N_V+N_Q+1] = -q.q3 * u.dt / 2;
    (*filtro).F[N_P+N_V+1][N_P+N_V+N_Q+1] = -q.q4 * u.dt / 2;
    (*filtro).F[N_P+N_V+2][N_P+N_V+N_Q+1] = q.q1 * u.dt / 2;
    (*filtro).F[N_P+N_V+3][N_P+N_V+N_Q+1] = q.q2 * u.dt / 2;
    //d(q) / d(wz)
    (*filtro).F[N_P+N_V][N_P+N_V+N_Q+2] = -q.q4 * u.dt / 2;
    (*filtro).F[N_P+N_V+1][N_P+N_V+N_Q+2] = q.q3 * u.dt / 2;
    (*filtro).F[N_P+N_V+2][N_P+N_V+N_Q+2] = -q.q2 * u.dt / 2;
    (*filtro).F[N_P+N_V+3][N_P+N_V+N_Q+2] = q.q1 * u.dt / 2;
    //d(w) / d(q1)
    (*filtro).F[N_P+N_V+N_Q][N_P+N_V] = 2*q.q1*(u.wx-bw[0])+2*q.q4*(u.wy-bw[1])-2*q.q3*(u.wz-bw[2]);
    (*filtro).F[N_P+N_V+N_Q+1][N_P+N_V] = -2*q.q4*(u.wx-bw[0])+2*q.q1*(u.wy-bw[1])+2*q.q2*(u.wz-bw[2]);
    (*filtro).F[N_P+N_V+N_Q+2][N_P+N_V] = 2*q.q3*(u.wx-bw[0])-2*q.q2*(u.wy-bw[1])+2*q.q1*(u.wz-bw[2]);
    //d(w) / d(q2)
    (*filtro).F[N_P+N_V+N_Q][N_P+N_V+1] = 2*q.q2*(u.wx-bw[0])+2*q.q3*(u.wy-bw[1])+2*q.q4*(u.wz-bw[2]);
    (*filtro).F[N_P+N_V+N_Q+1][N_P+N_V+1] = 2*q.q3*(u.wx-bw[0])-2*q.q2*(u.wy-bw[1])+2*q.q1*(u.wz-bw[2]);
    (*filtro).F[N_P+N_V+N_Q+2][N_P+N_V+1] = 2*q.q4*(u.wx-bw[0])-2*q.q1*(u.wy-bw[1])-2*q.q2*(u.wz-bw[2]);
    //d(w) / d(q3)
    (*filtro).F[N_P+N_V+N_Q][N_P+N_V+2] = -2*q.q3*(u.wx-bw[0])+2*q.q2*(u.wy-bw[1])-2*q.q1*(u.wz-bw[2]);
    (*filtro).F[N_P+N_V+N_Q+1][N_P+N_V+2] = 2*q.q2*(u.wx-bw[0])+2*q.q3*(u.wy-bw[1])+2*q.q4*(u.wz-bw[2]);
    (*filtro).F[N_P+N_V+N_Q+2][N_P+N_V+2] = 2*q.q1*(u.wx-bw[0])+2*q.q4*(u.wy-bw[1])-2*q.q3*(u.wz-bw[2]);
    //d(w) / d(q4)
    (*filtro).F[N_P+N_V+N_Q][N_P+N_V+3] = -2*q.q4*(u.wx-bw[0])+2*q.q1*(u.wy-bw[1])+2*q.q2*(u.wz-bw[2]);
    (*filtro).F[N_P+N_V+N_Q+1][N_P+N_V+3] = -2*q.q1*(u.wx-bw[0])-2*q.q4*(u.wy-bw[1])+2*q.q3*(u.wz-bw[2]);
    (*filtro).F[N_P+N_V+N_Q+2][N_P+N_V+3] = 2*q.q2*(u.wx-bw[0])+2*q.q3*(u.wy-bw[1])+2*q.q4*(u.wz-bw[2]);
    //d(w) / d(bwx)
    (*filtro).F[N_P+N_V+N_Q][N_P+N_V+N_Q+N_W] = -pow(q.q1, 2) - pow(q.q2, 2) + pow(q.q3, 2) + pow(q.q4, 2);
    (*filtro).F[N_P+N_V+N_Q+1][N_P+N_V+N_Q+N_W] = 2*q.q1*q.q4 - 2*q.q2*q.q3;
    (*filtro).F[N_P+N_V+N_Q+2][N_P+N_V+N_Q+N_W] = -2*q.q1*q.q3 - 2*q.q2*q.q4;
    //d(w) / d(bwy)
    (*filtro).F[N_P+N_V+N_Q][N_P+N_V+N_Q+N_W+1] = -2*q.q1*q.q4 - 2*q.q2*q.q3;
    (*filtro).F[N_P+N_V+N_Q+1][N_P+N_V+N_Q+N_W+1] = -pow(q.q1, 2) + pow(q.q2, 2) - pow(q.q3, 2) + pow(q.q4, 2);
    (*filtro).F[N_P+N_V+N_Q+2][N_P+N_V+N_Q+N_W+1] = 2*q.q1*q.q2 - 2*q.q3*q.q4;
    //d(w) / d(bwz)
    (*filtro).F[N_P+N_V+N_Q][N_P+N_V+N_Q+N_W+2] = 2*q.q1*q.q3 - 2*q.q2*q.q4;
    (*filtro).F[N_P+N_V+N_Q+1][N_P+N_V+N_Q+N_W+2] = -2*q.q1*q.q2 - 2*q.q3*q.q4;
    (*filtro).F[N_P+N_V+N_Q+2][N_P+N_V+N_Q+N_W+2] = -pow(q.q1, 2) + pow(q.q2, 2) + pow(q.q3, 2) - pow(q.q4, 2);

    /* Wk */
    // d(velocidad) / d(uax)
    (*filtro).W[N_P][0] = (-pow(q.q1, 2) - pow(q.q2, 2) + pow(q.q3, 2) + pow(q.q4, 2)) * u.dt;
    (*filtro).W[N_P+1][0] = (2*q.q1*q.q4 - 2*q.q2*q.q3) * u.dt;
    (*filtro).W[N_P+2][0] = (-2*q.q1*q.q3 - 2*q.q2*q.q4) * u.dt;
    // d(velocidad) / d(uay)
    (*filtro).W[N_P][1] = (-2*q.q1*q.q4 - 2*q.q2*q.q3) * u.dt;
    (*filtro).W[N_P+1][1] = (-pow(q.q1, 2) + pow(q.q2, 2) - pow(q.q3, 2) + pow(q.q4, 2)) * u.dt;
    (*filtro).W[N_P+2][1] = (2*q.q1*q.q2 - 2*q.q3*q.q4) * u.dt;
    // d(velocidad) / d(uaz)
    (*filtro).W[N_P][2] = (2*q.q1*q.q3 - 2*q.q2*q.q4) * u.dt;
    (*filtro).W[N_P+1][2] = (-2*q.q1*q.q2 - 2*q.q3*q.q4) * u.dt;
    (*filtro).W[N_P+2][2] = (-pow(q.q1, 2) + pow(q.q2, 2) + pow(q.q3, 2) - pow(q.q4, 2)) * u.dt;
    // d(w) / d(uwx)
    (*filtro).W[N_P+N_V+N_Q][3] = (-pow(q.q1, 2) - pow(q.q2, 2) + pow(q.q3, 2) + pow(q.q4, 2)) * u.dt;
    (*filtro).W[N_P+N_V+N_Q+1][3] = (2*q.q1*q.q4 - 2*q.q2*q.q3) * u.dt;
    (*filtro).W[N_P+N_V+N_Q+2][3] = (-2*q.q1*q.q3 - 2*q.q2*q.q4) * u.dt;
    // d(w) / d(uwy)
    (*filtro).W[N_P+N_V+N_Q][4] = (-2*q.q1*q.q4 - 2*q.q2*q.q3) * u.dt;
    (*filtro).W[N_P+N_V+N_Q+1][4] = (-pow(q.q1, 2) + pow(q.q2, 2) - pow(q.q3, 2) + pow(q.q4, 2)) * u.dt;
    (*filtro).W[N_P+N_V+N_Q+2][4] = (2*q.q1*q.q2 - 2*q.q3*q.q4) * u.dt;
    // d(w) / d(uwz)
    (*filtro).W[N_P+N_V+N_Q][5] = (2*q.q1*q.q3 - 2*q.q2*q.q4) * u.dt;
    (*filtro).W[N_P+N_V+N_Q+1][5] = (-2*q.q1*q.q2 - 2*q.q3*q.q4) * u.dt;
    (*filtro).W[N_P+N_V+N_Q+2][5] = (-pow(q.q1, 2) + pow(q.q2, 2) + pow(q.q3, 2) - pow(q.q4, 2)) * u.dt;
    // d(bw) / d(ubw)
    (*filtro).W[N_P+N_V+N_Q+N_W][6] = u.dt;
    (*filtro).W[N_P+N_V+N_Q+N_W+1][7] = u.dt;
    (*filtro).W[N_P+N_V+N_Q+N_W+2][8] = u.dt;
    // d(ba) / d(uba)
    (*filtro).W[N_P+N_V+N_Q+N_W+N_BW][9] = u.dt;
    (*filtro).W[N_P+N_V+N_Q+N_W+N_BW+1][10] = u.dt;
    (*filtro).W[N_P+N_V+N_Q+N_W+N_BW+2][11] = u.dt;

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
    quat_sub(&qa_meas, qa_meas, qba); // qa_meas - qba
    aux = quat_mult(qa_meas, quat_conjugate(q)); // (qa_meas - qba) * -q
    aux = quat_mult(q, qa_meas); // q * (qa_meas - qba) * -q
    quat_sub(&aux, aux, (*filtro).g); // q * (qa_meas - qba) * -q - g
    quat2vec(aux, aux3);
    matmul_scalar(aux3, N_A, 1, u.dt); // dt ((a_meas - ba) - g)
    add(v, aux3, (*filtro).states + N_P, N_V, 1); // Vk+1 = Vk + dt * (a - ba - g)
    // Quaternion
    aux = quat_mult(q, qw); // qk * qw
    quat_scalar(&aux, u.dt / 2); // (dt / 2) qk * qw
    quat_add(&aux, aux, q); // qk + (dt / 2) qk * qw
    (*filtro).states[N_P + N_V] = aux.q1;
    (*filtro).states[N_P + N_V + 1] = aux.q2;
    (*filtro).states[N_P + N_V + 2] = aux.q3;
    (*filtro).states[N_P + N_V + 3] = aux.q4;
    // Velocidad angular
    quat_sub(&qw_meas, qw_meas, qbw); // qw_meas - qbw
    aux = quat_mult(qw_meas, quat_conjugate(q)); // (qw_meas - qbw) * -q
    aux = quat_mult(q, qw_meas); // q * (qw_meas - qbw) * -q
    quat2vec(aux, aux3);
    (*filtro).states[N_P + N_V + N_Q] = aux3[0];
    (*filtro).states[N_P + N_V + N_Q + 1] = aux3[1];
    (*filtro).states[N_P + N_V + N_Q + 2] = aux3[2];
    //ba y bw no cambian
};

void correction_step(ofs_ekf_t* filtro, mediciones_t z, double dt){

};