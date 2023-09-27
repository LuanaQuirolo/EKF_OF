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
    (*filtro).Npix = 35; // Cantidad de píxeles
    (*filtro).FOV_OF = 4.2 * M_PI / 180; // FOV del sensor de OF
    (*filtro).f = (*filtro).Npix / (2 * tan((*filtro).FOV_OF / 2));  // Factor de conversión
};
    
void prediction_step(ofs_ekf_t* filtro, mediciones_t u){
    
    int puntero = 0;
    double p[] = {(*filtro).states[puntero], (*filtro).states[puntero+1], (*filtro).states[puntero+2]};
    puntero += N_P;
    double v[] = {(*filtro).states[puntero+0], (*filtro).states[puntero+1], (*filtro).states[puntero+2]};
    puntero += N_V;
    quaternion_t q = {(*filtro).states[puntero+0], (*filtro).states[puntero+1], (*filtro).states[puntero+2], (*filtro).states[puntero+3]};
    puntero += N_Q;
    quaternion_t qw = {0, (*filtro).states[puntero+0], (*filtro).states[puntero+1], (*filtro).states[puntero+2]};
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
    
    // Posicion
    matmul_scalar2(v, aux2, N_V, 1, u.dt); //Vk = dt * Vk
    add(p, aux2, (*filtro).states, N_P); // Pk+1 = Pk + dt * Vk
    // Velocidad
    quat_sub(&qa_meas, qa_meas, qba); // qa_meas - qba
    aux = quat_mult(qa_meas, quat_conjugate(q)); // (qa_meas - qba) * -q
    aux = quat_mult(q, qa_meas); // q * (qa_meas - qba) * -q
    quat_sub(&aux, aux, (*filtro).g); // q * (qa_meas - qba) * -q - g
    quat2vec(aux, aux3);
    matmul_scalar(aux3, N_A, 1, u.dt); // dt ((a_meas - ba) - g)
    add(v, aux3, (*filtro).states + N_P, N_V); // Vk+1 = Vk + dt * (a - ba - g)
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

    puntero = 0;
    p[0] = (*filtro).states[puntero];
    p[1] = (*filtro).states[puntero+1];
    p[2] = (*filtro).states[puntero+2];
    puntero += N_P;
    v[0] = (*filtro).states[puntero+0];
    v[1] = (*filtro).states[puntero+1]; 
    v[2] = (*filtro).states[puntero+2];
    puntero += N_V;
    q.q1 = (*filtro).states[puntero+0]; 
    q.q2 = (*filtro).states[puntero+1]; 
    q.q3 = (*filtro).states[puntero+2]; 
    q.q4 = (*filtro).states[puntero+3];
    puntero += N_Q;
    double w[] = {(*filtro).states[puntero+0], (*filtro).states[puntero+1], (*filtro).states[puntero+2]};
    puntero += N_W;
    ba[0] = (*filtro).states[puntero+0];
    ba[1] = (*filtro).states[puntero+1]; 
    ba[2] = (*filtro).states[puntero+2];
    puntero += N_BA;
    bw[0] = (*filtro).states[puntero+0];
    bw[1] = (*filtro).states[puntero+1]; 
    bw[2] = (*filtro).states[puntero+2];

    //FK
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
    for (int i = 0; i < N_STATES; i++){
        for (int j = 0; j < N_STATES; j++){
            std::cout << (((*filtro).F)[i][j]) << " ";
        };
        std::cout << std::endl;
    };  
    //WK

    //Covarianza a priori
};

void correction_step(ofs_ekf_t* filtro, mediciones_t z, double dt){

};