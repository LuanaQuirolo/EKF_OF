#include "EKF.h"

void ofs_ekf_init(ofs_ekf_t* filtro){
    filtro->N =  N_STATES; // Cantidad de estados
    filtro->beta = 0; // Indica si hay una lectura nueva del OFS
    filtro->gamma = 0; // Indica si hay una lectura nueva del sensor de distancia
    double temp[3] = {0, 0, 9.81};
    filtro->qg = vec2quat(temp);
    filtro->M00 = N_OBS_00; // Cantidad de observaciones con beta=gamma=0
    filtro->M01 = N_OBS_01; // Cantidad de observaciones con beta=0, gamma=1
    filtro->M10 = N_OBS_10; // Cantidad de observaciones con beta=1, gamma=0
    filtro->M11 = N_OBS_11; // Cantidad de observaciones con beta=1, gamma=1
    mat_zeros(filtro->states, N_STATES, 1); //p, v, q
    filtro->states[N_P + N_V] = 1; //q1 = 1;
    mat_zeros(*filtro->cov, N_STATES, N_STATES); // Matriz de covarianza de estados
    mat_addeye(*filtro->cov, N_STATES);
    mat_zeros(*filtro->F, N_STATES, N_STATES);
    mat_zeros(*filtro->W, N_STATES, N_PROC_NOISE);
    /* Wk */
    // d(velocidad) / d(uax)
    filtro->W[N_P][0] = 1;
    // d(velocidad) / d(uay)
    filtro->W[N_P+1][1] = 1;
    // d(velocidad) / d(uaz)
    filtro->W[N_P+2][2] = 1;
    mat_zeros(*filtro->Q, N_PROC_NOISE, N_PROC_NOISE);
    filtro->Q[0][0] = 0.1; //uax
    filtro->Q[1][1] = 0.1; //uay
    filtro->Q[2][2] = 0.1; //uaz
    filtro->Q[3][3] = 0.1; //uwx
    filtro->Q[4][4] = 0.1; //uwy
    filtro->Q[5][5] = 0.1; //uwz
    filtro->Q[6][6] = 0.1; //ubsz
    filtro->R[0][0] = 0.1; //uax
    filtro->R[1][1] = 0.1; //uay
    filtro->R[2][2] = 0.1; //uaz
    filtro->Npix = 35; // Cantidad de píxeles
    filtro->FOV_OF = 4.2 * M_PI / 180; // FOV del sensor de OF
    filtro->f = filtro->Npix / (2 * atan2f(filtro->FOV_OF, 2));  // Factor de conversión
};
    
void prediction_step(ofs_ekf_t* filtro, mediciones_t u){
    
    /* Auxiliares */
    filtro->puntero = 0;
    double p[] = {filtro->states[filtro->puntero], filtro->states[filtro->puntero+1], filtro->states[filtro->puntero+2]};
    filtro->puntero += N_P;
    double v[] = {filtro->states[filtro->puntero+0], filtro->states[filtro->puntero+1], filtro->states[filtro->puntero+2]};
    filtro->puntero += N_V;
    quaternion_t q = {filtro->states[filtro->puntero+0], filtro->states[filtro->puntero+1], filtro->states[filtro->puntero+2], filtro->states[filtro->puntero+3]};
    quaternion_t qa_meas = {0, u.ax, u.ay, u.az};
    quaternion_t qw_meas = {0, u.wx, u.wy, u.wz};
    double aux2[N_V];
    double aux3[N_V];
    double aux4[N_STATES][N_STATES];
    double aux5[N_STATES][N_STATES];
    double aux6[N_STATES][N_STATES];
    double aux7[N_STATES][N_STATES];
    double Wt[N_PROC_NOISE][N_STATES];
    double Ft[N_STATES][N_STATES];

    /* Fk */
    mat_addeye(*filtro->F, N_STATES);
    //d(posicion) / d(velocidad)
    filtro->F[0][N_P] = u.dt;
    filtro->F[1][N_P + 1] = u.dt;
    filtro->F[2][N_P + 2] = u.dt;
    //d(velocidad) / d(quaternion)
    filtro->F[N_P][N_P+N_V] = u.dt * (2*(u.ax) * q.q1 - 2*(u.ay)*q.q4 + 2*(u.az)*q.q3);
    filtro->F[N_P+1][N_P+N_V] = u.dt * (2*(u.ax) * q.q4 + 2*(u.ay)*q.q1 - 2*(u.az)*q.q2);
    filtro->F[N_P+2][N_P+N_V] = u.dt * (-2*(u.ax) * q.q3 + 2*(u.ay)*q.q2 + 2*(u.az)*q.q1);
    filtro->F[N_P][N_P+N_V+1] = u.dt * (2*(u.ax) * q.q2 + 2*(u.ay)*q.q3 + 2*(u.az)*q.q4);
    filtro->F[N_P+1][N_P+N_V+1] = u.dt * (2*(u.ax) * q.q3 - 2*(u.ay)*q.q2 - 2*(u.az)*q.q1);
    filtro->F[N_P+2][N_P+N_V+1] = u.dt * (2*(u.ax) * q.q4 + 2*(u.ay)*q.q1 - 2*(u.az)*q.q2);
    filtro->F[N_P][N_P+N_V+2] = u.dt * (-2*(u.ax) * q.q3 + 2*(u.ay)*q.q2 + 2*(u.az)*q.q1);
    filtro->F[N_P+1][N_P+N_V+2] = u.dt * (2*(u.ax) * q.q2 + 2*(u.ay)*q.q3 + 2*(u.az)*q.q4);
    filtro->F[N_P+2][N_P+N_V+2] = u.dt * (-2*(u.ax) * q.q1 + 2*(u.ay)*q.q4 - 2*(u.az)*q.q3);
    filtro->F[N_P][N_P+N_V+3] = u.dt * (-2*(u.ax) * q.q4 - 2*(u.ay)*q.q1 + 2*(u.az)*q.q2);
    filtro->F[N_P+1][N_P+N_V+3] = u.dt * (2*(u.ax) * q.q1 - 2*(u.ay)*q.q4 + 2*(u.az)*q.q3);
    filtro->F[N_P+2][N_P+N_V+3] = u.dt * (2*(u.ax) * q.q2 + 2*(u.ay)*q.q3 + 2*(u.az)*q.q4);
    //d(q) / d(q1)
    filtro->F[N_P+N_V][N_P+N_V] = 1;
    filtro->F[N_P+N_V+1][N_P+N_V] = -u.wx * u.dt / 2;
    filtro->F[N_P+N_V+2][N_P+N_V] = -u.wy * u.dt / 2;
    filtro->F[N_P+N_V+3][N_P+N_V] = -u.wz * u.dt / 2;
    //d(q) / d(q2)
    filtro->F[N_P+N_V][N_P+N_V+1] = u.wx * u.dt / 2;
    filtro->F[N_P+N_V+1][N_P+N_V+1] = 1;
    filtro->F[N_P+N_V+2][N_P+N_V+1] = u.wz * u.dt / 2;
    filtro->F[N_P+N_V+3][N_P+N_V+1] = -u.wy * u.dt / 2;
    //d(q) / d(q3)
    filtro->F[N_P+N_V][N_P+N_V+2] = u.wy * u.dt / 2;
    filtro->F[N_P+N_V+1][N_P+N_V+2] = -u.wz * u.dt / 2;
    filtro->F[N_P+N_V+2][N_P+N_V+2] = 1;
    filtro->F[N_P+N_V+3][N_P+N_V+2] = u.wx * u.dt / 2;
    //d(q) / d(q4)
    filtro->F[N_P+N_V][N_P+N_V+3] = u.wz * u.dt / 2;
    filtro->F[N_P+N_V+1][N_P+N_V+3] = u.wy * u.dt / 2;
    filtro->F[N_P+N_V+2][N_P+N_V+3] = -u.wx * u.dt / 2;
    filtro->F[N_P+N_V+3][N_P+N_V+3] = 1;

    /* Covarianza a priori */
    transpose(*filtro->W, *Wt, N_STATES, N_PROC_NOISE);
    mulmat(*filtro->Q, *Wt, *aux4, N_PROC_NOISE, N_PROC_NOISE, N_STATES);
    mulmat(*filtro->W, *aux4, *aux6, N_STATES, N_STATES, N_PROC_NOISE);
    transpose(*filtro->F, *Ft, N_STATES, N_STATES);
    mulmat(*filtro->cov, *Ft, *aux5, N_STATES, N_STATES, N_STATES);
    mulmat(*filtro->F, *aux5, *aux7, N_STATES, N_STATES, N_STATES);
    add(*aux6, *aux7, *filtro->cov, N_STATES, N_STATES); 

    /* Paso de prediccion */
    // Posicion
    matmul_scalar2(v, aux2, N_V, 1, u.dt); //Vk = dt * Vk
    add(p, aux2, filtro->states, N_P, 1); // Pk+1 = Pk + dt * Vk
    // Velocidad
    // Queremos pasar la acel de cuerpo a mundo
    filtro->aux = quat_mult(qa_meas, quat_conjugate(q)); // (qa_meas) * -q
    filtro->aux = quat_mult(q, filtro->aux); // q * (qa_meas) * -q
    quat_sub(&filtro->aux, filtro->aux, filtro->qg); // q * (qa_meas) * -q - g
    quat2vec(filtro->aux, aux3);
    matmul_scalar(aux3, N_V, 1, u.dt); // dt (a_meas - g)
    add(v, aux3, filtro->states + N_P, N_V, 1); // Vk+1 = Vk + dt * (a - g)
    // Quaternion
    filtro->aux = quat_mult(qw_meas, q); // qw * qk
    quat_scalar(&filtro->aux, u.dt / 2); // (dt / 2) qw * qk
    quat_add(&filtro->aux, filtro->aux, q); // qk + (dt / 2) qw * qk
    quat_Normalization(&filtro->aux);
    filtro->states[N_P + N_V] = filtro->aux.q1;
    filtro->states[N_P + N_V + 1] = filtro->aux.q2;
    filtro->states[N_P + N_V + 2] = filtro->aux.q3;
    filtro->states[N_P + N_V + 3] = filtro->aux.q4;
};

void correction_step(ofs_ekf_t* filtro, mediciones_t* z){

mat_zeros(*filtro->H, N_CORR_NOISE, N_STATES);
uint8_t meas_counter;
double meas[N_OBS_11];
/*************************** Jacobiano **************************/
if(filtro->beta == 1 && filtro->gamma == 1){
    IMU_states(filtro, 0);
    OFS_states(filtro, N_IMU, z);
    TOFS_states(filtro, N_IMU + N_OFS);
    meas_counter = N_OBS_11;
}
else if (filtro->beta == 0 && filtro->gamma == 1){
    IMU_states(filtro, 0);
    TOFS_states(filtro, N_IMU);
    meas_counter = N_OBS_10;
}
else if (filtro->beta == 1 && filtro->gamma == 0){
    IMU_states(filtro, 0);
    OFS_states(filtro, N_IMU, z);
    meas_counter = N_OBS_10;
}
else {
    IMU_states(filtro, 0);
    meas_counter = N_OBS_00;
}

/*************************** Mediciones predichas **************************/
mat_zeros(filtro->measurements, N_OBS_11, 1);
quaternion_t q = {filtro->states[N_P + N_V+0], filtro->states[N_P + N_V+1], \
                 filtro->states[N_P + N_V+2], filtro->states[N_P + N_V+3]};
// Aceleracion
filtro->aux = quat_mult(filtro->qg, q); // qg * qk (de mundo a cuerpo)
filtro->aux = quat_mult(quat_conjugate(q), filtro->aux);
quat2vec(filtro->aux, filtro->measurements);
meas[0] = z->ax;
meas[1] = z->ay;
meas[2] = z->az;
// Flujo optico (Si corresponde)
if(filtro->beta == 1){
    filtro->aux = vec2quat(filtro->states + 3); //Paso la velocidad respecto al mundo a quat para convertirlo a body
    filtro->aux = quat_mult(filtro->qg, q);
    filtro->aux = quat_mult(quat_conjugate(q), filtro->aux);
    double aux2[N_V];
    quat2vec(filtro->aux, aux2);
    filtro->measurements[N_OBS_00 + 0] = -(*z).dt * filtro->f * ((2 * (pow(q.q1, 2) + pow(q.q4, 2)) - 1) / \
                                         filtro->states[2] * aux2[0] + (*z).wy);
    filtro->measurements[N_OBS_00 + 1] = -(*z).dt * filtro->f * ((2 * (pow(q.q1, 2) + pow(q.q4, 2)) - 1) / \
                                         filtro->states[2] * aux2[1] - (*z).wx);  
    filtro->R[N_OBS_00][N_OBS_00] = 0.1; //u flowx    
    filtro->R[N_OBS_00 + 1][N_OBS_00 + 1] = 0.1; //u flowy  
    meas[N_OBS_00 + 0] = z->ofx;
    meas[N_OBS_00 + 1] = z->ofy;                                                                
}
// Sensor de distancia (Si corresponde)
if(filtro->beta == 1 && filtro->gamma == 1){
    filtro->measurements[N_OBS_11 - 1] = filtro->states[2] / (2 * (pow(q.q1, 2) + pow(q.q4, 2)) - 1);
    filtro->R[N_OBS_11 - 1][N_OBS_11 - 1] = 0.1; // distance
    meas[N_OBS_11 - 1] = z->range; 
}
else if(filtro->beta == 0 && filtro->gamma == 1){
    filtro->measurements[N_OBS_10 - 1] = filtro->states[2] / (2 * (pow(q.q1, 2) + pow(q.q4, 2)) - 1);
    filtro->R[N_OBS_10 - 1][N_OBS_10 - 1] = 0.1; // distance  
}
/*************************** Ganancia de Kalman **************************/
double Ht[N_STATES][meas_counter];
double aux3[N_STATES][meas_counter];
double aux4[meas_counter][meas_counter];
double aux5[meas_counter][meas_counter];
double aux6[meas_counter];
double aux7[N_STATES][N_STATES];
transpose(*filtro->H, *Ht, meas_counter, N_STATES); 
mulmat(*filtro->cov, *Ht, *aux3, N_STATES, N_STATES, meas_counter); // aux3 = cov * Ht
mulmat(*filtro->H, *aux3, *aux4, meas_counter, N_STATES, meas_counter); // aux4 = H * cov * Ht
accum(*aux4, *filtro->R, meas_counter, meas_counter); // aux4 = H * cov * Ht + R
cholsl(*aux4, *aux5, aux6, meas_counter); // aux5 = inv(H * cov * Ht + R)
mulmat(*aux3, *aux5, *filtro->G, N_STATES, meas_counter, meas_counter); // G = cov * Ht * inv(H * cov * Ht + R)

/*************************** Correccion **************************/
sub(meas, filtro->measurements, aux6, meas_counter); // aux6 = z_medido - z_esperado
mulvec(*filtro->G, aux6, meas, N_STATES, meas_counter); // meas = G(z_medido - z_esperado)
accum(filtro->states, aux6, N_STATES, 1); // estado = estado + G(z_medido - z_esperado)

/*************************** Covarianza **************************/
mulmat(*filtro->G, *filtro->H, *aux3, N_STATES, meas_counter, N_STATES); //  aux3 = G * H
mulmat(*aux3, *filtro->cov, *aux7, N_STATES, N_STATES, N_STATES); // aux7 = G * H * cov_priori
mat_negate(*aux7, N_STATES, N_STATES); // aux7 = - G * H * cov_priori
accum(*filtro->cov, *aux7, N_STATES, N_STATES); // cov = cov_priori - G * H * cov_priori
filtro->beta = 0;
filtro->gamma = 0;
}

void IMU_states(ofs_ekf_t* filtro, int8_t offset){

// partial_ax / partial_q
filtro->H[offset + 0][6] = -2 * g * (filtro->states[N_P+N_V+2]);     
filtro->H[offset + 0][7] = 2 * g * (filtro->states[N_P+N_V+3]);     
filtro->H[offset + 0][8] = -2 * g * (filtro->states[N_P+N_V]);     
filtro->H[offset + 0][9] = 2 * g * (filtro->states[N_P+N_V+1]);  
// partial_ay / partial_q   
filtro->H[offset + 1][6] = 2 * g * (filtro->states[N_P+N_V+1]);      
filtro->H[offset + 1][7] = 2 * g * (filtro->states[N_P+N_V]);      
filtro->H[offset + 1][8] = 2 * g * (filtro->states[N_P+N_V+3]);     
filtro->H[offset + 1][9] = 2 * g * (filtro->states[N_P+N_V+2]);
// partial_az / partial_q  
filtro->H[offset + 2][6] = 2 * g * (filtro->states[N_P+N_V]);      
filtro->H[offset + 2][7] = -2 * g * (filtro->states[N_P+N_V+1]);      
filtro->H[offset + 2][8] = -2 * g * (filtro->states[N_P+N_V+2]);     
filtro->H[offset + 2][9] = 2 * g * (filtro->states[N_P+N_V+3]);
}

void OFS_states(ofs_ekf_t* filtro, int8_t offset, mediciones_t *z){

double pz = filtro->states[2];
double vx = filtro->states[N_P];
double vy = filtro->states[N_P + 1];
double vz = filtro->states[N_P + 2];
quaternion_t q = {filtro->states[N_P + N_V + 0], filtro->states[N_P + N_V + 1], filtro->states[N_P + N_V + 2], filtro->states[N_P + N_V + 3]};
// partial_nx / partial_pz
filtro->H[offset][2] = ((*z).dt * filtro->f * (2 * pow(q.q1, 2) + 2 * pow(q.q4, 2) - 1)) * \
                        (q.q1 * (q.q1 * vx - q.q3 * vz + q.q4 * vy) - q.q2 * (-q.q2 * vx - q.q3 * vy - q.q4 * vz) \
                        -q.q3 * (q.q1 * vz - q.q2 * vy + q.q3 * vx) + q.q4 * (q.q1 * vy + q.q2 * vz - q.q2 * vx) \
                        / pow(pz, 2));     
// partial_nx / partial_v   
filtro->H[offset][3] = - ((*z).dt * filtro->f * (2 * pow(q.q1, 2) + 2 * pow(q.q4, 2) - 1)) * \
                         (pow(q.q1, 2) + pow(q.q2, 2) + pow(q.q3, 2) + pow(q.q4, 2)) / pz;       
filtro->H[offset][4] = - ((*z).dt * filtro->f * (2 * q.q1 * q.q4 + 2 * q.q2 * q.q3)) * \
                         (2 * pow(q.q1, 2) + 2 * pow(q.q4, 2) - 1) / pz;        
filtro->H[offset][5] = - ((*z).dt * filtro->f * (-2 * q.q1 * q.q3 + 2 * q.q2 * q.q4)) * \
                         (2 * pow(q.q1, 2) + 2 * pow(q.q4, 2) - 1) / pz;      
// partial_nx / partial_q  
filtro->H[offset][6] = - 2 * (*z).dt * filtro->f * (2 * q.q1 * (q.q1 * (q.q1 * vx - q.q3 * vz + q.q4 * vy) \
                         + q.q2 * (q.q2 * vx + q.q3 * vy + q.q4 * vz) - q.q3 * (q.q1 * vz - q.q2 * vy + q.q3 * vx) \
                         + q.q4 * (q.q1 * vy + q.q2 * vz - q.q4 * vx)) + (2 * pow(q.q1, 2) + 2 * pow(q.q4, 2) - 1) \
                         * (q.q1 * vx - q.q3 * vz + q.q4 * vy)) / pz;    
filtro->H[offset][7] = - 2 * (*z).dt * filtro->f * (2 * pow(q.q1, 2) + 2 * pow(q.q4, 2) - 1) \
                         * (q.q2 * vx + q.q3 * vy + q.q4 * vz) / pz;    
filtro->H[offset][8] = 2 * (*z).dt * filtro->f * (2 * pow(q.q1, 2) + 2 * pow(q.q4, 2) - 1) \
                         * (q.q1 * vz - q.q2 * vy + q.q3 * vx) / pz;   
filtro->H[offset][9] = - 2 * (*z).dt * filtro->f * (2 * q.q4 * (q.q1 * (q.q1 * vx - q.q3 * vz + q.q4 * vy) \
                         + q.q2 * (q.q2 * vx + q.q3 * vy + q.q4 * vz) - q.q3 * (q.q1 * vz - q.q2 * vy + q.q3 * vx) \
                         + q.q4 * (q.q1 * vy + q.q2 * vz - q.q4 * vx)) + (2 * pow(q.q1, 2) + 2 * pow(q.q4, 2) - 1) \
                         * (q.q1 * vy + q.q2 * vz - q.q4 * vx)) / pz;  
// partial_ny / partial_pz
filtro->H[offset + 1][2] = ((*z).dt * filtro->f * (2 * pow(q.q1, 2) + 2 * pow(q.q4, 2) - 1)) * \
                        (q.q1 * (q.q1 * vy + q.q2 * vz - q.q4 * vx) + q.q2 * (q.q1 * vz - q.q2 * vy + q.q3 * vx) \
                        +q.q3 * (q.q3 * vx + q.q3 * vy + q.q4 * vz) - q.q4 * (q.q1 * vx - q.q3 * vz + q.q4 * vy) \
                        / pow(pz, 2));    
// partial_ny / partial_v   
filtro->H[offset + 1][3] = ((*z).dt * filtro->f * (2 * q.q1 * q.q4 - 2 * q.q2 * q.q3)) * \
                         (2 * pow(q.q1, 2) + 2 * pow(q.q4, 2) - 1) / pz;           
filtro->H[offset + 1][4] = - ((*z).dt * filtro->f * (2 * pow(q.q1, 2) + 2 * pow(q.q4, 2) - 1)) * \
                         (pow(q.q1, 2) - pow(q.q2, 2) + pow(q.q3, 2) - pow(q.q4, 2)) / pz;       
filtro->H[offset + 1][5] = - ((*z).dt * filtro->f * (2 * q.q1 * q.q2 + 2 * q.q3 * q.q4)) * \
                         (2 * pow(q.q1, 2) + 2 * pow(q.q4, 2) - 1) / pz;   
// partial_ny / partial_q  
filtro->H[offset + 1][6] = - 2 * (*z).dt * filtro->f * (2 * q.q1 * (q.q1 * (q.q1 * vy + q.q2 * vz - q.q4 * vx) \
                         + q.q2 * (q.q1 * vz - q.q2 * vy + q.q3 * vx) + q.q3 * (q.q2 * vx + q.q3 * vy + q.q4 * vz) \
                         - q.q4 * (q.q1 * vx - q.q3 * vz + q.q4 * vy)) + (2 * pow(q.q1, 2) + 2 * pow(q.q4, 2) - 1) \
                         * (q.q1 * vy + q.q2 * vz - q.q4 * vx)) / pz;     
filtro->H[offset + 1][7] = - 2 * (*z).dt * filtro->f * (2 * pow(q.q1, 2) + 2 * pow(q.q4, 2) - 1) \
                         * (q.q1 * vz - q.q2 * vy + q.q3 * vx) / pz; 
filtro->H[offset + 1][8] = - 2 * (*z).dt * filtro->f * (2 * pow(q.q1, 2) + 2 * pow(q.q4, 2) - 1) \
                         * (q.q2 * vx + q.q3 * vy + q.q4 * vz) / pz; 
filtro->H[offset + 1][9] = - 2 * (*z).dt * filtro->f * (2 * q.q4 * (q.q1 * (q.q1 * vy + q.q2 * vz - q.q4 * vx) \
                         + q.q2 * (q.q1 * vz - q.q2 * vy + q.q3 * vx) + q.q3 * (q.q2 * vx + q.q3 * vy + q.q4 * vz) \
                         - q.q4 * (q.q1 * vx - q.q3 * vz + q.q4 * vy)) - (2 * pow(q.q1, 2) + 2 * pow(q.q4, 2) - 1) \
                         * (q.q1 * vx - q.q3 * vz + q.q4 * vy)) / pz;   
}

void TOFS_states(ofs_ekf_t* filtro, int8_t offset){
// partial_tofs / partial_pz
filtro->H[offset + 0][2] = 1/(2 * (pow(filtro->states[N_P+N_V], 2) + pow(filtro->states[N_P+N_V+3], 2)) - 1); 
// partial_tofs / partial_q1  
filtro->H[offset + 0][N_P + N_V] = -4 * filtro->states[2] * filtro->states[N_P+N_V] * pow(filtro->H[offset + 0][2], 2);
// partial_tofs / partial_q4
filtro->H[offset + 0][N_P + N_V + 3] = -4 * filtro->states[2] * filtro->states[N_P+N_V+3] * pow(filtro->H[offset + 0][2], 2);
}