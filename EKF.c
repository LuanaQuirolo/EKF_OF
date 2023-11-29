#include "EKF.h"
#include <stdio.h> //TODO: BORRAR
void print_vector(char* string, int size, double* vec){
    printf("%s: \n", string);
    for (int i = 0; i < size; i++) {
        printf("%f ", vec[i]);
    }
    printf("\n ------------------------------------- \n");
}
void print_gain(char* string, int m, int n, double* vec) {
    printf("%s: \n", string);
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            printf("%f ", vec[i*n+j]);
        }
        printf("\n");
    }
    printf("\n ------------------------------------- \n");
}
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
    matmul_scalar(*filtro->cov, N_STATES, N_STATES, 0.1);
    mat_addeye(*filtro->cov, N_STATES);
    mat_zeros(*filtro->F, N_STATES, N_STATES);
    mat_addeye(*filtro->F, N_STATES);
    mat_zeros(*filtro->W, N_STATES, N_PROC_NOISE);
    /* Wk */
    // d(velocidad) / d(uax)
    filtro->W[N_P][0] = 1;
    // d(velocidad) / d(uay)
    filtro->W[N_P+1][1] = 1;
    // d(velocidad) / d(uaz)
    filtro->W[N_P+2][2] = 1;
    mat_zeros(*filtro->Q, N_PROC_NOISE, N_PROC_NOISE);
    mat_zeros(*filtro->R, N_CORR_NOISE, N_CORR_NOISE);
    filtro->Q[0][0] = 0.001; //uax
    filtro->Q[1][1] = 0.001; //uay
    filtro->Q[2][2] = 0.001; //uaz
    filtro->Q[3][3] = 0.006; //uwx
    filtro->Q[4][4] = 0.006; //uwy
    filtro->Q[5][5] = 0.006; //uwz
    filtro->Q[6][6] = 0.0001; //ubsz
    filtro->R[0][0] = U_A; //uax
    filtro->R[1][1] = U_A; //uay
    filtro->R[2][2] = U_A; //uaz
    filtro->Npix = 35; // Cantidad de píxeles
    filtro->FOV_OF = 4.2 * M_PI / 180; // FOV del sensor de OF
    filtro->f = filtro->Npix / (2 * atan2f(filtro->FOV_OF, 2));  // Factor de conversión
};
    
void prediction_step(ofs_ekf_t* filtro, mediciones_t u){
    
    /* Auxiliares */
    filtro->p[0] = filtro->states[0];
    filtro->p[1] = filtro->states[1];
    filtro->p[2] = filtro->states[2];
    filtro->v[0] = filtro->states[N_P+0];
    filtro->v[1] = filtro->states[N_P+1];
    filtro->v[2] = filtro->states[N_P+2];
    filtro->q.q1 = filtro->states[N_P + N_V + 0];
    filtro->q.q2 = filtro->states[N_P + N_V + 1];
    filtro->q.q3 = filtro->states[N_P + N_V + 2];
    filtro->q.q4 = filtro->states[N_P + N_V + 3];
    filtro->qa_meas.q1 = 0;
    filtro->qa_meas.q2 = u.ax;
    filtro->qa_meas.q3 = u.ay;
    filtro->qa_meas.q4 = u.az;
    filtro->qw_meas.q1 = 0;
    filtro->qw_meas.q2 = u.wx;
    filtro->qw_meas.q3 = u.wy;
    filtro->qw_meas.q4 = u.wz;

    /* Fk */
    //d(posicion) / d(velocidad)
    filtro->F[0][N_P] = u.dt;
    filtro->F[1][N_P + 1] = u.dt;
    filtro->F[2][N_P + 2] = u.dt;
    //d(velocidad) / d(quaternion)
    filtro->F[N_P][N_P+N_V] = u.dt * (2*(u.ax) * filtro->q.q1 - 2*(u.ay)*filtro->q.q4 + 2*(u.az)*filtro->q.q3);
    filtro->F[N_P+1][N_P+N_V] = u.dt * (2*(u.ax) * filtro->q.q4 + 2*(u.ay)*filtro->q.q1 - 2*(u.az)*filtro->q.q2);
    filtro->F[N_P+2][N_P+N_V] = u.dt * (-2*(u.ax) * filtro->q.q3 + 2*(u.ay)*filtro->q.q2 + 2*(u.az)*filtro->q.q1);
    filtro->F[N_P][N_P+N_V+1] = u.dt * (2*(u.ax) * filtro->q.q2 + 2*(u.ay)*filtro->q.q3 + 2*(u.az)*filtro->q.q4);
    filtro->F[N_P+1][N_P+N_V+1] = u.dt * (2*(u.ax) * filtro->q.q3 - 2*(u.ay)*filtro->q.q2 - 2*(u.az)*filtro->q.q1);
    filtro->F[N_P+2][N_P+N_V+1] = u.dt * (2*(u.ax) * filtro->q.q4 + 2*(u.ay)*filtro->q.q1 - 2*(u.az)*filtro->q.q2);
    filtro->F[N_P][N_P+N_V+2] = u.dt * (-2*(u.ax) * filtro->q.q3 + 2*(u.ay)*filtro->q.q2 + 2*(u.az)*filtro->q.q1);
    filtro->F[N_P+1][N_P+N_V+2] = u.dt * (2*(u.ax) * filtro->q.q2 + 2*(u.ay)*filtro->q.q3 + 2*(u.az)*filtro->q.q4);
    filtro->F[N_P+2][N_P+N_V+2] = u.dt * (-2*(u.ax) * filtro->q.q1 + 2*(u.ay)*filtro->q.q4 - 2*(u.az)*filtro->q.q3);
    filtro->F[N_P][N_P+N_V+3] = u.dt * (-2*(u.ax) * filtro->q.q4 - 2*(u.ay)*filtro->q.q1 + 2*(u.az)*filtro->q.q2);
    filtro->F[N_P+1][N_P+N_V+3] = u.dt * (2*(u.ax) * filtro->q.q1 - 2*(u.ay)*filtro->q.q4 + 2*(u.az)*filtro->q.q3);
    filtro->F[N_P+2][N_P+N_V+3] = u.dt * (2*(u.ax) * filtro->q.q2 + 2*(u.ay)*filtro->q.q3 + 2*(u.az)*filtro->q.q4);
    //d(q) / d(q1)
    //filtro->F[N_P+N_V][N_P+N_V] = 1;
    filtro->F[N_P+N_V+1][N_P+N_V] = u.wx * u.dt / 2;
    filtro->F[N_P+N_V+2][N_P+N_V] = u.wy * u.dt / 2;
    filtro->F[N_P+N_V+3][N_P+N_V] = u.wz * u.dt / 2;
    //d(q) / d(q2)
    filtro->F[N_P+N_V][N_P+N_V+1] = -u.wx * u.dt / 2;
    //filtro->F[N_P+N_V+1][N_P+N_V+1] = 1;
    filtro->F[N_P+N_V+2][N_P+N_V+1] = -u.wz * u.dt / 2;
    filtro->F[N_P+N_V+3][N_P+N_V+1] = u.wy * u.dt / 2;
    //d(q) / d(q3)
    filtro->F[N_P+N_V][N_P+N_V+2] = -u.wy * u.dt / 2;
    filtro->F[N_P+N_V+1][N_P+N_V+2] = u.wz * u.dt / 2;
    //filtro->F[N_P+N_V+2][N_P+N_V+2] = 1;
    filtro->F[N_P+N_V+3][N_P+N_V+2] = -u.wx * u.dt / 2;
    //d(q) / d(q4)
    filtro->F[N_P+N_V][N_P+N_V+3] = -u.wz * u.dt / 2;
    filtro->F[N_P+N_V+1][N_P+N_V+3] = -u.wy * u.dt / 2;
    filtro->F[N_P+N_V+2][N_P+N_V+3] = u.wx * u.dt / 2;
    //filtro->F[N_P+N_V+3][N_P+N_V+3] = 1;
    /* Covarianza a priori */
    transpose(*filtro->W, *filtro->Wt, N_STATES, N_PROC_NOISE);
    mulmat(*filtro->Q, *filtro->Wt, *filtro->aux4, N_PROC_NOISE, N_PROC_NOISE, N_STATES); // aux4 = Q * Wt
    mulmat(*filtro->W, *filtro->aux4, *filtro->aux6, N_STATES, N_PROC_NOISE, N_STATES); // aux6 = W * Q * Wt
    transpose(*filtro->F, *filtro->Ft, N_STATES, N_STATES);
    mulmat(*filtro->cov, *filtro->Ft, *filtro->aux5, N_STATES, N_STATES, N_STATES); // aux5 = cov * Ft
    mulmat(*filtro->F, *filtro->aux5, *filtro->aux7, N_STATES, N_STATES, N_STATES); // aux7 = F * cov * Ft
    add(*filtro->aux6, *filtro->aux7, *filtro->cov, N_STATES, N_STATES); // cov = F * cov * Ft +  W * Q * Wt
    /* Paso de prediccion */
    // Posicion
    matmul_scalar2(filtro->v, filtro->aux2, N_V, 1, u.dt); // aux2 = Vk = dt * Vk
    add(filtro->p, filtro->aux2, filtro->states, N_P, 1); // Pk+1 = Pk + dt * Vk
    // Velocidad
    // Queremos pasar la acel de cuerpo a mundo
    filtro->aux = quat_mult(filtro->qa_meas, quat_conjugate(filtro->q)); // aux = (filtro->qa_meas) * -q
    filtro->aux = quat_mult(filtro->q, filtro->aux); // aux = q * (filtro->qa_meas) * -q
    quat_sub(&filtro->aux, filtro->aux, filtro->qg); // aux = q * (filtro->qa_meas) * -q - g
    quat2vec(filtro->aux, filtro->aux3); //aux3 = aux
    matmul_scalar(filtro->aux3, N_V, 1, u.dt); // aux3 = dt (a_meas - g)
    add(filtro->v, filtro->aux3, filtro->states + N_P, N_V, 1); // Vk+1 = Vk + dt * (a - g)
    // Quaternion
    filtro->aux = quat_mult(filtro->qw_meas, filtro->q); // aux = qw * qk
    quat_scalar(&filtro->aux, u.dt / 2); // aux = (dt / 2) qw * qk
    quat_add(&filtro->aux, filtro->aux, filtro->q); // aux = qk + (dt / 2) qw * qk
    quat_Normalization(&filtro->aux);
    filtro->states[N_P + N_V] = filtro->aux.q1;
    filtro->states[N_P + N_V + 1] = filtro->aux.q2;
    filtro->states[N_P + N_V + 2] = filtro->aux.q3;
    filtro->states[N_P + N_V + 3] = filtro->aux.q4;
};

void correction_step(ofs_ekf_t* filtro, mediciones_t* z){

mat_zeros(*filtro->H, N_OBS_11, N_STATES);
filtro->p[0] = filtro->states[0];
filtro->p[1] = filtro->states[1];
filtro->p[2] = filtro->states[2];
filtro->v[0] = filtro->states[N_P+0];
filtro->v[1] = filtro->states[N_P+1];
filtro->v[2] = filtro->states[N_P+2];
filtro->q.q1 = filtro->states[N_P + N_V + 0];
filtro->q.q2 = filtro->states[N_P + N_V + 1];
filtro->q.q3 = filtro->states[N_P + N_V + 2];
filtro->q.q4 = filtro->states[N_P + N_V + 3];

/*************************** Jacobiano **************************/
if(filtro->beta == 1 && filtro->gamma == 1){
    IMU_states(filtro, 0);
    OFS_states(filtro, N_IMU, z);
    TOFS_states(filtro, N_IMU + N_OFS);
    filtro->meas_counter = N_OBS_11;
}
else if (filtro->beta == 0 && filtro->gamma == 1){
    IMU_states(filtro, 0);
    TOFS_states(filtro, N_IMU);
    filtro->meas_counter = N_OBS_10;
}
else if (filtro->beta == 1 && filtro->gamma == 0){
    IMU_states(filtro, 0);
    OFS_states(filtro, N_IMU, z);
    filtro->meas_counter = N_OBS_01;
}
else {
    IMU_states(filtro, 0);
    filtro->meas_counter = N_OBS_00;
}
//print_gain("H", filtro->meas_counter, N_STATES, *filtro->H);
/*************************** Mediciones predichas **************************/
mat_zeros(filtro->exp_meas, N_OBS_11, 1);
// Aceleracion
filtro->aux = quat_mult(filtro->qg, filtro->q); // aux = qg * qk (de mundo a cuerpo)
filtro->aux = quat_mult(quat_conjugate(filtro->q), filtro->aux); // aux = q- * qg * q (de mundo a cuerpo)
quat2vec(filtro->aux, filtro->exp_meas); // measurements = q- * qg * q
filtro->meas[0] = z->ax;
filtro->meas[1] = z->ay;
filtro->meas[2] = z->az;
// Flujo optico (Si corresponde)
if(filtro->beta == 1){
    filtro->aux = vec2quat(filtro->states + 3); //Paso la velocidad respecto al mundo a quat para convertirlo a body
    filtro->aux = quat_mult(filtro->aux, filtro->q); // aux = v * q
    filtro->aux = quat_mult(quat_conjugate(filtro->q), filtro->aux); // aux = q- * v * q

    filtro->exp_meas[N_OBS_00 + 0] = -(*z).dt * filtro->f * ((2 * (pow(filtro->q.q1, 2) + pow(filtro->q.q4, 2)) - 1) / \
                                         filtro->states[2] * filtro->aux.q2 + (*z).wy);
    filtro->exp_meas[N_OBS_00 + 1] = -(*z).dt * filtro->f * ((2 * (pow(filtro->q.q1, 2) + pow(filtro->q.q4, 2)) - 1) / \
                                         filtro->states[2] * filtro->aux.q3 - (*z).wx);  
    filtro->R[N_OBS_00][N_OBS_00] = U_FLOW; //u flowx    
    filtro->R[N_OBS_00 + 1][N_OBS_00 + 1] = U_FLOW; //u flowy  
    filtro->meas[N_OBS_00 + 0] = z->ofx;
    filtro->meas[N_OBS_00 + 1] = z->ofy;                                                                
}
// Sensor de distancia (Si corresponde)
if(filtro->beta == 1 && filtro->gamma == 1){
    filtro->exp_meas[N_OBS_11 - 1] = filtro->states[2] / (2 * (pow(filtro->q.q1, 2) + pow(filtro->q.q4, 2)) - 1);
    filtro->R[N_OBS_11 - 1][N_OBS_11 - 1] = U_RANGE; // distance
    filtro->meas[N_OBS_11 - 1] = z->range; 
}
else if(filtro->beta == 0 && filtro->gamma == 1){
    filtro->exp_meas[N_OBS_10 - 1] = filtro->states[2] / (2 * (pow(filtro->q.q1, 2) + pow(filtro->q.q4, 2)) - 1);
    filtro->R[N_OBS_10 - 1][N_OBS_10 - 1] = U_RANGE; // distance  
    filtro->meas[N_OBS_10 - 1] = z->range;
}
/*************************** Ganancia de Kalman **************************/

transpose(*filtro->H, *filtro->Ht, filtro->meas_counter, N_STATES); 
mulmat(*filtro->cov, *filtro->Ht, *filtro->aux9, N_STATES, N_STATES, filtro->meas_counter); // aux9 = cov * Ht
mulmat(*filtro->H, *filtro->aux9, *filtro->aux10, filtro->meas_counter, N_STATES, filtro->meas_counter); // aux10 = H * cov * Ht
accum(*filtro->aux10, *filtro->R, filtro->meas_counter, filtro->meas_counter); // aux10 = H * cov * Ht + R
cholsl(*filtro->aux10, *filtro->aux11, filtro->aux12, filtro->meas_counter); // aux11 = inv(H * cov * Ht + R)
mulmat(*filtro->aux9, *filtro->aux11, *filtro->G, N_STATES, filtro->meas_counter, filtro->meas_counter); // G = cov * Ht * inv(H * cov * Ht + R)

/*************************** Correccion **************************/
sub(filtro->meas, filtro->exp_meas, filtro->aux12, filtro->meas_counter); // aux12 = z_medido - z_esperado
mulvec(*filtro->G, filtro->aux12, filtro->aux8, N_STATES, filtro->meas_counter); // aux8 = G(z_medido - z_esperado)
//print_vector("Mediciones", filtro->meas_counter, filtro->meas);
print_vector("Mediciones esperadas", filtro->meas_counter, filtro->exp_meas);
print_vector("Innovacion", filtro->meas_counter, filtro->aux8);
accum(filtro->states, filtro->aux8, N_STATES, 1); // estado = estado + G(z_medido - z_esperado)
filtro->q.q1 = filtro->states[N_P + N_V + 0];
filtro->q.q2 = filtro->states[N_P + N_V + 1];
filtro->q.q3 = filtro->states[N_P + N_V + 2];
filtro->q.q4 = filtro->states[N_P + N_V + 3];
quat_Normalization(&filtro->q);
filtro->states[N_P + N_V] = filtro->q.q1;
filtro->states[N_P + N_V + 1] = filtro->q.q2;
filtro->states[N_P + N_V + 2] = filtro->q.q3;
filtro->states[N_P + N_V + 3] = filtro->q.q4;
/*************************** Covarianza **************************/
mulmat(*filtro->G, *filtro->H, *filtro->aux9, N_STATES, filtro->meas_counter, N_STATES); //  aux9 = G * H
mulmat(*filtro->aux9, *filtro->cov, *filtro->aux7, N_STATES, N_STATES, N_STATES); // aux7 = G * H * cov_priori
mat_negate(*filtro->aux7, N_STATES, N_STATES); // aux7 = - G * H * cov_priori
accum(*filtro->cov, *filtro->aux7, N_STATES, N_STATES); // cov = cov_priori - G * H * cov_priori
(*filtro).beta = 0;
(*filtro).gamma = 0;
}

void IMU_states(ofs_ekf_t* filtro, int8_t offset){

// partial_ax / partial_q
filtro->H[offset + 0][N_P + N_V + 0] = -2 * g * (filtro->q.q3);     
filtro->H[offset + 0][N_P + N_V + 1] = 2 * g * (filtro->q.q4);     
filtro->H[offset + 0][N_P + N_V + 2] = -2 * g * (filtro->q.q1);     
filtro->H[offset + 0][N_P + N_V + 3] = 2 * g * (filtro->q.q2);  
// partial_ay / partial_q   
filtro->H[offset + 1][N_P + N_V + 0] = 2 * g * (filtro->q.q2);      
filtro->H[offset + 1][N_P + N_V + 1] = 2 * g * (filtro->q.q1);      
filtro->H[offset + 1][N_P + N_V + 2] = 2 * g * (filtro->q.q4);     
filtro->H[offset + 1][N_P + N_V + 3] = 2 * g * (filtro->q.q3);
// partial_az / partial_q  
filtro->H[offset + 2][N_P + N_V + 0] = 2 * g * (filtro->q.q1);      
filtro->H[offset + 2][N_P + N_V + 1] = -2 * g * (filtro->q.q2);      
filtro->H[offset + 2][N_P + N_V + 2] = -2 * g * (filtro->q.q3);     
filtro->H[offset + 2][N_P + N_V + 3] = 2 * g * (filtro->q.q4);
}

void OFS_states(ofs_ekf_t* filtro, int8_t offset, mediciones_t *z){

// partial_nx / partial_pz
filtro->H[offset][2] = ((*z).dt * filtro->f * (2 * pow(filtro->q.q1, 2) + 2 * pow(filtro->q.q4, 2) - 1)) * \
                        (filtro->q.q1 * (filtro->q.q1 * filtro->v[0] - filtro->q.q3 * filtro->v[2] + \
                        filtro->q.q4 * filtro->v[1]) - filtro->q.q2 * (-filtro->q.q2 * filtro->v[0] - \
                        filtro->q.q3 * filtro->v[1] - filtro->q.q4 * filtro->v[2]) \
                        -filtro->q.q3 * (filtro->q.q1 * filtro->v[2] - filtro->q.q2 * filtro->v[1] + \
                        filtro->q.q3 * filtro->v[0]) + filtro->q.q4 * (filtro->q.q1 * filtro->v[1] + \
                        filtro->q.q2 * filtro->v[2] - filtro->q.q4 * filtro->v[0]) \
                        / pow(filtro->p[2], 2));     
// partial_nx / partial_v   
filtro->H[offset][N_P + 0] = - ((*z).dt * filtro->f * (2 * pow(filtro->q.q1, 2) + 2 * pow(filtro->q.q4, 2) - 1)) * \
                         (pow(filtro->q.q1, 2) + pow(filtro->q.q2, 2) - pow(filtro->q.q3, 2) - pow(filtro->q.q4, 2)) / filtro->p[2];       
filtro->H[offset][N_P + 1] = - ((*z).dt * filtro->f * (2 *filtro->q.q1 *filtro->q.q4 + 2 *filtro->q.q2 *filtro->q.q3)) * \
                         (2 * pow(filtro->q.q1, 2) + 2 * pow(filtro->q.q4, 2) - 1) / filtro->p[2];        
filtro->H[offset][N_P + 2] = - ((*z).dt * filtro->f * (-2 *filtro->q.q1 *filtro->q.q3 + 2 *filtro->q.q2 *filtro->q.q4)) * \
                         (2 * pow(filtro->q.q1, 2) + 2 * pow(filtro->q.q4, 2) - 1) / filtro->p[2];      
// partial_nx / partial_q  
filtro->H[offset][N_P + N_V + 0] = - 2 * (*z).dt * filtro->f * (2 *filtro->q.q1 * (filtro->q.q1 * (filtro->q.q1 * filtro->v[0] -\
                        filtro->q.q3 * filtro->v[2] +filtro->q.q4 * filtro->v[1]) \
                        +filtro->q.q2 * (filtro->q.q2 * filtro->v[0] +filtro->q.q3 * filtro->v[1] +filtro->q.q4 * filtro->v[2])\
                        -filtro->q.q3 * (filtro->q.q1 * filtro->v[2] -filtro->q.q2 * filtro->v[1] +filtro->q.q3 * filtro->v[0]) \
                        +filtro->q.q4 * (filtro->q.q1 * filtro->v[1] +filtro->q.q2 * filtro->v[2] -filtro->q.q4 * filtro->v[0])) \
                        + (2 * pow(filtro->q.q1, 2) + 2 * pow(filtro->q.q4, 2) - 1) \
                        * (filtro->q.q1 * filtro->v[0] -filtro->q.q3 * filtro->v[2] +filtro->q.q4 * filtro->v[1])) / filtro->p[2];    
filtro->H[offset][N_P + N_V + 1] = - 2 * (*z).dt * filtro->f * (2 * pow(filtro->q.q1, 2) + 2 * pow(filtro->q.q4, 2) - 1) \
                         * (filtro->q.q2 * filtro->v[0] +filtro->q.q3 * filtro->v[1] +filtro->q.q4 * filtro->v[2]) / filtro->p[2];    
filtro->H[offset][N_P + N_V + 2] = 2 * (*z).dt * filtro->f * (2 * pow(filtro->q.q1, 2) + 2 * pow(filtro->q.q4, 2) - 1) \
                         * (filtro->q.q1 * filtro->v[2] -filtro->q.q2 * filtro->v[1] +filtro->q.q3 * filtro->v[0]) / filtro->p[2];   
filtro->H[offset][N_P + N_V + 3] = - 2 * (*z).dt * filtro->f * (2 *filtro->q.q4 * (filtro->q.q1 * (filtro->q.q1 * filtro->v[0] \
                        - filtro->q.q3 * filtro->v[2] +filtro->q.q4 * filtro->v[1]) +filtro->q.q2 * (filtro->q.q2 * filtro->v[0] +\
                        filtro->q.q3 * filtro->v[1] +filtro->q.q4 * filtro->v[2]) -filtro->q.q3 * (filtro->q.q1 * filtro->v[2] \
                        - filtro->q.q2 * filtro->v[1] +filtro->q.q3 * filtro->v[0]) +filtro->q.q4 * (filtro->q.q1 * filtro->v[1] \
                        + filtro->q.q2 * filtro->v[2] -filtro->q.q4 * filtro->v[0])) + (2 * pow(filtro->q.q1, 2) + 2 * pow(filtro->q.q4, 2) - 1) \
                        * (filtro->q.q1 * filtro->v[1] +filtro->q.q2 * filtro->v[2] -filtro->q.q4 * filtro->v[0])) / filtro->p[2];  
// partial_ny / partial_pz
filtro->H[offset + 1][2] = ((*z).dt * filtro->f * (2 * pow(filtro->q.q1, 2) + 2 * pow(filtro->q.q4, 2) - 1)) * \
                        (filtro->q.q1 * (filtro->q.q1 * filtro->v[1] +filtro->q.q2 * filtro->v[2] -filtro->q.q4 * filtro->v[0]) \
                        +filtro->q.q2 * (filtro->q.q1 * filtro->v[2] -filtro->q.q2 * filtro->v[1] +filtro->q.q3 * filtro->v[0]) \
                        + filtro->q.q3 * (filtro->q.q2 * filtro->v[0] +filtro->q.q3 * filtro->v[1] +filtro->q.q4 * filtro->v[2]) \
                        -filtro->q.q4 * (filtro->q.q1 * filtro->v[0] -filtro->q.q3 * filtro->v[2] +filtro->q.q4 * filtro->v[1]) \
                        / pow(filtro->p[2], 2));    
// partial_ny / partial_v   
filtro->H[offset + 1][N_P + 0] = ((*z).dt * filtro->f * (2 *filtro->q.q1 *filtro->q.q4 - 2 *filtro->q.q2 *filtro->q.q3)) * \
                         (2 * pow(filtro->q.q1, 2) + 2 * pow(filtro->q.q4, 2) - 1) / filtro->p[2];           
filtro->H[offset + 1][N_P + 1] = - ((*z).dt * filtro->f * (2 * pow(filtro->q.q1, 2) + 2 * pow(filtro->q.q4, 2) - 1)) * \
                         (pow(filtro->q.q1, 2) - pow(filtro->q.q2, 2) + pow(filtro->q.q3, 2) - pow(filtro->q.q4, 2)) / filtro->p[2];       
filtro->H[offset + 1][N_P + 2] = - ((*z).dt * filtro->f * (2 *filtro->q.q1 *filtro->q.q2 + 2 *filtro->q.q3 *filtro->q.q4)) * \
                         (2 * pow(filtro->q.q1, 2) + 2 * pow(filtro->q.q4, 2) - 1) / filtro->p[2];   
// partial_ny / partial_q  
filtro->H[offset + 1][N_P + N_V + 0] = - 2 * (*z).dt * filtro->f * (2 *filtro->q.q1 * (filtro->q.q1 * (filtro->q.q1 * filtro->v[1] \
                           +filtro->q.q2 * filtro->v[2] -filtro->q.q4 * filtro->v[0]) \
                         +filtro->q.q2 * (filtro->q.q1 * filtro->v[2] -filtro->q.q2 * filtro->v[1] +filtro->q.q3 * filtro->v[0]) \
                         +filtro->q.q3 * (filtro->q.q2 * filtro->v[0] +filtro->q.q3 * filtro->v[1] +filtro->q.q4 * filtro->v[2]) \
                         -filtro->q.q4 * (filtro->q.q1 * filtro->v[0] -filtro->q.q3 * filtro->v[2] +filtro->q.q4 * filtro->v[1])) \
                         + (2 * pow(filtro->q.q1, 2) + 2 * pow(filtro->q.q4, 2) - 1) \
                         * (filtro->q.q1 * filtro->v[1] +filtro->q.q2 * filtro->v[2] -filtro->q.q4 * filtro->v[0])) / filtro->p[2];     
filtro->H[offset + 1][N_P + N_V + 1] = - 2 * (*z).dt * filtro->f * (2 * pow(filtro->q.q1, 2) + 2 * pow(filtro->q.q4, 2) - 1) \
                         * (filtro->q.q1 * filtro->v[2] -filtro->q.q2 * filtro->v[1] +filtro->q.q3 * filtro->v[0]) / filtro->p[2]; 
filtro->H[offset + 1][N_P + N_V + 2] = - 2 * (*z).dt * filtro->f * (2 * pow(filtro->q.q1, 2) + 2 * pow(filtro->q.q4, 2) - 1) \
                         * (filtro->q.q2 * filtro->v[0] +filtro->q.q3 * filtro->v[1] +filtro->q.q4 * filtro->v[2]) / filtro->p[2]; 
filtro->H[offset + 1][N_P + N_V + 3] = - 2 * (*z).dt * filtro->f * (2 *filtro->q.q4 * (filtro->q.q1 * (filtro->q.q1 * filtro->v[1] \
                            +filtro->q.q2 * filtro->v[2] -filtro->q.q4 * filtro->v[0]) \
                         +filtro->q.q2 * (filtro->q.q1 * filtro->v[2] -filtro->q.q2 * filtro->v[1] +filtro->q.q3 * filtro->v[0]) \
                         +filtro->q.q3 * (filtro->q.q2 * filtro->v[0] +filtro->q.q3 * filtro->v[1] +filtro->q.q4 * filtro->v[2]) \
                         -filtro->q.q4 * (filtro->q.q1 * filtro->v[0] -filtro->q.q3 * filtro->v[2] +filtro->q.q4 * filtro->v[1])) \
                         - (2 * pow(filtro->q.q1, 2) + 2 * pow(filtro->q.q4, 2) - 1) \
                         * (filtro->q.q1 * filtro->v[0] -filtro->q.q3 * filtro->v[2] +filtro->q.q4 * filtro->v[1])) / filtro->p[2]; 
}

void TOFS_states(ofs_ekf_t* filtro, int8_t offset){
// partial_tofs / partial_pz
filtro->H[offset + 0][2] = 1/(2 * (pow(filtro->q.q1, 2) + pow(filtro->q.q4, 2)) - 1); 
// partial_tofs / partial_q1  
filtro->H[offset + 0][N_P + N_V] = -4 * filtro->p[2] * filtro->q.q1 * pow(filtro->H[offset + 0][2], 2);
// partial_tofs / partial_q4
filtro->H[offset + 0][N_P + N_V + 3] = -4 * filtro->p[2] * filtro->q.q4 * pow(filtro->H[offset + 0][2], 2);
}