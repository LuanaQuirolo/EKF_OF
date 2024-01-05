#include "EKF.h"
#include <stdio.h> //TODO: BORRAR

void print_vector(char* string, int size, double* vec){
    printf("%s: \n", string);
    for (int i = 0; i < size; i++) {
        printf("%f ", vec[i]);
    }
    printf("\n ------------------------------------- \n");
}
void print_gain(char* string, int m, int n, double vec[m][n]) {
    printf("%s: \n", string);
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            printf("%.9f ", vec[i][j]);
        }
        printf("\n");
    }
    printf("\n ------------------------------------- \n");
}
void ofs_ekf_init(ofs_ekf_t* filtro, double m[N_MAG]){
    filtro->N =  N_STATES; // Cantidad de estados
    filtro->m[0] = m[0];
    filtro->m[1] = m[1];
    filtro->m[2] = m[2];
    filtro->qm = vec2quat(filtro->m);
    filtro->beta = 0; // Indica si hay una lectura nueva del OFS
    filtro->gamma = 0; // Indica si hay una lectura nueva del sensor de distancia
    double temp[3] = {0, 0, -g};
    filtro->qg = vec2quat(temp);
    vec_zeros(N_STATES, filtro->states); //p, v, q
    filtro->states[N_P + N_V] = 1; //q1 = 1;
    mat_zeros(N_STATES, N_STATES, filtro->cov); // Matriz de covarianza de estados
    mat_addeye(N_STATES, filtro->cov);
    matmul_scalar(N_STATES, N_STATES, filtro->cov, 10);
    mat_zeros(N_STATES, N_STATES, filtro->F);
    mat_addeye(N_STATES, filtro->F);
    mat_zeros(N_STATES, N_PROC_NOISE, filtro->W);
    /* Wk */
    // d(velocidad) / d(uax)
    filtro->W[N_P][0] = 1;
    // d(velocidad) / d(uay)
    filtro->W[N_P+1][1] = 1;
    // d(velocidad) / d(uaz)
    filtro->W[N_P+2][2] = 1;
    // d(q1) / d(uw)
    filtro->W[N_P+N_V][3] = 1;
    // d(q2) / d(uw)
    filtro->W[N_P+N_V+1][4] = 1;
    // d(q3) / d(uw)
    filtro->W[N_P+N_V+2][5] = 1;
    // d(q4) / d(uw)
    filtro->W[N_P+N_V+3][6] = 1;
    transpose(N_STATES, N_PROC_NOISE, filtro->W, filtro->Wt);
    mat_zeros(N_PROC_NOISE, N_PROC_NOISE, filtro->Q);
    filtro->Q[0][0] = U_A; //uax
    filtro->Q[1][1] = U_A; //uay
    filtro->Q[2][2] = U_A; //uaz
    filtro->Q[3][3] = U_W; //uwx
    filtro->Q[4][4] = U_W; //uwy
    filtro->Q[5][5] = U_W; //uwz
    filtro->Q[6][6] = U_W; //uw
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
    //Como W no es variable, se calcula la transpuesta al inicio por unica vez
    mulmat(N_PROC_NOISE, N_PROC_NOISE, N_STATES, filtro->Q, filtro->Wt, filtro->aux4); // aux4 = Q * Wt
    mulmat(N_STATES, N_PROC_NOISE, N_STATES, filtro->W, filtro->aux4, filtro->aux6); // aux6 = W * Q * Wt
    transpose(N_STATES, N_STATES, filtro->F, filtro->Ft);
    mulmat(N_STATES, N_STATES, N_STATES, filtro->cov, filtro->Ft, filtro->aux5); // aux5 = cov * Ft
    mulmat(N_STATES, N_STATES, N_STATES, filtro->F, filtro->aux5, filtro->aux7); // aux7 = F * cov * Ft
    add(N_STATES, N_STATES, filtro->aux6, filtro->aux7, filtro->cov); // cov = F * cov * Ft +  W * Q * Wt
    /* Paso de prediccion */
    // Posicion
    vecmul_scalar2(N_V, filtro->v, filtro->aux2, u.dt); // aux2 = Vk = dt * Vk
    add_vec(N_P, filtro->p, filtro->aux2, filtro->states); // Pk+1 = Pk + dt * Vk
    if(filtro->states[2] < MIN_HEIGHT){
    	filtro->states[2] = MIN_HEIGHT;
    }
    else if (filtro->states[2] > MAX_HEIGHT){
    	filtro->states[2] = MAX_HEIGHT;
    }
    // Velocidad
    // Queremos pasar la acel de cuerpo a mundo
    filtro->aux13 = quat_mult(filtro->qa_meas, quat_conjugate(filtro->q)); // aux13 = (filtro->qa_meas) * -q
    filtro->aux = quat_mult(filtro->q, filtro->aux13); // aux13 = q * (filtro->qa_meas) * -q
    quat_sub(&filtro->aux13, filtro->aux, filtro->qg); // aux13 = q * (filtro->qa_meas) * -q - g
    quat2vec(filtro->aux13, filtro->aux3); //aux3 = aux13
    vecmul_scalar(N_V, filtro->aux3, u.dt); // aux3 = dt (a_meas - g)
    add_vec(N_V, filtro->v, filtro->aux3, &filtro->states[N_P]); // Vk+1 = Vk + dt * (a - g)
    // Quaternion
    filtro->aux = quat_mult(filtro->q, filtro->qw_meas); // aux = qk * qw
    quat_scalar(&filtro->aux, u.dt / 2); // aux = (dt / 2) qk * qw
    quat_add(&filtro->aux13, filtro->aux, filtro->q); // aux = qk + (dt / 2) qk * qw
    quat_Normalization(&filtro->aux13);
    filtro->states[N_P + N_V] = filtro->aux13.q1;
    filtro->states[N_P + N_V + 1] = filtro->aux13.q2;
    filtro->states[N_P + N_V + 2] = filtro->aux13.q3;
    filtro->states[N_P + N_V + 3] = filtro->aux13.q4;
};

void correction_step(ofs_ekf_t* filtro, mediciones_t* z){

mat_zeros(N_OBS_11, N_STATES, filtro->H);
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
if (filtro->beta == 1 && filtro->gamma == 1){
    IMU_states(filtro);
    MAG_states(filtro);
    OFS_states(filtro, z);
    TOFS_states(filtro, N_IMU + N_MAG + N_OFS);
    filtro->meas_counter = N_OBS_11;
}
else if (filtro->beta == 0 && filtro->gamma == 1){
    IMU_states(filtro);
    MAG_states(filtro);
    TOFS_states(filtro, N_IMU);
    filtro->meas_counter = N_OBS_10;
}
else if (filtro->beta == 1 && filtro->gamma == 0){
    IMU_states(filtro);
    MAG_states(filtro);
    OFS_states(filtro, z);
    filtro->meas_counter = N_OBS_01;
}
else {
    IMU_states(filtro);
    MAG_states(filtro);
    filtro->meas_counter = N_OBS_00;
}
/*************************** Mediciones predichas **************************/
vec_zeros(N_OBS_11, filtro->exp_meas);
// Aceleracion
filtro->aux13 = quat_mult(filtro->qg, filtro->q); // aux13 = qg * qk (de mundo a cuerpo)
filtro->aux = quat_mult(quat_conjugate(filtro->q), filtro->aux13); // aux = q- * qg * q (de mundo a cuerpo)
quat2vec(filtro->aux, filtro->exp_meas); // measurements = q- * qg * q
filtro->meas[0] = z->ax;
filtro->meas[1] = z->ay;
filtro->meas[2] = z->az;
//Campo magnetico
filtro->aux13 = quat_mult(filtro->qm, filtro->q); // aux = qm * qk (de mundo a cuerpo)
filtro->aux = quat_mult(quat_conjugate(filtro->q), filtro->aux13); // aux = q- * qm * q (de mundo a cuerpo)
quat2vec(filtro->aux, &filtro->exp_meas[N_IMU]); // measurements = q- * qm * q
filtro->meas[N_IMU + 0] = z->mx;
filtro->meas[N_IMU + 1] = z->my;
filtro->meas[N_IMU + 2] = z->mz;
// Redefino R
mat_zeros(N_CORR_NOISE, N_CORR_NOISE, filtro->R);
filtro->R[0][0] = U_A; //uax
filtro->R[1][1] = U_A; //uay
filtro->R[2][2] = U_A; //uaz
filtro->R[3][3] = U_MAG; //umagx
filtro->R[4][4] = U_MAG; //umagy
filtro->R[5][5] = U_MAG; //umagz
// Flujo optico (Si corresponde)
if(filtro->beta == 1){
    filtro->aux13 = vec2quat(&filtro->states[N_P]); //Paso la velocidad respecto al mundo a quat para convertirlo a body
    filtro->aux = quat_mult(filtro->aux13, filtro->q); // aux = v * q
    filtro->aux13 = quat_mult(quat_conjugate(filtro->q), filtro->aux); // aux13 = q- * v * q
    filtro->exp_meas[N_OBS_00 + 0] = -(*z).tau * filtro->f * (filtro->aux13.q2 * (2 * (pow(filtro->q.q1, 2) + pow(filtro->q.q4, 2)) - 1) / \
                                         filtro->p[2]  + (*z).wy);
    filtro->exp_meas[N_OBS_00 + 1] = -(*z).tau * filtro->f * (filtro->aux13.q3 * (2 * (pow(filtro->q.q1, 2) + pow(filtro->q.q4, 2)) - 1) / \
                                         filtro->p[2]  - (*z).wx);  
    filtro->R[N_OBS_00][N_OBS_00] = U_FLOW; //u flowx    
    filtro->R[N_OBS_00 + 1][N_OBS_00 + 1] = U_FLOW; //u flowy  
    filtro->meas[N_OBS_00 + 0] = z->ofx;
    filtro->meas[N_OBS_00 + 1] = z->ofy;                                                                
}

// Sensor de distancia (Si corresponde)
if(filtro->beta == 1 && filtro->gamma == 1){
    filtro->exp_meas[N_OBS_11 - 1] = filtro->p[2] / (2 * (pow(filtro->q.q1, 2) + pow(filtro->q.q4, 2)) - 1);
    filtro->R[N_OBS_11 - 1][N_OBS_11 - 1] = U_RANGE; // distance
    filtro->meas[N_OBS_11 - 1] = z->range; 
}
else if(filtro->beta == 0 && filtro->gamma == 1){
    filtro->exp_meas[N_OBS_10 - 1] = filtro->p[2] / (2 * (pow(filtro->q.q1, 2) + pow(filtro->q.q4, 2)) - 1);
    filtro->R[N_OBS_00][N_OBS_00] = U_RANGE; // distance  
    filtro->meas[N_OBS_00] = z->range;
}
/*************************** Ganancia de Kalman **************************/
transpose(filtro->meas_counter, N_STATES, filtro->H, filtro->Ht); 
mulmat(N_STATES, N_STATES, filtro->meas_counter, filtro->cov, filtro->Ht, filtro->aux9); // aux9 = cov * Ht
mulmat(filtro->meas_counter, N_STATES, filtro->meas_counter, filtro->H, filtro->aux9, filtro->aux10); // aux10 = H * cov * Ht
accum(filtro->meas_counter, filtro->meas_counter, filtro->aux10, filtro->R); // aux10 = H * cov * Ht + R
//print_gain("(H * cov * Ht + R): ",filtro->meas_counter, filtro->meas_counter,  filtro->aux10);
cholsl(filtro->meas_counter, filtro->aux10, filtro->aux11, filtro->aux12); // aux11 = inv(H * cov * Ht + R)
//print_gain("inv(H * cov * Ht + R): ",filtro->meas_counter, filtro->meas_counter,  filtro->aux11);
mulmat(N_STATES, filtro->meas_counter, filtro->meas_counter, filtro->aux9, filtro->aux11, filtro->G); // G = cov * Ht * inv(H * cov * Ht + R) //TODO: EN VEZ DE 10 era 11
//print_gain("G: ", N_STATES, filtro->meas_counter,  filtro->G);
/*************************** Correccion **************************/
sub_vec(filtro->meas, filtro->exp_meas, filtro->aux12, filtro->meas_counter); // aux12 = z_medido - z_esperado
mulvec(N_STATES, filtro->meas_counter, filtro->G, filtro->aux12, filtro->aux8); // aux8 = G(z_medido - z_esperado)
//print_vector("Mediciones", filtro->meas_counter, filtro->meas);
//print_vector("Mediciones esperadas", filtro->meas_counter, filtro->exp_meas);
//print_vector("Innovacion", filtro->meas_counter, filtro->aux8);
accum_vec(N_STATES, filtro->states, filtro->aux8); // estado = estado + G(z_medido - z_esperado)
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
mulmat(N_STATES, filtro->meas_counter, N_STATES, filtro->G, filtro->H, filtro->aux9); //  aux9 = G * H
mulmat(N_STATES, N_STATES, N_STATES, filtro->aux9, filtro->cov, filtro->aux7); // aux7 = G * H * cov_priori
mat_negate(N_STATES, N_STATES, filtro->aux7); // aux7 = - G * H * cov_priori
accum(N_STATES, N_STATES, filtro->cov, filtro->aux7); // cov = cov_priori - G * H * cov_priori
//print_gain("cov: ", N_STATES, N_STATES,  filtro->cov);
if(filtro->states[2] < MIN_HEIGHT){
	filtro->states[2] = MIN_HEIGHT;
}
else if (filtro->states[2] > MAX_HEIGHT){
	filtro->states[2] = MAX_HEIGHT;
}
(*filtro).beta = 0;
(*filtro).gamma = 0;
}

void IMU_states(ofs_ekf_t* filtro){

// partial_ax / partial_q
filtro->H[0][N_P + N_V + 0] = -2 * (-g) * (filtro->q.q3);     
filtro->H[0][N_P + N_V + 1] = 2 * (-g) * (filtro->q.q4);     
filtro->H[0][N_P + N_V + 2] = -2 * (-g) * (filtro->q.q1);     
filtro->H[0][N_P + N_V + 3] = 2 * (-g) * (filtro->q.q2);  
// partial_ay / partial_q   
filtro->H[1][N_P + N_V + 0] = 2 * (-g) * (filtro->q.q2);      
filtro->H[1][N_P + N_V + 1] = 2 * (-g) * (filtro->q.q1);      
filtro->H[1][N_P + N_V + 2] = 2 * (-g) * (filtro->q.q4);     
filtro->H[1][N_P + N_V + 3] = 2 * (-g) * (filtro->q.q3);
// partial_az / partial_q  
filtro->H[2][N_P + N_V + 0] = 2 * (-g) * (filtro->q.q1);      
filtro->H[2][N_P + N_V + 1] = -2 * (-g) * (filtro->q.q2);      
filtro->H[2][N_P + N_V + 2] = -2 * (-g) * (filtro->q.q3);     
filtro->H[2][N_P + N_V + 3] = 2 * (-g) * (filtro->q.q4);
}

void MAG_states(ofs_ekf_t* filtro){

// partial_magx / partial_q
filtro->H[N_IMU][N_P + N_V + 0] = 2 * filtro->m[0] * filtro->q.q1 + 2 * filtro->m[1] * filtro->q.q4 - 2 * filtro->m[2] * filtro->q.q3;     
filtro->H[N_IMU][N_P + N_V + 1] = 2 * filtro->m[0] * filtro->q.q2 + 2 * filtro->m[1] * filtro->q.q3 + 2 * filtro->m[2] * filtro->q.q4;     
filtro->H[N_IMU][N_P + N_V + 2] = -2 * filtro->m[0] * filtro->q.q3 + 2 * filtro->m[1] * filtro->q.q2 - 2 * filtro->m[2] * filtro->q.q1;     
filtro->H[N_IMU][N_P + N_V + 3] = -2 * filtro->m[0] * filtro->q.q4 + 2 * filtro->m[1] * filtro->q.q1 + 2 * filtro->m[2] * filtro->q.q2;  
// partial_magy / partial_q   
filtro->H[N_IMU + 1][N_P + N_V + 0] = -2 * filtro->m[0] * filtro->q.q4 + 2 * filtro->m[1] * filtro->q.q1 + 2 * filtro->m[2] * filtro->q.q2;      
filtro->H[N_IMU + 1][N_P + N_V + 1] = 2 * filtro->m[0] * filtro->q.q3 - 2 * filtro->m[1] * filtro->q.q2 + 2 * filtro->m[2] * filtro->q.q1;      
filtro->H[N_IMU + 1][N_P + N_V + 2] = 2 * filtro->m[0] * filtro->q.q2 + 2 * filtro->m[1] * filtro->q.q3 + 2 * filtro->m[2] * filtro->q.q4;     
filtro->H[N_IMU + 1][N_P + N_V + 3] = -2 * filtro->m[0] * filtro->q.q1 - 2 * filtro->m[1] * filtro->q.q4 + 2 * filtro->m[2] * filtro->q.q3;
// partial_magz / partial_q  
filtro->H[N_IMU + 2][N_P + N_V + 0] = 2 * filtro->m[0] * filtro->q.q3 - 2 * filtro->m[1] * filtro->q.q2 + 2 * filtro->m[2] * filtro->q.q1;      
filtro->H[N_IMU + 2][N_P + N_V + 1] = 2 * filtro->m[0] * filtro->q.q4 - 2 * filtro->m[1] * filtro->q.q1 - 2 * filtro->m[2] * filtro->q.q2;      
filtro->H[N_IMU + 2][N_P + N_V + 2] = 2 * filtro->m[0] * filtro->q.q1 + 2 * filtro->m[1] * filtro->q.q4 - 2 * filtro->m[2] * filtro->q.q3;     
filtro->H[N_IMU + 2][N_P + N_V + 3] = 2 * filtro->m[0] * filtro->q.q2 + 2 * filtro->m[1] * filtro->q.q3 + 2 * filtro->m[2] * filtro->q.q4;
}

void OFS_states(ofs_ekf_t* filtro, mediciones_t *z){
// partial_nx / partial_pz
filtro->H[N_IMU + N_MAG][2] = ((*z).tau * filtro->f * (2 * pow(filtro->q.q1, 2) + 2 * pow(filtro->q.q4, 2) - 1)) * \
                        (filtro->q.q1 * (filtro->q.q1 * filtro->v[0] - filtro->q.q3 * filtro->v[2] + \
                        filtro->q.q4 * filtro->v[1]) - filtro->q.q2 * (-filtro->q.q2 * filtro->v[0] - \
                        filtro->q.q3 * filtro->v[1] - filtro->q.q4 * filtro->v[2]) \
                        -filtro->q.q3 * (filtro->q.q1 * filtro->v[2] - filtro->q.q2 * filtro->v[1] + \
                        filtro->q.q3 * filtro->v[0]) + filtro->q.q4 * (filtro->q.q1 * filtro->v[1] + \
                        filtro->q.q2 * filtro->v[2] - filtro->q.q4 * filtro->v[0]) \
                        )/ pow(filtro->p[2], 2);     
// partial_nx / partial_v   
filtro->H[N_IMU + N_MAG][N_P + 0] = - ((*z).tau * filtro->f * (2 * pow(filtro->q.q1, 2) + 2 * pow(filtro->q.q4, 2) - 1)) * \
                         (pow(filtro->q.q1, 2) + pow(filtro->q.q2, 2) - pow(filtro->q.q3, 2) - pow(filtro->q.q4, 2)) / filtro->p[2];       
filtro->H[N_IMU + N_MAG][N_P + 1] = - ((*z).tau * filtro->f * (2 *filtro->q.q1 *filtro->q.q4 + 2 * filtro->q.q2 * filtro->q.q3)) * \
                         (2 * pow(filtro->q.q1, 2) + 2 * pow(filtro->q.q4, 2) - 1) / filtro->p[2];        
filtro->H[N_IMU + N_MAG][N_P + 2] = - ((*z).tau * filtro->f * (-2 *filtro->q.q1 *filtro->q.q3 + 2 * filtro->q.q2 * filtro->q.q4)) * \
                         (2 * pow(filtro->q.q1, 2) + 2 * pow(filtro->q.q4, 2) - 1) / filtro->p[2];      
// partial_nx / partial_q  
filtro->H[N_IMU + N_MAG][N_P + N_V + 0] = (*z).tau * filtro->f * (-4* filtro->q.q1 * (filtro->q.q1 * (filtro->q.q1 * filtro->v[0] -\
                        filtro->q.q3 * filtro->v[2] +filtro->q.q4 * filtro->v[1]) \
                        +filtro->q.q2 * (filtro->q.q2 * filtro->v[0] +filtro->q.q3 * filtro->v[1] +filtro->q.q4 * filtro->v[2])\
                        -filtro->q.q3 * (filtro->q.q1 * filtro->v[2] -filtro->q.q2 * filtro->v[1] +filtro->q.q3 * filtro->v[0]) \
                        +filtro->q.q4 * (filtro->q.q1 * filtro->v[1] +filtro->q.q2 * filtro->v[2] -filtro->q.q4 * filtro->v[0])) \
                        -2* (2 * pow(filtro->q.q1, 2) + 2 * pow(filtro->q.q4, 2) - 1) \
                        * (filtro->q.q1 * filtro->v[0] -filtro->q.q3 * filtro->v[2] +filtro->q.q4 * filtro->v[1])) / filtro->p[2];    
filtro->H[N_IMU + N_MAG][N_P + N_V + 1] = - 2 * (*z).tau * filtro->f * (2 * pow(filtro->q.q1, 2) + 2 * pow(filtro->q.q4, 2) - 1) \
                         * (filtro->q.q2 * filtro->v[0] +filtro->q.q3 * filtro->v[1] +filtro->q.q4 * filtro->v[2]) / filtro->p[2];    
filtro->H[N_IMU + N_MAG][N_P + N_V + 2] = 2 * (*z).tau * filtro->f * (2 * pow(filtro->q.q1, 2) + 2 * pow(filtro->q.q4, 2) - 1) \
                         * (filtro->q.q1 * filtro->v[2] -filtro->q.q2 * filtro->v[1] +filtro->q.q3 * filtro->v[0]) / filtro->p[2];   
filtro->H[N_IMU + N_MAG][N_P + N_V + 3] = (*z).tau * filtro->f * (-4 *filtro->q.q4 * (filtro->q.q1 * (filtro->q.q1 * filtro->v[0] \
                        - filtro->q.q3 * filtro->v[2] +filtro->q.q4 * filtro->v[1]) +filtro->q.q2 * (filtro->q.q2 * filtro->v[0] +\
                        filtro->q.q3 * filtro->v[1] +filtro->q.q4 * filtro->v[2]) -filtro->q.q3 * (filtro->q.q1 * filtro->v[2] \
                        - filtro->q.q2 * filtro->v[1] +filtro->q.q3 * filtro->v[0]) +filtro->q.q4 * (filtro->q.q1 * filtro->v[1] \
                        + filtro->q.q2 * filtro->v[2] -filtro->q.q4 * filtro->v[0])) -2* (2 * pow(filtro->q.q1, 2) + 2 * pow(filtro->q.q4, 2) - 1) \
                        * (filtro->q.q1 * filtro->v[1] +filtro->q.q2 * filtro->v[2] -filtro->q.q4 * filtro->v[0])) / filtro->p[2];  
// partial_ny / partial_pz
filtro->H[N_IMU + N_MAG + 1][2] = ((*z).tau * filtro->f * (2 * pow(filtro->q.q1, 2) + 2 * pow(filtro->q.q4, 2) - 1)) * \
                        (filtro->q.q1 * (filtro->q.q1 * filtro->v[1] +filtro->q.q2 * filtro->v[2] -filtro->q.q4 * filtro->v[0]) \
                        +filtro->q.q2 * (filtro->q.q1 * filtro->v[2] -filtro->q.q2 * filtro->v[1] +filtro->q.q3 * filtro->v[0]) \
                        + filtro->q.q3 * (filtro->q.q2 * filtro->v[0] +filtro->q.q3 * filtro->v[1] +filtro->q.q4 * filtro->v[2]) \
                        -filtro->q.q4 * (filtro->q.q1 * filtro->v[0] -filtro->q.q3 * filtro->v[2] +filtro->q.q4 * filtro->v[1]) \
                        ) / pow(filtro->p[2], 2);    
// partial_ny / partial_v   
filtro->H[N_IMU + N_MAG + 1][N_P + 0] = ((*z).tau * filtro->f * (2 * filtro->q.q1 *filtro->q.q4 - 2 * filtro->q.q2 * filtro->q.q3)) * \
                         (2 * pow(filtro->q.q1, 2) + 2 * pow(filtro->q.q4, 2) - 1) / filtro->p[2];           
filtro->H[N_IMU + N_MAG + 1][N_P + 1] = - ((*z).tau * filtro->f * (2 * pow(filtro->q.q1, 2) + 2 * pow(filtro->q.q4, 2) - 1)) * \
                         (pow(filtro->q.q1, 2) - pow(filtro->q.q2, 2) + pow(filtro->q.q3, 2) - pow(filtro->q.q4, 2)) / filtro->p[2];       
filtro->H[N_IMU + N_MAG + 1][N_P + 2] = - ((*z).tau * filtro->f * (2 *filtro->q.q1 *filtro->q.q2 + 2 * filtro->q.q3 * filtro->q.q4)) * \
                         (2 * pow(filtro->q.q1, 2) + 2 * pow(filtro->q.q4, 2) - 1) / filtro->p[2];   
// partial_ny / partial_q  
filtro->H[N_IMU + N_MAG + 1][N_P + N_V + 0] = - 2 * (*z).tau * filtro->f * (2 *filtro->q.q1 * (filtro->q.q1 * (filtro->q.q1 * filtro->v[1] \
                           +filtro->q.q2 * filtro->v[2] -filtro->q.q4 * filtro->v[0]) \
                         +filtro->q.q2 * (filtro->q.q1 * filtro->v[2] -filtro->q.q2 * filtro->v[1] +filtro->q.q3 * filtro->v[0]) \
                         +filtro->q.q3 * (filtro->q.q2 * filtro->v[0] +filtro->q.q3 * filtro->v[1] +filtro->q.q4 * filtro->v[2]) \
                         -filtro->q.q4 * (filtro->q.q1 * filtro->v[0] -filtro->q.q3 * filtro->v[2] +filtro->q.q4 * filtro->v[1])) \
                         + (2 * pow(filtro->q.q1, 2) + 2 * pow(filtro->q.q4, 2) - 1) \
                         * (filtro->q.q1 * filtro->v[1] +filtro->q.q2 * filtro->v[2] -filtro->q.q4 * filtro->v[0])) / filtro->p[2];     
filtro->H[N_IMU + N_MAG + 1][N_P + N_V + 1] = - 2 * (*z).tau * filtro->f * (2 * pow(filtro->q.q1, 2) + 2 * pow(filtro->q.q4, 2) - 1) \
                         * (filtro->q.q1 * filtro->v[2] -filtro->q.q2 * filtro->v[1] +filtro->q.q3 * filtro->v[0]) / filtro->p[2]; 
filtro->H[N_IMU + N_MAG + 1][N_P + N_V + 2] = - 2 * (*z).tau * filtro->f * (2 * pow(filtro->q.q1, 2) + 2 * pow(filtro->q.q4, 2) - 1) \
                         * (filtro->q.q2 * filtro->v[0] +filtro->q.q3 * filtro->v[1] +filtro->q.q4 * filtro->v[2]) / filtro->p[2]; 
filtro->H[N_IMU + N_MAG + 1][N_P + N_V + 3] = - 2 * (*z).tau * filtro->f * (2 * filtro->q.q4 * (filtro->q.q1 * (filtro->q.q1 * filtro->v[1] \
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

double calc_trace_cov(ofs_ekf_t *filtro) {
    double suma = 0;
    for (int i = 0; i < N_STATES; i++) {
        suma += filtro->cov[i][i];
    }
    return suma;
}