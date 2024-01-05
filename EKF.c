#include "ekf.h"
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
    /* Inicializo el campo magnetico respecto a mundo */
    filtro->m[0] = m[0];
    filtro->m[1] = m[1];
    filtro->m[2] = m[2];
    filtro->qm = vec2quat(filtro->m);
    /* Inicializo el valor de la gravedad respecto a mundo */
    filtro->qg.q1 = 0;
    filtro->qg.q2 = 0;
    filtro->qg.q3 = 0;
    filtro->qg.q4 = -g;
    /* Inicializo estados */
    vec_zeros(N_STATES, filtro->states); //p, v, q
    filtro->states[N_P + N_V] = 1; // Los ejes del vehiculo comienzan alineados respecto al mundo //TODO: ESTO CAMBIA CON EL MAGNETOMETRO 
    /* Inicializo covarianza */
    mat_zeros(N_STATES, N_STATES, filtro->cov); // Matriz de covarianza de estados
    mat_addeye(N_STATES, filtro->cov);
    matmul_scalar(N_STATES, N_STATES, filtro->cov, 10); 
    /* Inicializo F */
    mat_zeros(N_STATES, N_STATES, filtro->F);
    mat_addeye(N_STATES, filtro->F);
    /* Inicializo W */
    double W[N_STATES][N_PROC_NOISE]; // Derivada de vector de estados respecto de ruidos
    mat_zeros(N_STATES, N_PROC_NOISE, W);
    // d(velocidad) / d(uax)
    W[N_P][0] = 1;
    // d(velocidad) / d(uay)
    W[N_P+1][1] = 1;
    // d(velocidad) / d(uaz)
    W[N_P+2][2] = 1;
    // d(q1) / d(uw)
    W[N_P+N_V][3] = 1;
    // d(q2) / d(uw)
    W[N_P+N_V+1][4] = 1;
    // d(q3) / d(uw)
    W[N_P+N_V+2][5] = 1;
    // d(q4) / d(uw)
    W[N_P+N_V+3][6] = 1;
    /* Inicializo Wt */
    double Wt[N_PROC_NOISE][N_STATES]; // Derivada de vector de estados respecto de ruidos
    transpose(N_STATES, N_PROC_NOISE, W, Wt);
    /* Inicializo Q */
    double Q[N_PROC_NOISE][N_PROC_NOISE]; // Ruido de proceso
    mat_zeros(N_PROC_NOISE, N_PROC_NOISE, Q);
    Q[0][0] = U_P_A; //uax
    Q[1][1] = U_P_A; //uay
    Q[2][2] = U_P_A; //uaz
    Q[3][3] = U_P_W; //uw
    Q[4][4] = U_P_W; //uw
    Q[5][5] = U_P_W; //uw
    Q[6][6] = U_P_W; //uw
    /* Calculo W Q Wt dado que no varia segun la iteracion */
    mulmat(N_PROC_NOISE, N_PROC_NOISE, N_STATES, Q, Wt, filtro->aux); // aux4 = Q * Wt
    mulmat(N_STATES, N_PROC_NOISE, N_STATES, W, filtro->aux, filtro->WQWt); //  W * Q * Wt
    /* Inicializo R */
    mat_zeros(N_OBS, N_OBS, filtro->R);
    filtro->R[0][0] = U_A; //uax
    filtro->R[1][1] = U_A; //uay
    filtro->R[2][2] = U_A; //uaz
    filtro->R[3][3] = U_MAG; //umagx
    filtro->R[4][4] = U_MAG; //umagy
    filtro->R[5][5] = U_MAG; //umagz
    filtro->R[6][6] = U_FLOW; //u flowx    
    filtro->R[7][7] = U_FLOW; //u flowy  
    filtro->R[8][8] = U_RANGE; // distance
    /* Limpio H */
    mat_zeros(N_OBS, N_STATES, filtro->H); // Jacobiano de mediciones respecto a estados
    /* Calculo factor de conversion del flujo optico */
    filtro->f = 35.0 / (2 * atan2f(4.2 * M_PI / 180, 2));  // Factor de conversión
}

void prediction_step(ofs_ekf_t* filtro, mediciones_t u){
    /* Asigno valores a las variables auxiliares */
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

    /* Calculo Fk */
    // F tiene la diagonal en unos, se asigna al inicializar el filtro
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
    filtro->F[N_P+N_V+1][N_P+N_V] = u.wx * u.dt / 2;
    filtro->F[N_P+N_V+2][N_P+N_V] = u.wy * u.dt / 2;
    filtro->F[N_P+N_V+3][N_P+N_V] = u.wz * u.dt / 2;
    //d(q) / d(q2)
    filtro->F[N_P+N_V][N_P+N_V+1] = -u.wx * u.dt / 2;
    filtro->F[N_P+N_V+2][N_P+N_V+1] = -u.wz * u.dt / 2;
    filtro->F[N_P+N_V+3][N_P+N_V+1] = u.wy * u.dt / 2;
    //d(q) / d(q3)
    filtro->F[N_P+N_V][N_P+N_V+2] = -u.wy * u.dt / 2;
    filtro->F[N_P+N_V+1][N_P+N_V+2] = u.wz * u.dt / 2;
    filtro->F[N_P+N_V+3][N_P+N_V+2] = -u.wx * u.dt / 2;
    //d(q) / d(q4)
    filtro->F[N_P+N_V][N_P+N_V+3] = -u.wz * u.dt / 2;
    filtro->F[N_P+N_V+1][N_P+N_V+3] = -u.wy * u.dt / 2;
    filtro->F[N_P+N_V+2][N_P+N_V+3] = u.wx * u.dt / 2;

    /* Calculo de Covarianza a priori */
    //Como W y Q no son variables, se calculan al inicializar el filtro por unica vez -> WQWt
    transpose(N_STATES, N_STATES, filtro->F, filtro->aux2); //aux2 = Ft
    mulmat(N_STATES, N_STATES, N_STATES, filtro->cov, filtro->aux2, filtro->aux3); // aux3 = cov * Ft
    mulmat(N_STATES, N_STATES, N_STATES, filtro->F, filtro->aux3, filtro->aux2); // aux2 = F * cov * Ft
    add(N_STATES, N_STATES, filtro->aux2, filtro->WQWt, filtro->cov); // cov = F * cov * Ft +  W * Q * Wt

    /* Calculo de estados */
    // Posicion
    vecmul_scalar2(N_V, filtro->v, filtro->temp, u.dt); // temp = Vk = dt * Vk
    add_vec(N_P, filtro->p, filtro->temp, filtro->states); // Pk+1 = Pk + dt * Vk
    if(filtro->states[2] < MIN_HEIGHT){
    	filtro->states[2] = MIN_HEIGHT;
    }
    else if (filtro->states[2] > MAX_HEIGHT){
    	filtro->states[2] = MAX_HEIGHT;
    }
    // Velocidad
    // Queremos pasar la acel de cuerpo a mundo
    filtro->prov = quat_mult(filtro->qa_meas, quat_conjugate(filtro->q)); // prov = (filtro->qa_meas) * -q
    filtro->prov2 = quat_mult(filtro->q, filtro->prov); // prov2 = q * (filtro->qa_meas) * -q
    quat_sub(&filtro->prov, filtro->prov2, filtro->qg); // prov = q * (filtro->qa_meas) * -q - g
    quat2vec(filtro->prov, filtro->temp); //temp = prov
    vecmul_scalar(N_V, filtro->temp, u.dt); // temp = dt (a_meas - g)
    add_vec(N_V, filtro->v, filtro->temp, &filtro->states[N_P]); // Vk+1 = Vk + dt * (a - g)
    // Cuaterniones
    filtro->prov = quat_mult(filtro->q, filtro->qw_meas); // prov = qk * qw
    quat_scalar(&filtro->prov, u.dt / 2); // prov = (dt / 2) qk * qw
    quat_add(&filtro->prov2, filtro->prov, filtro->q); // q = qk + (dt / 2) qk * qw
    quat_Normalization(&filtro->prov2);
    filtro->states[N_P + N_V] = filtro->prov2.q1;
    filtro->states[N_P + N_V + 1] = filtro->prov2.q2;
    filtro->states[N_P + N_V + 2] = filtro->prov2.q3;
    filtro->states[N_P + N_V + 3] = filtro->prov2.q4;

}

void correction_step(ofs_ekf_t* filtro, mediciones_t* z){

    /* Asigno variables auxiliares */
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

    /* Calculo del jacobiano H */
    /*************** IMU ***************/
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

    /*************** MAG ***************/
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

    /*************** FLUJO OPTICO ***************/
    // partial_nx / partial_pz
    filtro->H[N_IMU + N_MAG][2] = ((*z).dt * filtro->f * (2 * pow(filtro->q.q1, 2) + 2 * pow(filtro->q.q4, 2) - 1)) * \
                            (filtro->q.q1 * (filtro->q.q1 * filtro->v[0] - filtro->q.q3 * filtro->v[2] + \
                            filtro->q.q4 * filtro->v[1]) - filtro->q.q2 * (-filtro->q.q2 * filtro->v[0] - \
                            filtro->q.q3 * filtro->v[1] - filtro->q.q4 * filtro->v[2]) \
                            -filtro->q.q3 * (filtro->q.q1 * filtro->v[2] - filtro->q.q2 * filtro->v[1] + \
                            filtro->q.q3 * filtro->v[0]) + filtro->q.q4 * (filtro->q.q1 * filtro->v[1] + \
                            filtro->q.q2 * filtro->v[2] - filtro->q.q4 * filtro->v[0]) \
                            )/ pow(filtro->p[2], 2);     
    // partial_nx / partial_v   
    filtro->H[N_IMU + N_MAG][N_P + 0] = - ((*z).dt * filtro->f * (2 * pow(filtro->q.q1, 2) + 2 * pow(filtro->q.q4, 2) - 1)) * \
                            (pow(filtro->q.q1, 2) + pow(filtro->q.q2, 2) - pow(filtro->q.q3, 2) - pow(filtro->q.q4, 2)) / filtro->p[2];       
    filtro->H[N_IMU + N_MAG][N_P + 1] = - ((*z).dt * filtro->f * (2 *filtro->q.q1 *filtro->q.q4 + 2 * filtro->q.q2 * filtro->q.q3)) * \
                            (2 * pow(filtro->q.q1, 2) + 2 * pow(filtro->q.q4, 2) - 1) / filtro->p[2];        
    filtro->H[N_IMU + N_MAG][N_P + 2] = - ((*z).dt * filtro->f * (-2 *filtro->q.q1 *filtro->q.q3 + 2 * filtro->q.q2 * filtro->q.q4)) * \
                            (2 * pow(filtro->q.q1, 2) + 2 * pow(filtro->q.q4, 2) - 1) / filtro->p[2];      
    // partial_nx / partial_q  
    filtro->H[N_IMU + N_MAG][N_P + N_V + 0] = (*z).dt * filtro->f * (-4* filtro->q.q1 * (filtro->q.q1 * (filtro->q.q1 * filtro->v[0] -\
                            filtro->q.q3 * filtro->v[2] +filtro->q.q4 * filtro->v[1]) \
                            +filtro->q.q2 * (filtro->q.q2 * filtro->v[0] +filtro->q.q3 * filtro->v[1] +filtro->q.q4 * filtro->v[2])\
                            -filtro->q.q3 * (filtro->q.q1 * filtro->v[2] -filtro->q.q2 * filtro->v[1] +filtro->q.q3 * filtro->v[0]) \
                            +filtro->q.q4 * (filtro->q.q1 * filtro->v[1] +filtro->q.q2 * filtro->v[2] -filtro->q.q4 * filtro->v[0])) \
                            -2* (2 * pow(filtro->q.q1, 2) + 2 * pow(filtro->q.q4, 2) - 1) \
                            * (filtro->q.q1 * filtro->v[0] -filtro->q.q3 * filtro->v[2] +filtro->q.q4 * filtro->v[1])) / filtro->p[2];    
    filtro->H[N_IMU + N_MAG][N_P + N_V + 1] = - 2 * (*z).dt * filtro->f * (2 * pow(filtro->q.q1, 2) + 2 * pow(filtro->q.q4, 2) - 1) \
                            * (filtro->q.q2 * filtro->v[0] +filtro->q.q3 * filtro->v[1] +filtro->q.q4 * filtro->v[2]) / filtro->p[2];    
    filtro->H[N_IMU + N_MAG][N_P + N_V + 2] = 2 * (*z).dt * filtro->f * (2 * pow(filtro->q.q1, 2) + 2 * pow(filtro->q.q4, 2) - 1) \
                            * (filtro->q.q1 * filtro->v[2] -filtro->q.q2 * filtro->v[1] +filtro->q.q3 * filtro->v[0]) / filtro->p[2];   
    filtro->H[N_IMU + N_MAG][N_P + N_V + 3] = (*z).dt * filtro->f * (-4 *filtro->q.q4 * (filtro->q.q1 * (filtro->q.q1 * filtro->v[0] \
                            - filtro->q.q3 * filtro->v[2] +filtro->q.q4 * filtro->v[1]) +filtro->q.q2 * (filtro->q.q2 * filtro->v[0] +\
                            filtro->q.q3 * filtro->v[1] +filtro->q.q4 * filtro->v[2]) -filtro->q.q3 * (filtro->q.q1 * filtro->v[2] \
                            - filtro->q.q2 * filtro->v[1] +filtro->q.q3 * filtro->v[0]) +filtro->q.q4 * (filtro->q.q1 * filtro->v[1] \
                            + filtro->q.q2 * filtro->v[2] -filtro->q.q4 * filtro->v[0])) -2* (2 * pow(filtro->q.q1, 2) + 2 * pow(filtro->q.q4, 2) - 1) \
                            * (filtro->q.q1 * filtro->v[1] +filtro->q.q2 * filtro->v[2] -filtro->q.q4 * filtro->v[0])) / filtro->p[2];  
    // partial_ny / partial_pz
    filtro->H[N_IMU + N_MAG + 1][2] = ((*z).dt * filtro->f * (2 * pow(filtro->q.q1, 2) + 2 * pow(filtro->q.q4, 2) - 1)) * \
                            (filtro->q.q1 * (filtro->q.q1 * filtro->v[1] +filtro->q.q2 * filtro->v[2] -filtro->q.q4 * filtro->v[0]) \
                            +filtro->q.q2 * (filtro->q.q1 * filtro->v[2] -filtro->q.q2 * filtro->v[1] +filtro->q.q3 * filtro->v[0]) \
                            + filtro->q.q3 * (filtro->q.q2 * filtro->v[0] +filtro->q.q3 * filtro->v[1] +filtro->q.q4 * filtro->v[2]) \
                            -filtro->q.q4 * (filtro->q.q1 * filtro->v[0] -filtro->q.q3 * filtro->v[2] +filtro->q.q4 * filtro->v[1]) \
                            ) / pow(filtro->p[2], 2);    
    // partial_ny / partial_v   
    filtro->H[N_IMU + N_MAG + 1][N_P + 0] = ((*z).dt * filtro->f * (2 * filtro->q.q1 *filtro->q.q4 - 2 * filtro->q.q2 * filtro->q.q3)) * \
                            (2 * pow(filtro->q.q1, 2) + 2 * pow(filtro->q.q4, 2) - 1) / filtro->p[2];           
    filtro->H[N_IMU + N_MAG + 1][N_P + 1] = - ((*z).dt * filtro->f * (2 * pow(filtro->q.q1, 2) + 2 * pow(filtro->q.q4, 2) - 1)) * \
                            (pow(filtro->q.q1, 2) - pow(filtro->q.q2, 2) + pow(filtro->q.q3, 2) - pow(filtro->q.q4, 2)) / filtro->p[2];       
    filtro->H[N_IMU + N_MAG + 1][N_P + 2] = - ((*z).dt * filtro->f * (2 *filtro->q.q1 *filtro->q.q2 + 2 * filtro->q.q3 * filtro->q.q4)) * \
                            (2 * pow(filtro->q.q1, 2) + 2 * pow(filtro->q.q4, 2) - 1) / filtro->p[2];   
    // partial_ny / partial_q  
    filtro->H[N_IMU + N_MAG + 1][N_P + N_V + 0] = - 2 * (*z).dt * filtro->f * (2 *filtro->q.q1 * (filtro->q.q1 * (filtro->q.q1 * filtro->v[1] \
                            +filtro->q.q2 * filtro->v[2] -filtro->q.q4 * filtro->v[0]) \
                            +filtro->q.q2 * (filtro->q.q1 * filtro->v[2] -filtro->q.q2 * filtro->v[1] +filtro->q.q3 * filtro->v[0]) \
                            +filtro->q.q3 * (filtro->q.q2 * filtro->v[0] +filtro->q.q3 * filtro->v[1] +filtro->q.q4 * filtro->v[2]) \
                            -filtro->q.q4 * (filtro->q.q1 * filtro->v[0] -filtro->q.q3 * filtro->v[2] +filtro->q.q4 * filtro->v[1])) \
                            + (2 * pow(filtro->q.q1, 2) + 2 * pow(filtro->q.q4, 2) - 1) \
                            * (filtro->q.q1 * filtro->v[1] +filtro->q.q2 * filtro->v[2] -filtro->q.q4 * filtro->v[0])) / filtro->p[2];     
    filtro->H[N_IMU + N_MAG + 1][N_P + N_V + 1] = - 2 * (*z).dt * filtro->f * (2 * pow(filtro->q.q1, 2) + 2 * pow(filtro->q.q4, 2) - 1) \
                            * (filtro->q.q1 * filtro->v[2] -filtro->q.q2 * filtro->v[1] +filtro->q.q3 * filtro->v[0]) / filtro->p[2]; 
    filtro->H[N_IMU + N_MAG + 1][N_P + N_V + 2] = - 2 * (*z).dt * filtro->f * (2 * pow(filtro->q.q1, 2) + 2 * pow(filtro->q.q4, 2) - 1) \
                            * (filtro->q.q2 * filtro->v[0] +filtro->q.q3 * filtro->v[1] +filtro->q.q4 * filtro->v[2]) / filtro->p[2]; 
    filtro->H[N_IMU + N_MAG + 1][N_P + N_V + 3] = - 2 * (*z).dt * filtro->f * (2 * filtro->q.q4 * (filtro->q.q1 * (filtro->q.q1 * filtro->v[1] \
                                +filtro->q.q2 * filtro->v[2] -filtro->q.q4 * filtro->v[0]) \
                            +filtro->q.q2 * (filtro->q.q1 * filtro->v[2] -filtro->q.q2 * filtro->v[1] +filtro->q.q3 * filtro->v[0]) \
                            +filtro->q.q3 * (filtro->q.q2 * filtro->v[0] +filtro->q.q3 * filtro->v[1] +filtro->q.q4 * filtro->v[2]) \
                            -filtro->q.q4 * (filtro->q.q1 * filtro->v[0] -filtro->q.q3 * filtro->v[2] +filtro->q.q4 * filtro->v[1])) \
                            - (2 * pow(filtro->q.q1, 2) + 2 * pow(filtro->q.q4, 2) - 1) \
                            * (filtro->q.q1 * filtro->v[0] -filtro->q.q3 * filtro->v[2] +filtro->q.q4 * filtro->v[1])) / filtro->p[2]; 
    
    /*************** DISTANCIA ***************/
    // partial_tofs / partial_pz
    filtro->H[N_IMU + N_MAG + N_OFS + 0][2] = 1/(2 * (pow(filtro->q.q1, 2) + pow(filtro->q.q4, 2)) - 1); 
    // partial_tofs / partial_q1  
    filtro->H[N_IMU + N_MAG + N_OFS + 0][N_P + N_V] = -4 * filtro->p[2] * filtro->q.q1 * pow(filtro->H[N_IMU + N_MAG + N_OFS + 0][2], 2);
    // partial_tofs / partial_q4
    filtro->H[N_IMU + N_MAG + N_OFS + 0][N_P + N_V + 3] = -4 * filtro->p[2] * filtro->q.q4 * pow(filtro->H[N_IMU + N_MAG + N_OFS + 0][2], 2);

    /* Mediciones */
    filtro->meas[0] = z->ax;
    filtro->meas[1] = z->ay;
    filtro->meas[2] = z->az;
    filtro->meas[N_IMU + 0] = z->mx;
    filtro->meas[N_IMU + 1] = z->my;
    filtro->meas[N_IMU + 2] = z->mz;
    filtro->meas[N_IMU + N_MAG + 0] = z->ofx;
    filtro->meas[N_IMU + N_MAG + 1] = z->ofy;
    filtro->meas[N_IMU + N_MAG + N_OFS] = z->range;

    /* Mediciones predichas */
    // Aceleracion
    filtro->prov = quat_mult(filtro->qg, filtro->q); // prov = qg * qk (de mundo a cuerpo)
    filtro->prov2 = quat_mult(quat_conjugate(filtro->q), filtro->prov); // prov = q- * qg * q (de mundo a cuerpo)
    quat2vec(filtro->prov2, filtro->exp_meas); // measurements = q- * qg * q
    //Campo magnetico
    filtro->prov = quat_mult(filtro->qm, filtro->q); // prov = qm * qk (de mundo a cuerpo)
    filtro->prov2 = quat_mult(quat_conjugate(filtro->q), filtro->prov); // prov = q- * qm * q (de mundo a cuerpo)
    quat2vec(filtro->prov2, &filtro->exp_meas[N_IMU]); // measurements = q- * qm * q
    //Flujo optico
    filtro->prov = vec2quat(&filtro->states[N_P]); //Paso la velocidad respecto al mundo a quat para convertirlo a body
    filtro->prov2 = quat_mult(filtro->prov, filtro->q); // prov = v * q
    filtro->prov = quat_mult(quat_conjugate(filtro->q), filtro->prov2); // prov = q- * v * q
    filtro->exp_meas[N_IMU + N_MAG + 0] = -(*z).dt * filtro->f * (filtro->prov.q2 * (2 * (pow(filtro->q.q1, 2) + pow(filtro->q.q4, 2)) - 1) / \
                                         filtro->p[2]  + (*z).wy);
    filtro->exp_meas[N_IMU + N_MAG + 1] = -(*z).dt * filtro->f * (filtro->prov.q3 * (2 * (pow(filtro->q.q1, 2) + pow(filtro->q.q4, 2)) - 1) / \
                                         filtro->p[2]  - (*z).wx);  
    //Distancia
    filtro->exp_meas[N_IMU + N_MAG + N_OFS] = filtro->p[2] / (2 * (pow(filtro->q.q1, 2) + pow(filtro->q.q4, 2)) - 1);

    /* Paso secuencial para actualizacion de estados y covarianza */
    for (int i = 0; i < N_OBS; i++){
        mat_copy(N_STATES, N_STATES, filtro->cov, filtro->aux4); //aux4 = cov = Pi
        mat_getrow(N_OBS, N_STATES, filtro->H, i, filtro->temp3); // temp3 = H[i][:]
        mulvec(N_STATES, N_STATES, filtro->aux4, filtro->temp3, filtro->temp4); // temp4 = Pi * H[i][:].t
        vec_dot(N_STATES, filtro->temp3, filtro->temp4, &filtro->S); // S = H[i][:] * Pi * H[i][:].t
        filtro->S += filtro->R[i][i]; // S = H[i][:] * Pi * H[i][:].t + Ri
        vecmul_scalar(N_STATES, filtro->temp4, 1 / filtro->S); //Ki = Pi * H[i][:].t / (H[i][:] * Pi * H[i][:].t + Ri)
        // Correccion de estados
        vecmul_scalar2(N_STATES, filtro->temp4, filtro->temp5, filtro->meas[i] - filtro->exp_meas[i]); // temp5 = Ki * (y[i] - y_est[i])
        accum_vec(N_STATES, filtro->states, filtro->temp5); // states = states + Ki * (y[i] - y_est[i])
        // Correccion de covarianza
        mat_zeros(N_STATES, N_STATES, filtro->aux3); // aux3 = 0
        mat_addeye(N_STATES, filtro->aux3); // aux3 = I
        vec_outer(N_STATES, filtro->temp4, filtro->temp3, filtro->aux2); // aux2 = Ki * H[i][:]
        mat_negate(N_STATES, N_STATES, filtro->aux2); // aux2 = -Ki * H[i][:]
        accum(N_STATES, N_STATES, filtro->aux3, filtro->aux2); //aux3 = I - Ki * H[i][:]
        //printf("Aca\n");
        //print_gain("cov", N_STATES, N_STATES, filtro->cov);
        //print_gain("I - K H", N_STATES, N_STATES, filtro->aux3);
        mulmat(N_STATES, N_STATES, N_STATES, filtro->aux3, filtro->aux4, filtro->cov); //cov = (I - Ki * H[i][:]) Pi
        //print_gain("cov", N_STATES, N_STATES, filtro->cov);
        //print_gain("H", N_STATES, N_STATES, filtro->cov);
    }
    /*
    transpose(N_OBS, N_STATES, filtro->H, filtro->Ht); 
    mulmat(N_STATES, N_STATES, N_OBS, filtro->cov, filtro->Ht, filtro->aux9); // aux9 = cov * Ht
    mulmat(N_OBS, N_STATES, N_OBS, filtro->H, filtro->aux9, filtro->aux10); // aux10 = H * cov * Ht
    accum(N_OBS, N_OBS, filtro->aux10, filtro->R); // aux10 = H * cov * Ht + R
    //print_gain("(H * cov * Ht + R): ",N_OBS, N_OBS,  filtro->aux10);
    cholsl(N_OBS, filtro->aux10, filtro->aux11, filtro->temp2); // aux11 = inv(H * cov * Ht + R)
    //print_gain("inv(H * cov * Ht + R): ",N_OBS, N_OBS,  filtro->aux11);
    mulmat(N_STATES, N_OBS, N_OBS, filtro->aux9, filtro->aux11, filtro->G); // G = cov * Ht * inv(H * cov * Ht + R) //TODO: EN VEZ DE 10 era 11
    //print_gain("G: ", N_STATES, N_OBS,  filtro->G);
    sub_vec(filtro->meas, filtro->exp_meas, filtro->temp2, N_OBS); // temp2 = z_medido - z_esperado
    mulvec(N_STATES, N_OBS, filtro->G, filtro->temp2, filtro->temp3); // aux8 = G(z_medido - z_esperado)
    //print_vector("Mediciones", N_OBS, filtro->meas);
    //print_vector("Mediciones esperadas", N_OBS, filtro->exp_meas);
    //print_vector("Innovacion", N_OBS, filtro->aux8);
    accum_vec(N_STATES, filtro->states, filtro->temp3); // estado = estado + G(z_medido - z_esperado)*/
    // Normalizo q
    filtro->q.q1 = filtro->states[N_P + N_V + 0];
    filtro->q.q2 = filtro->states[N_P + N_V + 1];
    filtro->q.q3 = filtro->states[N_P + N_V + 2];
    filtro->q.q4 = filtro->states[N_P + N_V + 3];
    quat_Normalization(&filtro->q);
    filtro->states[N_P + N_V] = filtro->q.q1;
    filtro->states[N_P + N_V + 1] = filtro->q.q2;
    filtro->states[N_P + N_V + 2] = filtro->q.q3;
    filtro->states[N_P + N_V + 3] = filtro->q.q4;
    /*
    mulmat(N_STATES, N_OBS, N_STATES, filtro->G, filtro->H, filtro->aux2); //  aux2 = G * H
    mulmat(N_STATES, N_STATES, N_STATES, filtro->aux2, filtro->cov, filtro->aux3); // aux3 = G * H * cov_priori
    mat_negate(N_STATES, N_STATES, filtro->aux3); // aux3 = - G * H * cov_priori
    accum(N_STATES, N_STATES, filtro->cov, filtro->aux3); // cov = cov_priori - G * H * cov_priori
    */
    if(filtro->states[2] < MIN_HEIGHT){
	filtro->states[2] = MIN_HEIGHT;
    }
    else if (filtro->states[2] > MAX_HEIGHT){
        filtro->states[2] = MAX_HEIGHT;
    }
}

double calc_trace_cov(ofs_ekf_t *filtro){
    double suma = 0;
    for (int i = 0; i < N_STATES; i++) {
        suma += filtro->cov[i][i];
    }
    return suma;
}