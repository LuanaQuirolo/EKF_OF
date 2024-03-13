#include "ekf.h"

void ofs_ekf_init(ofs_ekf_t* filtro, mediciones_t *meas){
    /* Inicializo el campo magnetico respecto a mundo */
    filtro->m[0] = (float) meas->mx;
    filtro->m[1] = (float) meas->my;
    filtro->m[2] = (float) meas->mz;
    filtro->qm = vec2quat(filtro->m);
    /* Inicializo el valor de la gravedad respecto a mundo */
    filtro->qg.q1 = (float) 0.0f;
    filtro->qg.q2 = (float) 0.0f;
    filtro->qg.q3 = (float) 0.0f;
    filtro->qg.q4 = (float) -g;
    /* Inicializo estados */
    vec_zeros(N_STATES, filtro->states); //p, v, q
    filtro->states[2] = (float) meas->range;
    float roll = (float) atan2f(meas->ay, meas->az);
    float pitch = (float) atan2f(-meas->ax, meas->ay * sinf(roll) + meas->az + cosf(roll));
    float yaw = (float) atan2f(-meas->my * cosf(roll) + meas->mz * sinf(roll), \
                meas->mx * cosf(pitch) + meas->my * sinf(roll) * sinf(pitch) + meas->mz * cosf(roll) * sinf(pitch));
    //yaw= (float)0.0f;//TODO BORRAR
    quat_euler(&filtro->q, &roll, &pitch, &yaw);
    filtro->states[N_P + N_V] = (float) filtro->q.q1;
    filtro->states[N_P + N_V + 1] = (float) filtro->q.q2;
    filtro->states[N_P + N_V + 2] = (float) filtro->q.q3;
    filtro->states[N_P + N_V + 3] = (float) filtro->q.q4;
    filtro->q_aux2 = quat_mult(filtro->qm, quat_conjugate(filtro->q)); // aux = qm * qk- (de cuerpo a mundo)
    filtro->q_aux = quat_mult(filtro->q, filtro->q_aux2); // aux = q * qm * q- (de cuerpo a mundo)
    filtro->qm.q1 = (float) filtro->q_aux.q1;
    filtro->qm.q2 = (float) filtro->q_aux.q2;
    filtro->qm.q3 = (float) filtro->q_aux.q3;
    filtro->qm.q4 = (float) filtro->q_aux.q4;
    filtro->m[0] = (float) filtro->qm.q2;
    filtro->m[1] = (float) filtro->qm.q3;
    filtro->m[2] = (float) filtro->qm.q4;
    /* Inicializo covarianza */
    mat_zeros(N_STATES, N_STATES, filtro->cov); // Matriz de covarianza de estados
    mat_addeye(N_STATES, filtro->cov);
    matmul_scalar(N_STATES, N_STATES, filtro->cov, 5.0f);
    /* Inicializo F */
    mat_zeros(N_STATES, N_STATES, filtro->F);
    mat_addeye(N_STATES, filtro->F);
    /* Inicializo W */
    float W[N_STATES][N_PROC_NOISE]; // Derivada de vector de estados respecto de ruidos
    mat_zeros(N_STATES, N_PROC_NOISE, W);
    /* Wk */
    // d(velocidad) / d(uax)
    W[N_P][0] = (float) 1.0f;
    // d(velocidad) / d(uay)
    W[N_P+1][1] = (float) 1.0f;
    // d(velocidad) / d(uaz)
    W[N_P+2][2] = (float) 1.0f;
    // d(q1) / d(uw)
    W[N_P+N_V][3] = (float) 1.0f;
    // d(q2) / d(uw)
    W[N_P+N_V+1][4] = (float) 1.0f;
    // d(q3) / d(uw)
    W[N_P+N_V+2][5] = (float) 1.0f;
    // d(q4) / d(uw)
    W[N_P+N_V+3][6] = (float) 1.0f;
    /* Inicializo Wt */
    float Wt[N_PROC_NOISE][N_STATES]; // Derivada de vector de estados respecto de ruidos
    transpose(N_STATES, N_PROC_NOISE, W, Wt);
    /* Inicializo Q */
    float Q[N_PROC_NOISE][N_PROC_NOISE]; // Ruido de proceso
    mat_zeros(N_PROC_NOISE, N_PROC_NOISE, Q);
    Q[0][0] = (float) U_P_AX; //uax
    Q[1][1] = (float) U_P_AY; //uay
    Q[2][2] = (float) U_P_AZ; //uaz
    Q[3][3] = (float) U_P_W; //uw
    Q[4][4] = (float) U_P_W; //uw
    Q[5][5] = (float) U_P_W; //uw
    Q[6][6] = (float) U_P_W; //uw
    /* Calculo W Q Wt dado que no varia segun la iteracion */
    mulmat(N_PROC_NOISE, N_PROC_NOISE, N_STATES, Q, Wt, filtro->m_aux); // aux4 = (float) Q * Wt
    mulmat(N_STATES, N_PROC_NOISE, N_STATES, W, filtro->m_aux, filtro->WQWt); //  W * Q * Wt
    /* Limpio H */
    mat_zeros(N_OBS, N_STATES, filtro->H);
    /* Defino R */
    mat_zeros(N_OBS, N_OBS, filtro->R);
    filtro->R[0][0] = (float) U_AX; //uax
    filtro->R[1][1] = (float) U_AY; //uay
    filtro->R[2][2] = (float) U_AZ; //uaz
    filtro->R[N_IMU][N_IMU] = (float) U_MAGX; //umagx
    filtro->R[N_IMU + 1][N_IMU + 1] = (float) U_MAGY; //umagy
    filtro->R[N_IMU + 2][N_IMU + 2] = (float) U_MAGZ; //umagz
    filtro->R[N_IMU + N_MAG][N_IMU + N_MAG] = (float) U_FLOWX; //u flowx
    filtro->R[N_IMU + N_MAG + 1][N_IMU + N_MAG + 1] = (float) U_FLOWY; //u flowy
    filtro->R[N_OBS - 1][N_OBS - 1] = (float) U_RANGE; // distance
    /* Calculo factor de conversion del flujo optico */
    filtro->f = (float) 35.0f / (2.0f * tanf(4.2f * (float)M_PI / (180.0f * 2.0f)));  // Factor de conversión
};

void prediction_step(ofs_ekf_t* filtro, mediciones_t *u){
    /* Asigno valores a las variables auxiliares */
    filtro->p[0] = (float) filtro->states[0];
    filtro->p[1] = (float) filtro->states[1];
    filtro->p[2] = (float) filtro->states[2];
    filtro->v[0] = (float) filtro->states[N_P+0];
    filtro->v[1] = (float) filtro->states[N_P+1];
    filtro->v[2] = (float) filtro->states[N_P+2];
    filtro->q.q1 = (float) filtro->states[N_P + N_V + 0];
    filtro->q.q2 = (float) filtro->states[N_P + N_V + 1];
    filtro->q.q3 = (float) filtro->states[N_P + N_V + 2];
    filtro->q.q4 = (float) filtro->states[N_P + N_V + 3];
    filtro->qa_meas.q1 = (float) 0.0f;
    filtro->qa_meas.q2 = (float) u->ax;
    filtro->qa_meas.q3 = (float) u->ay;
    filtro->qa_meas.q4 = (float) u->az;
    filtro->qw_meas.q1 = (float) 0.0f;
    filtro->qw_meas.q2 = (float) u->wx;
    filtro->qw_meas.q3 = (float) u->wy;
    filtro->qw_meas.q4 = (float) u->wz;

    /* Calculo Fk */
    // F tiene la diagonal en unos, se asigna al inicializar el filtro
    //d(posicion) / d(velocidad)
    filtro->F[0][N_P] = (float) u->dt;
    filtro->F[1][N_P + 1] = (float) u->dt;
    filtro->F[2][N_P + 2] = (float) u->dt;
    //d(velocidad) / d(quaternion)
    filtro->F[N_P][N_P+N_V] = (float) u->dt * (2.0f*(u->ax) * filtro->q.q1 - 2.0f*(u->ay)*filtro->q.q4 + 2.0f*(u->az)*filtro->q.q3);
    filtro->F[N_P+1][N_P+N_V] = (float) u->dt * (2.0f*(u->ax) * filtro->q.q4 + 2.0f*(u->ay)*filtro->q.q1 - 2.0f*(u->az)*filtro->q.q2);
    filtro->F[N_P+2][N_P+N_V] = (float) u->dt * (-2.0f*(u->ax) * filtro->q.q3 + 2.0f*(u->ay)*filtro->q.q2 + 2.0f*(u->az)*filtro->q.q1);
    filtro->F[N_P][N_P+N_V+1] = (float) u->dt * (2.0f*(u->ax) * filtro->q.q2 + 2.0f*(u->ay)*filtro->q.q3 + 2.0f*(u->az)*filtro->q.q4);
    filtro->F[N_P+1][N_P+N_V+1] = (float) u->dt * (2.0f*(u->ax) * filtro->q.q3 - 2.0f*(u->ay)*filtro->q.q2 - 2.0f*(u->az)*filtro->q.q1);
    filtro->F[N_P+2][N_P+N_V+1] = (float) u->dt * (2.0f*(u->ax) * filtro->q.q4 + 2.0f*(u->ay)*filtro->q.q1 - 2.0f*(u->az)*filtro->q.q2);
    filtro->F[N_P][N_P+N_V+2] = (float) u->dt * (-2.0f*(u->ax) * filtro->q.q3 + 2.0f*(u->ay)*filtro->q.q2 + 2.0f*(u->az)*filtro->q.q1);
    filtro->F[N_P+1][N_P+N_V+2] = (float) u->dt * (2.0f*(u->ax) * filtro->q.q2 + 2.0f*(u->ay)*filtro->q.q3 + 2.0f*(u->az)*filtro->q.q4);
    filtro->F[N_P+2][N_P+N_V+2] = (float) u->dt * (-2.0f*(u->ax) * filtro->q.q1 + 2.0f*(u->ay)*filtro->q.q4 - 2.0f*(u->az)*filtro->q.q3);
    filtro->F[N_P][N_P+N_V+3] = (float) u->dt * (-2.0f*(u->ax) * filtro->q.q4 - 2.0f*(u->ay)*filtro->q.q1 + 2.0f*(u->az)*filtro->q.q2);
    filtro->F[N_P+1][N_P+N_V+3] = (float) u->dt * (2.0f*(u->ax) * filtro->q.q1 - 2.0f*(u->ay)*filtro->q.q4 + 2.0f*(u->az)*filtro->q.q3);
    filtro->F[N_P+2][N_P+N_V+3] = (float) u->dt * (2.0f*(u->ax) * filtro->q.q2 + 2.0f*(u->ay)*filtro->q.q3 + 2.0f*(u->az)*filtro->q.q4);
    //d(q) / d(q1)
    //filtro->F[N_P+N_V][N_P+N_V] = (float) 1;
    filtro->F[N_P+N_V+1][N_P+N_V] = (float) u->wx * u->dt / 2.0f;
    filtro->F[N_P+N_V+2][N_P+N_V] = (float) u->wy * u->dt / 2.0f;
    filtro->F[N_P+N_V+3][N_P+N_V] = (float) u->wz * u->dt / 2.0f;
    //d(q) / d(q2)
    filtro->F[N_P+N_V][N_P+N_V+1] = (float) -u->wx * u->dt / 2.0f;
    //filtro->F[N_P+N_V+1][N_P+N_V+1] = (float) 1;
    filtro->F[N_P+N_V+2][N_P+N_V+1] = (float) -u->wz * u->dt / 2.0f;
    filtro->F[N_P+N_V+3][N_P+N_V+1] = (float) u->wy * u->dt / 2.0f;
    //d(q) / d(q3)
    filtro->F[N_P+N_V][N_P+N_V+2] = (float) -u->wy * u->dt / 2.0f;
    filtro->F[N_P+N_V+1][N_P+N_V+2] = (float) u->wz * u->dt / 2.0f;
    //filtro->F[N_P+N_V+2][N_P+N_V+2] = (float) 1;
    filtro->F[N_P+N_V+3][N_P+N_V+2] = (float) -u->wx * u->dt / 2.0f;
    //d(q) / d(q4)
    filtro->F[N_P+N_V][N_P+N_V+3] = (float) -u->wz * u->dt / 2.0f;
    filtro->F[N_P+N_V+1][N_P+N_V+3] = (float) -u->wy * u->dt / 2.0f;
    filtro->F[N_P+N_V+2][N_P+N_V+3] = (float) u->wx * u->dt / 2.0f;
    //filtro->F[N_P+N_V+3][N_P+N_V+3] = (float) 1;

    /* Calculo de Covarianza a priori */
    //Como W y Q no son variables, se calculan al inicializar el filtro por unica vez -> WQWt
    transpose(N_STATES, N_STATES, filtro->F, filtro->m_aux3);
    mulmat(N_STATES, N_STATES, N_STATES, filtro->cov, filtro->m_aux3, filtro->m_aux2); // aux5 = (float) cov * Ft
    mulmat(N_STATES, N_STATES, N_STATES, filtro->F, filtro->m_aux2, filtro->m_aux3); // aux7 = (float) F * cov * Ft
    add(N_STATES, N_STATES, filtro->WQWt, filtro->m_aux3, filtro->cov); // cov = (float) F * cov * Ft +  W * Q * Wt

    /* Calculo de estados */
    // Posicion
    vecmul_scalar2(N_V, filtro->v, filtro->v_aux, u->dt); // aux2 = (float) Vk = (float) dt * Vk
    add_vec(N_P, filtro->p, filtro->v_aux, filtro->states); // Pk+1 = (float) Pk + dt * Vk
    if(filtro->states[2] < (float)MIN_HEIGHT){
    	filtro->states[2] = (float) MIN_HEIGHT;
    }
    else if (filtro->states[2] > (float)MAX_HEIGHT){
    	filtro->states[2] = (float) MAX_HEIGHT;
    }
    // Velocidad
    // Queremos pasar la acel de cuerpo a mundo
    filtro->q_aux2 =  quat_mult(filtro->qa_meas, quat_conjugate(filtro->q)); // aux13 =  (filtro->qa_meas) * -q
    filtro->q_aux =  quat_mult(filtro->q, filtro->q_aux2); // aux13 = (float) q * (filtro->qa_meas) * -q
    quat_sub(&filtro->q_aux2, filtro->q_aux, filtro->qg); // aux13 = (float) q * (filtro->qa_meas) * -q - g
    quat2vec(filtro->q_aux2, filtro->v_aux2); //aux3 = (float) aux13
    vecmul_scalar(N_V, filtro->v_aux2, u->dt); // aux3 = (float) dt (a_meas - g)
    add_vec(N_V, filtro->v, filtro->v_aux2, filtro->states + N_P); // Vk+1 = (float) Vk + dt * (a - g)
    // Quaternion
    filtro->q_aux = quat_mult(filtro->q, filtro->qw_meas); // aux = (float) qk * qw
    quat_scalar(&filtro->q_aux, (float)u->dt / 2.0f); // aux = (float) (dt / 2) qk * qw
    quat_add(&filtro->q_aux2, filtro->q_aux, filtro->q); // aux = (float) qk + (dt / 2) qk * qw
    quat_Normalization(&filtro->q_aux2);
    filtro->states[N_P + N_V] = (float) filtro->q_aux2.q1;
    filtro->states[N_P + N_V + 1] = (float) filtro->q_aux2.q2;
    filtro->states[N_P + N_V + 2] = (float) filtro->q_aux2.q3;
    filtro->states[N_P + N_V + 3] = (float) filtro->q_aux2.q4;
};

void correction_step(ofs_ekf_t* filtro, mediciones_t* z){

    /* Asigno variables auxiliares */
    filtro->p[0] = (float) filtro->states[0];
    filtro->p[1] = (float) filtro->states[1];
    filtro->p[2] = (float) filtro->states[2];
    filtro->v[0] = (float) filtro->states[N_P+0];
    filtro->v[1] = (float) filtro->states[N_P+1];
    filtro->v[2] = (float) filtro->states[N_P+2];
    filtro->q.q1 = (float) filtro->states[N_P + N_V + 0];
    filtro->q.q2 = (float) filtro->states[N_P + N_V + 1];
    filtro->q.q3 = (float) filtro->states[N_P + N_V + 2];
    filtro->q.q4 = (float) filtro->states[N_P + N_V + 3];

    /* Calculo H */
    /************* IMU *************/
    // partial_a / partial_q1
    filtro->H[0][N_P + N_V] = (float) 2.0f * (float)g * (filtro->q.q3);
    filtro->H[1][N_P + N_V] = (float) - 2.0f * (float)g * (filtro->q.q2);
    filtro->H[2][N_P + N_V] = (float) - 2.0f * (float)g * (filtro->q.q1);
    // partial_a / partial_q2
    filtro->H[0][N_P + N_V + 1] = (float) -2.0f * (float)g * (filtro->q.q4);
    filtro->H[1][N_P + N_V + 1] = (float) -2.0f * (float)g * (filtro->q.q1);
    filtro->H[2][N_P + N_V + 1] = (float) 2.0f * (float)g * (filtro->q.q2);
    // partial_a / partial_q3
    filtro->H[0][N_P + N_V + 2] = (float) 2.0f * (float)g * (filtro->q.q1);
    filtro->H[1][N_P + N_V + 2] = (float) -2.0f * (float)g * (filtro->q.q4);
    filtro->H[2][N_P + N_V + 2] = (float) 2.0f * (float)g * (filtro->q.q3);
    // partial_a / partial_q3
    filtro->H[0][N_P + N_V + 3] = (float) -2.0f * (float)g * (filtro->q.q2);
    filtro->H[1][N_P + N_V + 3] = (float) -2.0f * (float)g * (filtro->q.q3);
    filtro->H[2][N_P + N_V + 3] = (float) -2.0f * (float)g * (filtro->q.q4);
    /************* MAG *************/
    // partial_mag / partial_q1
    filtro->H[N_IMU + 0][N_P + N_V] = (float) 2.0f * filtro->m[0] * filtro->q.q1 + 2.0f * filtro->m[1] * filtro->q.q4 - 2.0f * filtro->m[2] * filtro->q.q3;
    filtro->H[N_IMU + 1][N_P + N_V] = (float) -2.0f * filtro->m[0] * filtro->q.q4 + 2.0f * filtro->m[1] * filtro->q.q1 + 2.0f * filtro->m[2] * filtro->q.q2;
    filtro->H[N_IMU + 2][N_P + N_V] = (float) 2.0f * filtro->m[0] * filtro->q.q3 - 2.0f * filtro->m[1] * filtro->q.q2 + 2.0f * filtro->m[2] * filtro->q.q1;
    // partial_mag / partial_q2
    filtro->H[N_IMU + 0][N_P + N_V + 1] = (float) 2.0f * filtro->m[0] * filtro->q.q2 + 2.0f * filtro->m[1] * filtro->q.q3 + 2.0f * filtro->m[2] * filtro->q.q4;
    filtro->H[N_IMU + 1][N_P + N_V + 1] = (float) 2.0f * filtro->m[0] * filtro->q.q3 - 2.0f * filtro->m[1] * filtro->q.q2 + 2.0f * filtro->m[2] * filtro->q.q1;
    filtro->H[N_IMU + 2][N_P + N_V + 1] = (float) 2.0f * filtro->m[0] * filtro->q.q4 - 2.0f * filtro->m[1] * filtro->q.q1 - 2.0f * filtro->m[2] * filtro->q.q2;
    // partial_mag / partial_q3
    filtro->H[N_IMU + 0][N_P + N_V + 2] = (float) -2.0f * filtro->m[0] * filtro->q.q3 + 2.0f * filtro->m[1] * filtro->q.q2 - 2.0f * filtro->m[2] * filtro->q.q1;
    filtro->H[N_IMU + 1][N_P + N_V + 2] = (float) 2.0f * filtro->m[0] * filtro->q.q2 + 2.0f * filtro->m[1] * filtro->q.q3 + 2.0f * filtro->m[2] * filtro->q.q4;
    filtro->H[N_IMU + 2][N_P + N_V + 2] = (float) 2.0f * filtro->m[0] * filtro->q.q1 + 2.0f * filtro->m[1] * filtro->q.q4 - 2.0f * filtro->m[2] * filtro->q.q3;
    // partial_mag / partial_q4
    filtro->H[N_IMU + 0][N_P + N_V + 3] = (float) -2.0f * filtro->m[0] * filtro->q.q4 + 2.0f * filtro->m[1] * filtro->q.q1 + 2.0f * filtro->m[2] * filtro->q.q2;
    filtro->H[N_IMU + 1][N_P + N_V + 3] = (float) -2.0f * filtro->m[0] * filtro->q.q1 - 2.0f * filtro->m[1] * filtro->q.q4 + 2.0f * filtro->m[2] * filtro->q.q3;
    filtro->H[N_IMU + 2][N_P + N_V + 3] = (float) 2.0f * filtro->m[0] * filtro->q.q2 + 2.0f * filtro->m[1] * filtro->q.q3 + 2.0f * filtro->m[2] * filtro->q.q4;
    /************* FLUJO OPTICO *************/
    // partial_nx / partial_pz
    filtro->H[N_IMU + N_MAG][2] = (float) (z->dt * filtro->f * (2.0f * powf(filtro->q.q1, 2.0f) + 2.0f * powf(filtro->q.q4, 2.0f) - 1.0f)) * \
                            (filtro->q.q1 * (filtro->q.q1 * filtro->v[0] - filtro->q.q3 * filtro->v[2] + \
                            filtro->q.q4 * filtro->v[1]) - filtro->q.q2 * (-filtro->q.q2 * filtro->v[0] - \
                            filtro->q.q3 * filtro->v[1] - filtro->q.q4 * filtro->v[2]) \
                            -filtro->q.q3 * (filtro->q.q1 * filtro->v[2] - filtro->q.q2 * filtro->v[1] + \
                            filtro->q.q3 * filtro->v[0]) + filtro->q.q4 * (filtro->q.q1 * filtro->v[1] + \
                            filtro->q.q2 * filtro->v[2] - filtro->q.q4 * filtro->v[0]) \
                            )/ powf(filtro->p[2], 2.0f);
    // partial_nx / partial_v
    filtro->H[N_IMU + N_MAG][N_P + 0] = (float) - (z->dt * filtro->f * (2.0f * powf(filtro->q.q1, 2.0f) + 2.0f * powf(filtro->q.q4, 2.0f) - 1.0f)) * \
                            (powf(filtro->q.q1, 2.0f) + powf(filtro->q.q2, 2.0f) - powf(filtro->q.q3, 2.0f) - powf(filtro->q.q4, 2.0f)) / filtro->p[2];
    filtro->H[N_IMU + N_MAG][N_P + 1] = (float) - (z->dt * filtro->f * (2.0f *filtro->q.q1 *filtro->q.q4 + 2.0f * filtro->q.q2 * filtro->q.q3)) * \
                            (2.0f * powf(filtro->q.q1, 2.0f) + 2.0f * powf(filtro->q.q4, 2.0f) - 1.0f) / filtro->p[2];
    filtro->H[N_IMU + N_MAG][N_P + 2] = (float) - (z->dt * filtro->f * (-2.0f *filtro->q.q1 *filtro->q.q3 + 2.0f * filtro->q.q2 * filtro->q.q4)) * \
                            (2.0f * powf(filtro->q.q1, 2.0f) + 2.0f * powf(filtro->q.q4, 2.0f) - 1.0f) / filtro->p[2];
    // partial_nx / partial_q
    filtro->H[N_IMU + N_MAG][N_P + N_V + 0] = (float) z->dt * filtro->f * (-4* filtro->q.q1 * (filtro->q.q1 * (filtro->q.q1 * filtro->v[0] -\
                            filtro->q.q3 * filtro->v[2] +filtro->q.q4 * filtro->v[1]) \
                            +filtro->q.q2 * (filtro->q.q2 * filtro->v[0] +filtro->q.q3 * filtro->v[1] +filtro->q.q4 * filtro->v[2])\
                            -filtro->q.q3 * (filtro->q.q1 * filtro->v[2] -filtro->q.q2 * filtro->v[1] +filtro->q.q3 * filtro->v[0]) \
                            +filtro->q.q4 * (filtro->q.q1 * filtro->v[1] +filtro->q.q2 * filtro->v[2] -filtro->q.q4 * filtro->v[0])) \
                            -2.0f * (2.0f * powf(filtro->q.q1, 2.0f) + 2.0f * powf(filtro->q.q4, 2.0f) - 1.0f) \
                            * (filtro->q.q1 * filtro->v[0] -filtro->q.q3 * filtro->v[2] +filtro->q.q4 * filtro->v[1])) / filtro->p[2];
    filtro->H[N_IMU + N_MAG][N_P + N_V + 1] = (float) - 2.0f * z->dt * filtro->f * (2.0f * powf(filtro->q.q1, 2.0f) + 2.0f * powf(filtro->q.q4, 2.0f) - 1.0f) \
                            * (filtro->q.q2 * filtro->v[0] +filtro->q.q3 * filtro->v[1] +filtro->q.q4 * filtro->v[2]) / filtro->p[2];
    filtro->H[N_IMU + N_MAG][N_P + N_V + 2] = (float) 2.0f * z->dt * filtro->f * (2.0f * powf(filtro->q.q1, 2.0f) + 2.0f * powf(filtro->q.q4, 2.0f) - 1.0f) \
                            * (filtro->q.q1 * filtro->v[2] -filtro->q.q2 * filtro->v[1] +filtro->q.q3 * filtro->v[0]) / filtro->p[2];
    filtro->H[N_IMU + N_MAG][N_P + N_V + 3] = (float) z->dt * filtro->f * (-4.0f *filtro->q.q4 * (filtro->q.q1 * (filtro->q.q1 * filtro->v[0] \
                            - filtro->q.q3 * filtro->v[2] +filtro->q.q4 * filtro->v[1]) +filtro->q.q2 * (filtro->q.q2 * filtro->v[0] +\
                            filtro->q.q3 * filtro->v[1] +filtro->q.q4 * filtro->v[2]) -filtro->q.q3 * (filtro->q.q1 * filtro->v[2] \
                            - filtro->q.q2 * filtro->v[1] +filtro->q.q3 * filtro->v[0]) +filtro->q.q4 * (filtro->q.q1 * filtro->v[1] \
                            + filtro->q.q2 * filtro->v[2] -filtro->q.q4 * filtro->v[0])) -2.0f * (2.0f * powf(filtro->q.q1, 2) + 2.0f * powf(filtro->q.q4, 2) - 1.0f) \
                            * (filtro->q.q1 * filtro->v[1] +filtro->q.q2 * filtro->v[2] -filtro->q.q4 * filtro->v[0])) / filtro->p[2];
    // partial_ny / partial_pz
    filtro->H[N_IMU + N_MAG + 1][2] = (float) (z->dt * filtro->f * (2.0f * powf(filtro->q.q1, 2.0f) + 2.0f * powf(filtro->q.q4, 2.0f) - 1.0f)) * \
                            (filtro->q.q1 * (filtro->q.q1 * filtro->v[1] +filtro->q.q2 * filtro->v[2] -filtro->q.q4 * filtro->v[0]) \
                            +filtro->q.q2 * (filtro->q.q1 * filtro->v[2] -filtro->q.q2 * filtro->v[1] +filtro->q.q3 * filtro->v[0]) \
                            + filtro->q.q3 * (filtro->q.q2 * filtro->v[0] +filtro->q.q3 * filtro->v[1] +filtro->q.q4 * filtro->v[2]) \
                            -filtro->q.q4 * (filtro->q.q1 * filtro->v[0] -filtro->q.q3 * filtro->v[2] +filtro->q.q4 * filtro->v[1]) \
                            ) / powf(filtro->p[2], 2.0f);
    // partial_ny / partial_v
    filtro->H[N_IMU + N_MAG + 1][N_P + 0] = (float) (z->dt * filtro->f * (2.0f * filtro->q.q1 *filtro->q.q4 - 2.0f * filtro->q.q2 * filtro->q.q3)) * \
                            (2.0f * powf(filtro->q.q1, 2.0f) + 2.0f * powf(filtro->q.q4, 2.0f) - 1.0f) / filtro->p[2];
    filtro->H[N_IMU + N_MAG + 1][N_P + 1] = (float) - (z->dt * filtro->f * (2.0f * powf(filtro->q.q1, 2.0f) + 2.0f * powf(filtro->q.q4, 2.0f) - 1.0f)) * \
                            (powf(filtro->q.q1, 2.0f) - powf(filtro->q.q2, 2.0f) + powf(filtro->q.q3, 2.0f) - powf(filtro->q.q4, 2.0f)) / filtro->p[2];
    filtro->H[N_IMU + N_MAG + 1][N_P + 2] = (float) - (z->dt * filtro->f * (2.0f *filtro->q.q1 *filtro->q.q2 + 2.0f * filtro->q.q3 * filtro->q.q4)) * \
                            (2.0f * powf(filtro->q.q1, 2.0f) + 2.0f * powf(filtro->q.q4, 2.0f) - 1.0f) / filtro->p[2];
    // partial_ny / partial_q
    filtro->H[N_IMU + N_MAG + 1][N_P + N_V + 0] = (float) - 2.0f * z->dt * filtro->f * (2.0f *filtro->q.q1 * (filtro->q.q1 * (filtro->q.q1 * filtro->v[1] \
                            +filtro->q.q2 * filtro->v[2] -filtro->q.q4 * filtro->v[0]) \
                            +filtro->q.q2 * (filtro->q.q1 * filtro->v[2] -filtro->q.q2 * filtro->v[1] +filtro->q.q3 * filtro->v[0]) \
                            +filtro->q.q3 * (filtro->q.q2 * filtro->v[0] +filtro->q.q3 * filtro->v[1] +filtro->q.q4 * filtro->v[2]) \
                            -filtro->q.q4 * (filtro->q.q1 * filtro->v[0] -filtro->q.q3 * filtro->v[2] +filtro->q.q4 * filtro->v[1])) \
                            + (2.0f * powf(filtro->q.q1, 2.0f) + 2.0f * powf(filtro->q.q4, 2.0f) - 1.0f) \
                            * (filtro->q.q1 * filtro->v[1] +filtro->q.q2 * filtro->v[2] -filtro->q.q4 * filtro->v[0])) / filtro->p[2];
    filtro->H[N_IMU + N_MAG + 1][N_P + N_V + 1] = (float) - 2.0f * z->dt * filtro->f * (2.0f * powf(filtro->q.q1, 2.0f) + 2.0f * powf(filtro->q.q4, 2.0f) - 1.0f) \
                            * (filtro->q.q1 * filtro->v[2] -filtro->q.q2 * filtro->v[1] +filtro->q.q3 * filtro->v[0]) / filtro->p[2];
    filtro->H[N_IMU + N_MAG + 1][N_P + N_V + 2] = (float) - 2.0f * z->dt * filtro->f * (2.0f * powf(filtro->q.q1, 2.0f) + 2.0f * powf(filtro->q.q4, 2.0f) - 1.0f) \
                            * (filtro->q.q2 * filtro->v[0] +filtro->q.q3 * filtro->v[1] +filtro->q.q4 * filtro->v[2]) / filtro->p[2];
    filtro->H[N_IMU + N_MAG + 1][N_P + N_V + 3] = (float) - 2.0f * z->dt * filtro->f * (2.0f * filtro->q.q4 * (filtro->q.q1 * (filtro->q.q1 * filtro->v[1] \
                                +filtro->q.q2 * filtro->v[2] -filtro->q.q4 * filtro->v[0]) \
                            +filtro->q.q2 * (filtro->q.q1 * filtro->v[2] -filtro->q.q2 * filtro->v[1] +filtro->q.q3 * filtro->v[0]) \
                            +filtro->q.q3 * (filtro->q.q2 * filtro->v[0] +filtro->q.q3 * filtro->v[1] +filtro->q.q4 * filtro->v[2]) \
                            -filtro->q.q4 * (filtro->q.q1 * filtro->v[0] -filtro->q.q3 * filtro->v[2] +filtro->q.q4 * filtro->v[1])) \
                            - (2.0f * powf(filtro->q.q1, 2.0f) + 2.0f * powf(filtro->q.q4, 2.0f) - 1.0f) \
                            * (filtro->q.q1 * filtro->v[0] -filtro->q.q3 * filtro->v[2] +filtro->q.q4 * filtro->v[1])) / filtro->p[2];
    /************* DISTANCIA *************/
    // partial_tofs / partial_pz
    filtro->H[N_IMU + N_MAG + N_OFS][2] = (float) 1.0f/(2.0f * (powf(filtro->q.q1, 2.0f) + powf(filtro->q.q4, 2.0f)) - 1.0f);
    // partial_tofs / partial_q1
    filtro->H[N_IMU + N_MAG + N_OFS][N_P + N_V] = (float) -4.0f * filtro->p[2] * filtro->q.q1 * powf(filtro->H[N_IMU + N_MAG + N_OFS][2], 2.0f);
    // partial_tofs / partial_q4
    filtro->H[N_IMU + N_MAG + N_OFS][N_P + N_V + 3] = (float) -4.0f * filtro->p[2] * filtro->q.q4 * powf(filtro->H[N_IMU + N_MAG + N_OFS][2], 2.0f);

    /* Mediciones predichas */
    // Aceleracion
    filtro->q_aux2 = quat_mult(filtro->qg, filtro->q); // aux13 = qg * qk (de mundo a cuerpo)
    filtro->q_aux = quat_mult(quat_conjugate(filtro->q), filtro->q_aux2); // aux = (float) q- * qg * q (de mundo a cuerpo)
    quat2vec(filtro->q_aux, filtro->exp_meas); // measurements = (float) q- * qg * q
    filtro->meas[0] = (float) z->ax;
    filtro->meas[1] = (float) z->ay;
    filtro->meas[2] = (float) z->az;
    //Campo magnetico
    filtro->q_aux2 =  quat_mult(filtro->qm, filtro->q); // aux =  qm * qk (de mundo a cuerpo)
    filtro->q_aux =  quat_mult(quat_conjugate(filtro->q), filtro->q_aux2); // aux = (float) q- * qm * q (de mundo a cuerpo)
    quat2vec(filtro->q_aux, &filtro->exp_meas[N_IMU]); // measurements = (float) m_inicial
    filtro->meas[N_IMU + 0] = (float) z->mx;
    filtro->meas[N_IMU + 1] = (float) z->my;
    filtro->meas[N_IMU + 2] = (float) z->mz;
    // Flujo optico (Si corresponde)
    filtro->q_aux2 = vec2quat(&filtro->states[N_P]); //Paso la velocidad respecto al mundo a quat para convertirlo a body
    filtro->q_aux = quat_mult(filtro->q_aux2, filtro->q); // aux = v * q
    filtro->q_aux2 = quat_mult(quat_conjugate(filtro->q), filtro->q_aux); // aux13 = (float) q- * v * q
    filtro->exp_meas[N_IMU + N_MAG + 0] = (float) -z->dt * filtro->f * (filtro->q_aux2.q2 * (2.0f * (powf(filtro->q.q1, 2.0f) + powf(filtro->q.q4, 2.0f)) - 1.0f) / \
                                            filtro->p[2]  + z->wy);
    filtro->exp_meas[N_IMU + N_MAG + 1] = (float) -z->dt * filtro->f * (filtro->q_aux2.q3 * (2.0f * (powf(filtro->q.q1, 2.0f) + powf(filtro->q.q4, 2.0f)) - 1.0f) / \
                                            filtro->p[2]  - z->wx);
    filtro->meas[N_IMU + N_MAG + 0] = (float) z->ofx;
    filtro->meas[N_IMU + N_MAG + 1] = (float) z->ofy;
    // Sensor de distancia (Si corresponde)
    filtro->exp_meas[N_OBS - 1] = (float) filtro->p[2] / (2.0f * (powf(filtro->q.q1, 2.0f) + powf(filtro->q.q4, 2.0f)) - 1.0f);
    filtro->meas[N_OBS - 1] = (float) z->range;
    /*************************** Ganancia de Kalman **************************/ /*
    transpose(N_OBS, N_STATES, filtro->H, filtro->m_aux7);
    mulmat(N_STATES, N_STATES, N_OBS, filtro->cov, filtro->m_aux7, filtro->m_aux4); // aux9 = (float) cov * Ht
    mulmat(N_OBS, N_STATES, N_OBS, filtro->H, filtro->m_aux4, filtro->m_aux5); // aux10 = (float) H * cov * Ht
    accum(N_OBS, N_OBS, filtro->m_aux5, filtro->R); // aux10 = (float) H * cov * Ht + R
    //print_mat("(H * cov * Ht + R): ",N_OBS, N_OBS,  filtro->m_aux5);
    cholsl(N_OBS, filtro->m_aux5, filtro->m_aux6, filtro->v_aux4); // aux11 = (float) inv(H * cov * Ht + R)
    //print_mat("inv(H * cov * Ht + R): ",N_OBS, N_OBS,  filtro->m_aux6);
    mulmat(N_STATES, N_OBS, N_OBS, filtro->m_aux4, filtro->m_aux6, filtro->G); // G = (float) cov * Ht * inv(H * cov * Ht + R) //TODO: EN VEZ DE 10 era 11
    //print_mat("G: ", N_STATES, N_OBS,  filtro->G);
    */
    /*************************** Correccion **************************/ /*
    sub_vec(filtro->meas, filtro->exp_meas, filtro->v_aux4, N_OBS); // aux12 = (float) z_medido - z_esperado
    mulvec(N_STATES, N_OBS, filtro->G, filtro->v_aux4, filtro->v_aux3); // aux8 = (float) G(z_medido - z_esperado)
    //print_vector("Mediciones", N_OBS, filtro->meas);
    //print_vector("Mediciones esperadas", N_OBS, filtro->exp_meas);
    //print_vector("Innovacion", N_OBS, filtro->v_aux3);
    accum_vec(N_STATES, filtro->states, filtro->v_aux3); // estado = (float) estado + G(z_medido - z_esperado)*/

    /*************************** Covarianza **************************//*
    mulmat(N_STATES, N_OBS, N_STATES, filtro->G, filtro->H, filtro->m_aux4); //  aux9 = (float) G * H
    mulmat(N_STATES, N_STATES, N_STATES, filtro->m_aux4, filtro->cov, filtro->m_aux3); // aux7 = (float) G * H * cov_priori
    mat_negate(N_STATES, N_STATES, filtro->m_aux3); // aux7 = (float) - G * H * cov_priori
    accum(N_STATES, N_STATES, filtro->cov, filtro->m_aux3); // cov = (float) cov_priori - G * H * cov_priori
    */
    for (int i = 0; i < N_OBS; i++){
        mat_copy(N_STATES, N_STATES, filtro->cov, filtro->m_aux2); //m_aux2 = (float) cov = (float) Pi
        mat_getrow(N_OBS, N_STATES, filtro->H, i, filtro->v_aux3); // v_aux3 = (float) H[i][:]
        mulvec(N_STATES, N_STATES, filtro->m_aux2, filtro->v_aux3, filtro->v_aux5); // v_aux5 = (float) Pi * H[i][:].t
        vec_dot(N_STATES, filtro->v_aux3, filtro->v_aux5, &filtro->S); // S = (float) H[i][:] * Pi * H[i][:].t
        filtro->S += (float) filtro->R[i][i]; // S = (float) H[i][:] * Pi * H[i][:].t + Ri
        vecmul_scalar(N_STATES, filtro->v_aux5, 1.0f / filtro->S); //Ki = (float) Pi * H[i][:].t / (H[i][:] * Pi * H[i][:].t + Ri)
        // Correccion de estados
        vecmul_scalar2(N_STATES, filtro->v_aux5, filtro->v_aux6, filtro->meas[i] - filtro->exp_meas[i]); // v_aux6 = (float) Ki * (y[i] - y_est[i])
        accum_vec(N_STATES, filtro->states, filtro->v_aux6); // states = (float) states + Ki * (y[i] - y_est[i])
        // Correccion de covarianza
        mat_zeros(N_STATES, N_STATES, filtro->m_aux3); // m_aux3 = (float) 0
        mat_addeye(N_STATES, filtro->m_aux3); // m_aux3 = (float) I
        vec_outer(N_STATES, filtro->v_aux5, filtro->v_aux3, filtro->m_aux8); // m_aux8 = (float) Ki * H[i][:]
        mat_negate(N_STATES, N_STATES, filtro->m_aux8); // m_aux8 = (float) -Ki * H[i][:]
        accum(N_STATES, N_STATES, filtro->m_aux3, filtro->m_aux8); //m_aux3 = (float) I - Ki * H[i][:]
        //printf("Aca\n");
        //print_mat("cov", N_STATES, N_STATES, filtro->cov);
        //print_mat("I - K H", N_STATES, N_STATES, filtro->m_aux3);
        mulmat(N_STATES, N_STATES, N_STATES, filtro->m_aux3, filtro->m_aux2, filtro->cov); //cov = (float) (I - Ki * H[i][:]) Pi
        //print_mat("cov", N_STATES, N_STATES, filtro->cov);
        //print_mat("H", N_STATES, N_STATES, filtro->cov);
    }
    //print_mat("cov: ", N_STATES, N_STATES,  filtro->cov);
    filtro->q.q1 = (float) filtro->states[N_P + N_V + 0];
    filtro->q.q2 = (float) filtro->states[N_P + N_V + 1];
    filtro->q.q3 = (float) filtro->states[N_P + N_V + 2];
    filtro->q.q4 = (float) filtro->states[N_P + N_V + 3];
    quat_Normalization(&filtro->q);
    filtro->states[N_P + N_V] = (float) filtro->q.q1;
    filtro->states[N_P + N_V + 1] = (float) filtro->q.q2;
    filtro->states[N_P + N_V + 2] = (float) filtro->q.q3;
    filtro->states[N_P + N_V + 3] = (float) filtro->q.q4;
    if(filtro->states[2] < (float)MIN_HEIGHT){
        filtro->states[2] = (float) MIN_HEIGHT;
    }
    else if (filtro->states[2] > MAX_HEIGHT){
        filtro->states[2] = (float) MAX_HEIGHT;
    }
}

float calc_trace_cov(ofs_ekf_t *filtro) {
    float suma = 0.0f;
    for (int i = 0; i < N_STATES; i++) {
        suma += (float) filtro->cov[i][i];
    }
    return  (float)suma;
}
