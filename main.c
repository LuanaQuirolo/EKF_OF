#include "quaternions.h"
#include "matrix.h"
#include "EKF.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

void print_states(ofs_ekf_t filtro);
void print_cov(ofs_ekf_t *filtro);
void print_expected_measurements(ofs_ekf_t filtro);

void print_expected_measurements(ofs_ekf_t filtro){
    printf("Mediciones esperadas\n");
    for (int i = 0; i < N_OBS_11; i++) {
        printf("%f ", filtro.exp_meas[i]);
    }
    printf("\n");
}
void print_measurements(mediciones_t meas){
    printf("Mediciones\n");
    printf("ax: %f ", meas.ax);
    printf("ay: %f ", meas.ay);
    printf("az: %f ", meas.az);
    printf("wx: %f ", meas.wx);
    printf("wy: %f ", meas.wy);
    printf("wz: %f ", meas.wz);
    printf("ofx: %f ", meas.ofx);
    printf("ofy: %f ", meas.ofy);
    printf("z: %f ", meas.range);

    printf("\n");
}

void print_states(ofs_ekf_t filtro) {
    float roll = 0.0, pitch = 0.0, yaw = 0.0;
    printf("Posicion\n");
    for (int i = 0; i < N_P; i++) {
        printf("%f ", filtro.states[i]);
    }
    printf("\n");
    printf("Velocidad\n");
    for (int i = 0; i < N_V; i++) {
        printf("%f ", filtro.states[N_P + i]);
    }
    printf("\n");
    printf("Quaterniones\n");
    for (int i = 0; i < N_Q; i++) {
        printf("%f ", filtro.states[N_P + N_V + i]);
    }
    printf("\n");
    printf("Angulos Euler\n");
    quaternion_t aux = {filtro.states[N_P + N_V], filtro.states[N_P + N_V + 1],
                       filtro.states[N_P + N_V + 2], filtro.states[N_P + N_V + 3]};
    eulerAngles(aux, &roll, &pitch, &yaw);
    printf("%f %f %f\n", roll, pitch, yaw);
}

void print_cov(ofs_ekf_t *filtro) {
    printf("Covarianza\n");
    for (int i = 0; i < N_STATES; i++) {
        for (int j = 0; j < N_STATES; j++) {
            printf("%f ", filtro->cov[i][j]);
        }
        printf("\n");
    }
}

void print_trace_cov(ofs_ekf_t *filtro) {
    printf("Traza de covarianza\n");
    double suma = 0;
    for (int i = 0; i < N_STATES; i++) {
        suma += filtro->cov[i][i];
    }
    printf("%f \n", suma);
}

void print_jac(double* jacobiano) {
    printf("Jacobiano\n");
    for (int i = 0; i < N_CORR_NOISE; i++) {
        for (int j = 0; j < N_STATES; j++) {
            printf("%f ", jacobiano[i*N_STATES+j]);
        }
        printf("\n");
    }
}

int main() {
    ofs_ekf_t filtro;
    ofs_ekf_init(&filtro);
    filtro.states[2] = 1;
    filtro.states[3] = 0;
    filtro.states[4] = 0;
    filtro.states[5] = 0;
    //print_states(filtro);
    srand(time(NULL));
    filtro.beta = 1;
    filtro.gamma = 1;
    mediciones_t meas = {0.01, 0.01, 0, 0, -9.81, 0, 0, 0.1, 0, 0, 1}; //dt, tau, ax, ay, az, wx, wy, wz, ofx, ofy, range
    for (int i = 0; i < 1; i++){
        /* RUIDOS */ 
        meas.ax = 0 + 0.0001 * (2.0 * rand() / (RAND_MAX) - 1.0);
        meas.ay = 0 + 0.0001 * (2.0 * rand() / (RAND_MAX) - 1.0);
        meas.az = -9.81 + 0.0001 * (2.0 * rand() / (RAND_MAX) - 1.0);
        meas.wx = 0 + 0.00001 * (2.0 * rand() / (RAND_MAX) - 1.0);
        meas.wy = 0 + 0.00001 * (2.0 * rand() / (RAND_MAX) - 1.0);
        meas.wz = 0.1 + 0.00001 * (2.0 * rand() / (RAND_MAX) - 1.0);
        meas.ofx = 0 + 0.5 * (2.0 * rand() / (RAND_MAX) - 1.0);
        meas.ofy =  0 + 0.5 * (2.0 * rand() / (RAND_MAX) - 1.0);
        meas.range = 1 + 0.001 * (2.0 * rand() / (RAND_MAX) - 1.0);
        filtro.beta = 1;
        filtro.gamma = 1;
        //print_measurements(meas);
        prediction_step(&filtro, meas);
        //print_states(filtro);
        print_trace_cov(&filtro);
        correction_step(&filtro, &meas);
        print_states(filtro);
        //print_expected_measurements(filtro);
        //print_trace_cov(&filtro);
        //printf("----------------------------------------\n");
    }
    quaternion_t aux = {0.35, -0.019, 0.0009, 0.936};
    float roll, pitch, yaw;
    eulerAngles(aux, &roll, &pitch, &yaw);
    printf("%f, %f, %f \n", roll, pitch, yaw);
    //print_states(filtro);
    /*int size = 2;
    double matrix[size][size];
    double inv[size][size];
    double aux[size];
    mat_zeros(*matrix, size, size);
    mat_zeros(*inv, size, size);
    matrix[0][0] = 10;
    matrix[1][1] = 10;
    matrix[0][1] = 10;
    cholsl(*matrix, *inv, aux, size); // aux11 = inv(H * cov * Ht + R)
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            printf("%f ", matrix[i][j]);
        }
        printf("\n");
    }
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            printf("%f ", inv[i][j]);
        }
        printf("\n");
    }*/
    return 0;
}
