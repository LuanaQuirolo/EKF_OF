#include "quaternions.h"
#include "matrix.h"
#include "ekf.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

void print_states(ofs_ekf_t filtro);
void print_cov(ofs_ekf_t *filtro);
void print_expected_measurements(ofs_ekf_t filtro);

void print_expected_measurements(ofs_ekf_t filtro){
    printf("Mediciones esperadas\n");
    for (int i = 0; i < N_OBS; i++) {
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
    for (int i = 0; i < N_OBS; i++) {
        for (int j = 0; j < N_STATES; j++) {
            printf("%f ", jacobiano[i*N_STATES+j]);
        }
        printf("\n");
    }
}

int main() {
    ofs_ekf_t filtro;
    double temp[3] = {-189/1090, 313/1090, 12/1090};
    ofs_ekf_init(&filtro, temp, temp, 0.8);
    filtro.states[3] = 1;
    filtro.states[4] = 1;
    filtro.states[5] = 1;
    //print_states(filtro);
    srand(time(NULL));
    mediciones_t meas = {0.01, 0, 0, -9.81, 0, 0, 0, -189/1090, 313/1090, 12/1090, 0, 0, 0.8}; //dt, tau, ax, ay, az, wx, wy, wz, ofx, ofy, range
    for (int i = 0; i < 1000; i++){
        /* RUIDOS */ 
        meas.ax = 0  + 0.008 * (2.0 * rand() / (RAND_MAX) - 1.0);
        meas.ay = 0  + 0.008 * (2.0 * rand() / (RAND_MAX) - 1.0);
        meas.az = -9.81  + 0.008 * (2.0 * rand() / (RAND_MAX) - 1.0);
        meas.mx = -189/1090  + 0.0001 * (2.0 * rand() / (RAND_MAX) - 1.0);
        meas.my = 313/1090  + 0.0001 * (2.0 * rand() / (RAND_MAX) - 1.0);
        meas.mz = 12/1090  + 0.0001 * (2.0 * rand() / (RAND_MAX) - 1.0);
        meas.wx = 0  + 0.001 * (2.0 * rand() / (RAND_MAX) - 1.0);
        meas.wy = 0  + 0.001 * (2.0 * rand() / (RAND_MAX) - 1.0);
        meas.wz = 0  + 0.001 * (2.0 * rand() / (RAND_MAX) - 1.0);
        meas.ofx = 1 * meas.dt * filtro.f / meas.range + 0.2 * (2.0 * rand() / (RAND_MAX) -0.0);
        meas.ofy =  1 * meas.dt * filtro.f / meas.range + 0.2 * (2.0 * rand() / (RAND_MAX) -01.0);
        meas.range = 0.8 + 1 * meas.dt * i  + 0.00001 * (2.0 * rand() / (RAND_MAX) - 1.0);
        //print_measurements(meas);
        prediction_step(&filtro, &meas);
        //print_states(filtro);
        //print_trace_cov(&filtro);
        correction_step(&filtro, &meas);
        print_states(filtro);
        //print_expected_measurements(filtro);
        print_trace_cov(&filtro);
        //printf("----------------------------------------\n");
    }
    return 0;
}
