#include "quaternions.h"
#include "matrix.h"
#include "EKF.h"
#include <stdio.h>

void print_states(ofs_ekf_t filtro);
void print_cov(ofs_ekf_t *filtro);

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

int main() {
    ofs_ekf_t filtro;
    ofs_ekf_init(&filtro);
    filtro.states[3] = 0;
    filtro.states[4] = 0;
    filtro.states[5] = 0;
    mediciones_t meas = {0.1, 0, 0, 9.81, 0, 0, 0.1, 0, 0, 1};
    for (int i = 0; i < 2; i++) {
        prediction_step(&filtro, meas);
        print_states(filtro);
        // print_cov(&filtro);
    }
    return 0;
}
