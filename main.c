#include "quaternions.h"
#include "matrix.h"
#include "ekf.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <stdlib.h>

#define MAX_LINE_LENGTH 256
#define NUM_VALUES 14

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
    double suma = 0;
    for (int i = 0; i < N_STATES; i++) {
        suma += filtro.cov[i][i];
    }
    printf("trace: %f \n", suma);
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
    float roll, pitch, yaw;
    quaternion_t aux;
    char buffer_Tx_PC[256] = {0};
    char line[MAX_LINE_LENGTH];
    double values[NUM_VALUES];

    FILE * file_read = fopen("/Users/luana/Documents/Mediciones/rectangular.txt", "r"); // Reemplaza "datos.txt" por el nombre de tu archivo
    FILE * file_write = fopen("/Users/luana/Documents/Mediciones/rectangular_procesado.txt", "w");
    if (file_read == NULL) {
        perror("Error al abrir el archivo");
        return 1;
    }
    if (file_write == NULL) {
        perror("Error al abrir el archivo");
        return 1;
    }

    //Leo las primeras mediciones
    fgets(line, sizeof(line), file_read);
    char *token;
    int index = 0;
    token = strtok(line, ", ");
    while (token != NULL && index < NUM_VALUES) {
        values[index++] = atof(token);
        token = strtok(NULL, ", ");
    }
    ofs_ekf_t filtro;
    //double temp[3] = {-189/1090.0, 313/1090.0, 12/1090.0};
    values[10] = - values[10];
    values[7] = ((values[7] * 1090 / 1.2  + 726) + 1120) * 1.03 / 1090;
    values[8] = ((values[8] * 1090 / 1.48 - 588) + 1401) * 0.87/ 1090;
    values[9] = ((values[9] * 1090 / 0.67 - 1302) - 1078) * 1.13 / 1090;
    mediciones_t meas = {values[0], values[1], values[2], values[3], values[4], values[5], values[6], values[7], values[8], values[9], (double)values[10], (double)values[11], values[13]}; //dt, ax, ay, az, wx, wy, wz, mx, my, mz, ofx, ofy, range
    ofs_ekf_init(&filtro, &meas);
    print_states(filtro);
    srand(time(NULL));
    while (fgets(line, sizeof(line), file_read)){
        char *token;
        index = 0;

        token = strtok(line, ", ");
        while (token != NULL && index < NUM_VALUES) {
            values[index++] = atof(token);
            token = strtok(NULL, ", ");
        }
        /* RUIDOS */ 
        meas.dt = values[0];
        meas.ax = values[1];
        meas.ay = values[2];
        meas.az = values[3];
        meas.wx = values[4];
        meas.wy = values[5];
        meas.wz = values[6];
        values[7] = ((values[7] * 1090 / 1.2  + 726) + 1120) * 1.03 / 1090;
        values[8] = ((values[8] * 1090 / 1.48 - 588) + 1401) * 0.87/ 1090;
        values[9] = ((values[9] * 1090 / 0.67 - 1302) - 1078) * 1.13 / 1090;
        meas.mx = values[7];
        meas.my = values[8];
        meas.mz = values[9];
        if(values[12] != 0){
            meas.ofx = -(double)values[10];
            meas.ofy =  (double)values[11];
        }
        else{
            meas.ofx = (double)0;
            meas.ofy =  (double)0;
        }
        meas.range = values[13];
        //print_measurements(meas);
        prediction_step(&filtro, &meas);
        //print_states(filtro);
        //print_trace_cov(&filtro);
        correction_step(&filtro, &meas);
        aux.q1 = filtro.states[N_P + N_V];
        aux.q2 = filtro.states[N_P + N_V + 1];
        aux.q3 = filtro.states[N_P + N_V + 2];
        aux.q4 = filtro.states[N_P + N_V + 3];
        eulerAngles(aux, &roll, &pitch, &yaw);
        fprintf(file_write, "%lf %s %lf %s %lf %s %lf %s %lf %s %lf %s %lf %s %lf %s %lf %s %lf %s %lf \n", meas.dt, ", ", filtro.states[0], ", ", filtro.states[1], ", ", filtro.states[2], ", ",  //Posicion x, y, z
			  filtro.states[3], ", ", filtro.states[4], ", ", filtro.states[5], ", ",  // Velocidad vx, vy, vz
			  roll, ", ", pitch, ", ", yaw, ", ", calc_trace_cov(&filtro)); 
        print_states(filtro);
        //print_expected_measurements(filtro);
        //print_trace_cov(&filtro);
        //printf("----------------------------------------\n");
    }
    fclose(file_read);
    fclose(file_write);
    return 0;
}
