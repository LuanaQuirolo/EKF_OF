#include "quaternions.h"
#include "matrix.h"
#include "ekf.h"
#include <stdio.h>
#include <string.h>
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
    float roll = 0.0f, pitch = 0.0f, yaw = 0.0f;
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
    float suma = 0.0f;
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
    float suma = 0;
    for (int i = 0; i < N_STATES; i++) {
        suma += filtro->cov[i][i];
    }
    printf("%f \n", suma);
}

void print_jac(float* jacobiano) {
    printf("Jacobiano\n");
    for (int i = 0; i < N_OBS; i++) {
        for (int j = 0; j < N_STATES; j++) {
            printf("%f ", jacobiano[i*N_STATES+j]);
        }
        printf("\n");
    }
}

int main(int argc, char *argv[]) {
    float roll, pitch, yaw;
    quaternion_t aux;
    char buffer_Tx_PC[256] = {0};
    char line[MAX_LINE_LENGTH];
    float values[NUM_VALUES];
    char path_read[256] = "/Users/luana/Documents/Mediciones4/";
    char path_write[256] = "/Users/luana/Documents/Mediciones4/";
    char path_cov[256] = "/Users/luana/Documents/Mediciones4/";
    strcat(path_read, argv[1]);
    strcat(path_write, argv[1]);
    strcat(path_cov, argv[1]);
    strcat(path_read, ".txt");
    strcat(path_write, "_procesado.txt");
    strcat(path_cov, "_cov.txt");
    FILE * file_read = fopen(path_read, "r"); // Reemplaza "datos.txt" por el nombre de tu archivo
    FILE * file_write = fopen(path_write, "w");
    FILE * file_cov = fopen(path_cov, "w");
    if (file_read == NULL) {
        perror("Error al abrir el archivo");
        return 1;
    }
    if (file_write == NULL) {
        perror("Error al abrir el archivo");
        return 1;
    }
    if (file_cov == NULL) {
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
    //float temp[3] = {-189/1090.0, 313/1090.0, 12/1090.0};
    values[7] = (values[7] - 705.00) * 1.21f /1090.0f;//+ 178.0f) * 1.15f /1090.0f;// + 707) * 1.11 / 1090.0;
    values[8] = (values[8] - 1228.00) * 0.83f /1090.0f;//- 1153.0f) * 0.9f /1090.0f;// - 1200) * 0.96 / 1090.0;
    values[9] = (values[9] + 987.00) * 1.03f /1090.0f;//- 78.0f) * 0.98f /1090.0f;// + 117) * 0.95 / 1090.0;
    mediciones_t meas = {values[0], values[1], values[2], values[3], values[4], values[5], values[6], values[7], values[8], values[9], (float)values[10], (float)values[11], values[13]}; //dt, ax, ay, az, wx, wy, wz, mx, my, mz, ofx, ofy, range
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
        meas.mx = (values[7] - 705.00) * 1.21f /1090.0f;//+ 178.0f) * 1.15f /1090.0f;// + 707) * 1.11 / 1090.0;
        meas.my = (values[8] - 1228.00) * 0.83f /1090.0f;//- 1153.0f) * 0.9f /1090.0f;// - 1200) * 0.96 / 1090.0;
        meas.mz = (values[9] + 987.00) * 1.03f /1090.0f;//- 78.0f) * 0.98f /1090.0f;// + 117) * 0.95 / 1090.0;
        if(fabs(values[12]) > 25){
            meas.ofx = (float)values[10];
            meas.ofy = (float)values[11];
        }
        else{
            meas.ofx = 0.0f;
            meas.ofy =  0.0f;
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
        fprintf(file_write, "%.6f %s %.6f %s %.6f %s %.6f %s %.6f %s %.6f %s %.6f %s %.6f %s %.6f %s %.6f %s %.6f \n", meas.dt, ", ", filtro.states[0], ", ", filtro.states[1], ", ", filtro.states[2], ", ",  //Posicion x, y, z
			  filtro.states[3], ", ", filtro.states[4], ", ", filtro.states[5], ", ",  // Velocidad vx, vy, vz
			  roll, ", ", pitch, ", ", yaw, ", ", calc_trace_cov(&filtro)); 
        fprintf(file_cov, "%.6f %s %.6f %s %.6f %s %.6f %s %.6f %s %.6f %s %.6f %s %.6f %s %.6f %s %.6f \n", filtro.cov[0][0], ", ", filtro.cov[1][1], ", ", filtro.cov[2][2], ", ",
                filtro.cov[3][3], ", ", filtro.cov[4][4], ", ", filtro.cov[5][5], ", ",
                filtro.cov[6][6], ", ", filtro.cov[7][7], ", ", filtro.cov[8][8], ", ", filtro.cov[9][9]); 
        //print_states(filtro);
        //print_expected_measurements(filtro);
        //print_trace_cov(&filtro);
        //printf("----------------------------------------\n");
    }
    fclose(file_read);
    fclose(file_write);
    fclose(file_cov);
    return 0;
}
