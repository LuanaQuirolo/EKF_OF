#include "quaternions.hpp"
#include "matrix.hpp"
#include "EKF.hpp"
#include <stdio.h>
#include <iostream>
using namespace std;

void print_states(ofs_ekf_t filtro); 
void print_cov(ofs_ekf_t *filtro);
int main(){
   /*
    int size = 3;
    double P[size][size];
    double V[size];
    double R[size];
    zeros(*P, size, size);
    zeros(V, size, 1);
    zeros(R, size, 1);
    P[0][0] = 3;
    P[1][1] = 3;
    P[2][2] = 3;
    ofs_ekf_t filtro;
    ofs_ekf_init(&filtro);
    filtro.Npix = 35;
    cholsl(*P, *P, V, size);
    mulvec(*P, V, R, size, size);
    matmul_scalar(V, 3, 1, 0.5/0.866025);
    for (int i = 0; i < size; i++){
        for (int j = 0; j < size; j++){
            std::cout << ((filtro.cov)[i][j]) << " ";
        };
        std::cout << std::endl;
    };  
    for (int i = 0; i < size; i++){
        std::cout << V[i] << std::endl;
    }; 
    std::cout << (filtro.g).q1 << std::endl;
    std::cout << (filtro.g).q4 << std::endl;
    std::cout << (int)(filtro.Npix) << std::endl;
    for (int i = 0; i < N_STATES; i++){
        for (int j = 0; j < N_STATES; j++){
            std::cout << ((*filtro).cov[i][j]) << " ";
        };
        std::cout << std::endl;
    };     */
    ofs_ekf_t filtro;
    ofs_ekf_init(&filtro);
    filtro.states[3] = 0;
    filtro.states[4] = 0;
    filtro.states[5] = 0;
    mediciones_t meas = {0.1, 0, 0, 9.81, 0, 0, 0.1, 0, 0, 1};
    for (int i = 0; i<2; i++){
        prediction_step(&filtro, meas);
        print_states(filtro);
        //print_cov(&filtro);
    }
    return 0;
}


void print_states(ofs_ekf_t filtro){
    float roll = 0.0, pitch = 0.0, yaw = 0.0;
        std::cout << "Posicion" << std::endl;
    for (int i = 0; i < N_P; i++){
        std::cout << filtro.states[i] << " ";
    };
    std::cout << std::endl;
    std::cout << "Velocidad" << std::endl;
    for (int i = 0; i < N_V; i++){
        std::cout << filtro.states[N_P + i] << " ";
    }; 
    std::cout << std::endl;
    std::cout << "Quaterniones" << std::endl;
    for (int i = 0; i < N_Q; i++){
        std::cout << filtro.states[N_P + N_V + i] << " ";
    }; 
    std::cout << std::endl;
    std::cout << "Angulos Euler" << std::endl;
    quaternion_t aux = {filtro.states[N_P + N_V], filtro.states[N_P + N_V +1], filtro.states[N_P + N_V+2], filtro.states[N_P + N_V+3]};
    eulerAngles(aux, &roll, &pitch, &yaw);
    std::cout << roll << " " << pitch << " " << yaw << std::endl;
}

void print_cov(ofs_ekf_t *filtro){
    std::cout << "Covarianza" << std::endl;
    for (int i = 0; i < N_STATES; i++){
        for (int j = 0; j < N_STATES; j++){
            std::cout << ((*filtro).cov[i][j]) << " ";
        };
        std::cout << std::endl;
    };
}