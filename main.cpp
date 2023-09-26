#include "quaternions.hpp"
#include "matrix.hpp"
#include "EKF.hpp"
#include <stdio.h>
#include <iostream>
using namespace std;

int main(){
    /*
    float roll = 0.0, pitch = 0.0, yaw = 0.0;
    struct quaternion hola = {1, 0, 0, 0};
    struct quaternion chau = {1, 0, 0, 0};
    eulerAngles(hola, &roll, &pitch, &yaw);
    std::cout << quat_mult(hola, chau).q1 << quat_mult(hola, chau).q2<< quat_mult(hola, chau).q3<< quat_mult(hola, chau).q4 << std::endl;
    std::cout << hola.q1 << hola.q2<< hola.q3<< hola.q4 << std::endl;
    */
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
    std::cout << (int)(filtro.Npix) << std::endl;   */
    ofs_ekf_t filtro;
    ofs_ekf_init(&filtro);
    filtro.states[3] = 1;
    filtro.states[4] = 0.5;
    filtro.states[5] = 3;
    mediciones_t meas;
    prediction_step(&filtro, meas, 1);
    for (int i = 0; i < 3; i++){
        std::cout << filtro.states[i] << std::endl;
    }; 
    return 0;
}
