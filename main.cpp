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
    int size = 3;
    double P[size][size];
    //double V[size];
    double R[size];
    zeros(*P, size, size);
    //zeros(V, size, 1);
    zeros(R, size, 1);
    P[0][0] = 3;
    P[1][1] = 3;
    P[2][2] = 3;
    double V[] = { 1.0, 2.5, 3.7 }; 
    //cholsl(*P, *P, V, size);
    //mulvec(*P, V, R, size, size);
    for (int i = 0; i < size; i++){
        for (int j = 0; j < size; j++){
            std::cout << (P[i][j]) << " ";
        };
        std::cout << std::endl;
    };  
    for (int i = 0; i < size; i++){
        std::cout << V[i] << std::endl;
    };    
    return 0;
}
