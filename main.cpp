#include "quaternions.hpp"
#include <stdio.h>
#include <iostream>
using namespace std;

int main(){
    float roll = 0.0, pitch = 0.0, yaw = 0.0;
    struct quaternion hola = {1, 0, 0, 0};
    struct quaternion chau = {1, 0, 0, 0};
    eulerAngles(hola, &roll, &pitch, &yaw);
    std::cout << quat_mult(hola, chau).q1 << quat_mult(hola, chau).q2<< quat_mult(hola, chau).q3<< quat_mult(hola, chau).q4 << std::endl;
    std::cout << hola.q1 << hola.q2<< hola.q3<< hola.q4 << std::endl;
}
