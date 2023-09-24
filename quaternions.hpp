/***************************************************************
 *                  Cuaterniones
 *                  Creado: 10 sep. 2023
 *                  Autor: Luana Quirolo
 *                  Padron: 102102
****************************************************************/

#ifndef QUATERNIONS_H_
#define QUATERNIONS_H_

#ifndef PI  
    #define PI 3.14159265358979f
#endif

//#include <math.h>
//#include <stdio.h>
#include <cmath>

typedef struct quaternion{
    float q1;
    float q2;
    float q3;
    float q4;
} quaternion_t;


// Multiply two quaternions and return a copy of the result, prod = L * R
quaternion_t quat_mult(quaternion_t L, quaternion_t R);

// Multiply a reference of a quaternion by a scalar, q = s*q
static inline void quat_scalar(quaternion_t * q, float scalar);

// Adds two quaternions together and the sum is the pointer to another quaternion, Sum = L + R
static inline void quat_add(quaternion_t * Sum, quaternion_t L, quaternion_t R);
// Subtracts two quaternions together and the sum is the pointer to another quaternion, sum = L - R
static inline void quat_sub(quaternion_t * Sum, quaternion_t L, quaternion_t R);

// the conjugate of a quaternion is it's imaginary component sign changed  q* = [s, -v] if q = [s, v]
static inline quaternion_t quat_conjugate(quaternion_t q);

// norm of a quaternion is the same as a complex number
// sqrt( q1^2 + q2^2 + q3^2 + q4^2)
// the norm is also the sqrt(q * conjugate(q)), but thats a lot of operations in the quaternion multiplication
static inline float quat_Norm (quaternion_t q);

// Normalizes pointer q by calling quat_Norm(q),
static inline void quat_Normalization(quaternion_t * q);

/*
 returns as pointers, roll pitch and yaw from the quaternion
 Assume right hand system
 Roll is about the x axis, represented as phi
 Pitch is about the y axis, represented as theta
 Yaw is about the z axis, represented as psi 
 */
void eulerAngles(quaternion_t q, float* roll, float* pitch, float* yaw);

#endif /* QUATERNIONS_H_ */