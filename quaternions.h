/***************************************************************
 *                  Cuaterniones
 *                  Creado: 10 sep. 2023
 *                  Autor: Luana Quirolo
 *                  Padron: 102102
****************************************************************/

#ifndef QUATERNIONS_H_
#define QUATERNIONS_H_

#include <math.h>

typedef struct quaternion{
    float q1;
    float q2;
    float q3;
    float q4;
} quaternion_t;


// Multiply two quaternions and return a copy of the result, prod = L * R
quaternion_t quat_mult (struct quaternion L, struct quaternion R);

// Multiply a reference of a quaternion by a scalar, q = s*q
void quat_scalar(struct quaternion * q, float scalar);

// Adds two quaternions together and the sum is the pointer to another quaternion, Sum = L + R
void quat_add(struct quaternion * Sum, struct quaternion L, struct quaternion R);

// Subtracts two quaternions together and the sum is the pointer to another quaternion, sum = L - R
void quat_sub(struct quaternion * Sum, struct quaternion L, struct quaternion R);

// the conjugate of a quaternion is it's imaginary component sign changed  q* = [s, -v] if q = [s, v]
quaternion_t quat_conjugate(struct quaternion q);

// norm of a quaternion is the same as a complex number
// sqrt( q1^2 + q2^2 + q3^2 + q4^2)
// the norm is also the sqrt(q * conjugate(q)), but thats a lot of operations in the quaternion multiplication
float quat_Norm (quaternion_t q);

// Normalizes pointer q by calling quat_Norm(q),
void quat_Normalization(struct quaternion * q);

// Extends a vector 3D to a quaternion
quaternion_t vec2quat(float* vector);

// Returns vector from quaternion
void quat2vec(quaternion_t q, float* vector);

/*
 returns as pointers, roll pitch and yaw from the quaternion
 Assume right hand system
 Roll is about the x axis, represented as phi
 Pitch is about the y axis, represented as theta
 Yaw is about the z axis, represented as psi 
 */
void eulerAngles(struct quaternion q, float* roll, float* pitch, float* yaw);

/*
Secuencia XYZ. Se asumen angulos en radianes
 returns as pointers, roll pitch and yaw from the quaternion
 Assume right hand system
 Roll is about the x axis, represented as phi
 Pitch is about the y axis, represented as theta
 Yaw is about the z axis, represented as psi 
 */
void quat_euler(struct quaternion *q, float* roll, float* pitch, float* yaw);

#endif /* QUATERNIONS_H_ */