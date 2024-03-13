/***************************************************************
 *                  Cuaterniones
 *                  Creado: 10 sep. 2023
 *                  Autor: Luana Quirolo
 *                  Padron: 102102
****************************************************************/

#include "quaternions.h"

// Multiply two quaternions and return a copy of the result, prod = (float)L * R
quaternion_t quat_mult (struct quaternion L, struct quaternion R){
    
    quaternion_t product;
    product.q1 = (float)(L.q1 * R.q1) - (L.q2 * R.q2) - (L.q3 * R.q3) - (L.q4 * R.q4);
    product.q2 = (float)(L.q1 * R.q2) + (L.q2 * R.q1) + (L.q3 * R.q4) - (L.q4 * R.q3);
    product.q3 = (float)(L.q1 * R.q3) - (L.q2 * R.q4) + (L.q3 * R.q1) + (L.q4 * R.q2);
    product.q4 = (float)(L.q1 * R.q4) + (L.q2 * R.q3) - (L.q3 * R.q2) + (L.q4 * R.q1);
    
    return product;
}

// Multiply a reference of a quaternion by a scalar, q = (float)s*q
void quat_scalar(struct quaternion * q, float scalar){
    q -> q1 *= (float)scalar;
    q -> q2 *= (float)scalar;
    q -> q3 *= (float)scalar;
    q -> q4 *= (float)scalar;
}

// Adds two quaternions together and the sum is the pointer to another quaternion, Sum = (float)L + R
void quat_add(struct quaternion * Sum, struct quaternion L, struct quaternion R){
    Sum -> q1 = (float)L.q1 + R.q1;
    Sum -> q2 = (float)L.q2 + R.q2;
    Sum -> q3 = (float)L.q3 + R.q3;
    Sum -> q4 = (float)L.q4 + R.q4;
}

// Subtracts two quaternions together and the sum is the pointer to another quaternion, sum = (float)L - R
void quat_sub(struct quaternion * Sum, struct quaternion L, struct quaternion R){
    Sum -> q1 = (float)L.q1 -  (float)R.q1;
    Sum -> q2 = (float)L.q2 -  (float)R.q2;
    Sum -> q3 = (float)L.q3 -  (float)R.q3;
    Sum -> q4 = (float)L.q4 -  (float)R.q4;
}

// the conjugate of a quaternion is it's imaginary component sign changed  q* = (float)[s, -v] if q = (float)[s, v]
quaternion_t quat_conjugate(struct quaternion q){
    q.q2 = (float)-q.q2;
    q.q3 = (float)-q.q3;
    q.q4 = (float)-q.q4;
    return q;
}

// norm of a quaternion is the same as a complex number
// sqrt( q1^2 + q2^2 + q3^2 + q4^2)
// the norm is also the sqrt(q * conjugate(q)), but thats a lot of operations in the quaternion multiplication
float quat_Norm (quaternion_t q)
{
    return  (float)sqrtf((float)q.q1*q.q1 + (float)q.q2*q.q2 + (float)q.q3*q.q3 +(float)q.q4*q.q4);
}

// Normalizes pointer q by calling quat_Norm(q),
void quat_Normalization(struct quaternion * q){
    float norm = (float)quat_Norm(*q);
    q -> q1 /= (float)norm;
    q -> q2 /= (float)norm;
    q -> q3 /= (float)norm;
    q -> q4 /= (float)norm;
}

// Extends a vector 3D to a quaternion
quaternion_t vec2quat(float* vector){

    quaternion_t extended_vector = {0.0f, vector[0], vector[1], vector[2]};
    return extended_vector;
};

// Returns vector from quaternion
void quat2vec(quaternion_t q, float* vector){

    vector[0] = (float)q.q2;
    vector[1] = (float)q.q3;
    vector[2] = (float)q.q4;
};

/*
Secuencia XYZ
 returns as pointers, roll pitch and yaw from the quaternion
 Assume right hand system
 Roll is about the x axis, represented as phi
 Pitch is about the y axis, represented as theta
 Yaw is about the z axis, represented as psi 
 */
void eulerAngles(struct quaternion q, float* roll, float* pitch, float* yaw){
    
    *roll = (float)atan2f((2.0f*q.q1*q.q2 + 2*q.q3*q.q4), (1.0f - 2.0f*q.q2*q.q2 - 2.0f*q.q3*q.q3)); 
    *pitch = (float)asinf(2.0f*q.q1*q.q3 - 2*q.q2*q.q4);                                  
    *yaw  = (float)atan2f((2.0f*q.q1*q.q4 + 2*q.q2*q.q3), (1.0f - 2.0f*q.q3*q.q3 - 2.0f*q.q4*q.q4));
    
    *yaw *= (float)(180.0f / (float)M_PI);
    *pitch *= (float)(180.0f /(float) M_PI);
    *roll *= (float)(180.0f / (float)M_PI);

}

/*
Secuencia XYZ. Se asumen angulos en radianes
 returns as pointers, roll pitch and yaw from the quaternion
 Assume right hand system
 Roll is about the x axis, represented as phi
 Pitch is about the y axis, represented as theta
 Yaw is about the z axis, represented as psi 
 */
void quat_euler(struct quaternion *q, float* roll, float* pitch, float* yaw){
    
    q->q1 = (float)cosf(*roll/2.0f) * cosf(*pitch/2.0f) * cosf(*yaw/2.0f) + sinf(*roll/2.0f) * sinf(*pitch/2.0f) * sinf(*yaw/2.0f);
    q->q2 = (float)sinf(*roll/2.0f) * cosf(*pitch/2.0f) * cosf(*yaw/2.0f) - cosf(*roll/2.0f) * sinf(*pitch/2.0f) * sinf(*yaw/2.0f);
    q->q3 = (float)cosf(*roll/2.0f) * sinf(*pitch/2.0f) * cosf(*yaw/2.0f) + sinf(*roll/2.0f) * cosf(*pitch/2.0f) * sinf(*yaw/2.0f);
    q->q4 = (float)cosf(*roll/2.0f) * cosf(*pitch/2.0f) * sinf(*yaw/2.0f) - sinf(*roll/2.0f) * sinf(*pitch/2.0f) * cosf(*yaw/2.0f);
}
