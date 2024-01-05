/***************************************************************
 *                  Cuaterniones
 *                  Creado: 10 sep. 2023
 *                  Autor: Luana Quirolo
 *                  Padron: 102102
****************************************************************/

#include "quaternions.h"

// Multiply two quaternions and return a copy of the result, prod = L * R
quaternion_t quat_mult (struct quaternion L, struct quaternion R){
    
    quaternion_t product;
    product.q1 = (L.q1 * R.q1) - (L.q2 * R.q2) - (L.q3 * R.q3) - (L.q4 * R.q4);
    product.q2 = (L.q1 * R.q2) + (L.q2 * R.q1) + (L.q3 * R.q4) - (L.q4 * R.q3);
    product.q3 = (L.q1 * R.q3) - (L.q2 * R.q4) + (L.q3 * R.q1) + (L.q4 * R.q2);
    product.q4 = (L.q1 * R.q4) + (L.q2 * R.q3) - (L.q3 * R.q2) + (L.q4 * R.q1);
    
    return product;
}

// Multiply a reference of a quaternion by a scalar, q = s*q
void quat_scalar(struct quaternion * q, double scalar){
    q -> q1 *= scalar;
    q -> q2 *= scalar;
    q -> q3 *= scalar;
    q -> q4 *= scalar;
}

// Adds two quaternions together and the sum is the pointer to another quaternion, Sum = L + R
void quat_add(struct quaternion * Sum, struct quaternion L, struct quaternion R){
    Sum -> q1 = L.q1 + R.q1;
    Sum -> q2 = L.q2 + R.q2;
    Sum -> q3 = L.q3 + R.q3;
    Sum -> q4 = L.q4 + R.q4;
}

// Subtracts two quaternions together and the sum is the pointer to another quaternion, sum = L - R
void quat_sub(struct quaternion * Sum, struct quaternion L, struct quaternion R){
    Sum -> q1 = L.q1 - R.q1;
    Sum -> q2 = L.q2 - R.q2;
    Sum -> q3 = L.q3 - R.q3;
    Sum -> q4 = L.q4 - R.q4;
}

// the conjugate of a quaternion is it's imaginary component sign changed  q* = [s, -v] if q = [s, v]
quaternion_t quat_conjugate(struct quaternion q){
    q.q2 = -q.q2;
    q.q3 = -q.q3;
    q.q4 = -q.q4;
    return q;
}

// norm of a quaternion is the same as a complex number
// sqrt( q1^2 + q2^2 + q3^2 + q4^2)
// the norm is also the sqrt(q * conjugate(q)), but thats a lot of operations in the quaternion multiplication
double quat_Norm (quaternion_t q)
{
    return sqrt(q.q1*q.q1 + q.q2*q.q2 + q.q3*q.q3 +q.q4*q.q4);
}

// Normalizes pointer q by calling quat_Norm(q),
void quat_Normalization(struct quaternion * q){
    double norm = quat_Norm(*q);
    q -> q1 /= norm;
    q -> q2 /= norm;
    q -> q3 /= norm;
    q -> q4 /= norm;
}

// Extends a vector 3D to a quaternion
quaternion_t vec2quat(double* vector){

    quaternion_t extended_vector = {0, vector[0], vector[1], vector[2]};
    return extended_vector;
};

// Returns vector from quaternion
void quat2vec(quaternion_t q, double* vector){

    vector[0] = q.q2;
    vector[1] = q.q3;
    vector[2] = q.q4;
};

/*
 returns as pointers, roll pitch and yaw from the quaternion
 Assume right hand system
 Roll is about the x axis, represented as phi
 Pitch is about the y axis, represented as theta
 Yaw is about the z axis, represented as psi 
 */
void eulerAngles(struct quaternion q, float* roll, float* pitch, float* yaw){
    
    *roll = atan2f((2*q.q1*q.q2 + 2*q.q3*q.q4), (1 - 2*q.q2*q.q2 - 2*q.q3*q.q3)); 
    *pitch = asinf(2*q.q1*q.q3 - 2*q.q2*q.q4);                                  
    *yaw  = atan2f((2*q.q1*q.q4 + 2*q.q2*q.q3), (1 - 2*q.q3*q.q3 - 2*q.q4*q.q4));
    
    *yaw *= (180.0f / M_PI);
    *pitch *= (180.0f / M_PI);
    *roll *= (180.0f / M_PI);

}
