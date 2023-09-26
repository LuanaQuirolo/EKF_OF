quaternion_t quat_mult (struct quaternion L, struct quaternion R){
    
    quaternion_t product;
    product.q1 = (L.q1 * R.q1) - (L.q2 * R.q2) - (L.q3 * R.q3) - (L.q4 * R.q4);
    product.q2 = (L.q1 * R.q2) + (L.q2 * R.q1) + (L.q3 * R.q4) - (L.q4 * R.q3);
    product.q3 = (L.q1 * R.q3) - (L.q2 * R.q4) + (L.q3 * R.q1) + (L.q4 * R.q2);
    product.q4 = (L.q1 * R.q4) + (L.q2 * R.q3) - (L.q3 * R.q2) + (L.q4 * R.q1);
    
    return product;
}

static inline void quat_scalar(struct quaternion * q, float scalar){
    q -> q1 *= scalar;
    q -> q2 *= scalar;
    q -> q3 *= scalar;
    q -> q4 *= scalar;
}

static inline void quat_add(struct quaternion * Sum, struct quaternion L, struct quaternion R){
    Sum -> q1 = L.q1 + R.q1;
    Sum -> q2 = L.q2 + R.q2;
    Sum -> q3 = L.q3 + R.q3;
    Sum -> q4 = L.q4 + R.q4;
}

static inline void quat_sub(struct quaternion * Sum, struct quaternion L, struct quaternion R){
    Sum -> q1 = L.q1 - R.q1;
    Sum -> q2 = L.q2 - R.q2;
    Sum -> q3 = L.q3 - R.q3;
    Sum -> q4 = L.q4 - R.q4;
}

static inline quaternion_t quat_conjugate(struct quaternion q){
    q.q2 = -q.q2;
    q.q3 = -q.q3;
    q.q4 = -q.q4;
    return q;
}

static inline float quat_Norm (struct quaternion q)
{
    return sqrt(q.q1*q.q1 + q.q2*q.q2 + q.q3*q.q3 +q.q4*q.q4);
}

static inline void quat_Normalization(struct quaternion * q){
    float norm = quat_Norm(*q);
    q -> q1 /= norm;
    q -> q2 /= norm;
    q -> q3 /= norm;
    q -> q4 /= norm;
}

quaternion_t quat_extend(double* vector){

    quaternion_t extended_vector;
    extended_vector.q1 = 0;
    extended_vector.q2 = vector[0];
    extended_vector.q3 = vector[1];
    extended_vector.q4 = vector[2];
    
    return extended_vector;
};

void eulerAngles(struct quaternion q, float* roll, float* pitch, float* yaw){
    
    *yaw = atan2f((2*q.q1*q.q2 + 2*q.q3*q.q4), (1 - 2*q.q2*q.q2 - 2*q.q3*q.q3)); 
    *pitch = asinf(2*q.q1*q.q3 - 2*q.q2*q.q4);                                  
    *roll  = atan2f((2*q.q1*q.q4 + 2*q.q2*q.q3), (1 - 2*q.q3*q.q3 - 2*q.q4*q.q4));
    
    *yaw *= (180.0f / M_PI);
    *pitch *= (180.0f / M_PI);
    *roll *= (180.0f / M_PI);

}