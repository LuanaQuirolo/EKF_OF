#include <stdio.h> //TODO: BORRAR

void mulmat(int arows, int acols, int bcols, long double a[arows][acols], long double b[acols][bcols])
{
    int i, j,l;
    long double c[2][2];

    for(i=0; i<arows; i++)
        for(j=0; j<bcols; j++) {
            c[i][j] = 0;
            for(l=0; l<acols; l++){
                c[i][j] += a[i][l] * b[l][j];
            }
        }

    for(i=0; i<arows; i++)
        for(j=0; j<bcols; j++) {
            a[i][j] = c[i][j];
        }
}

int main() {
    // long double a[2][2] = {{0.49857234,0.1230985}, {0.12300980543,0.84359023490874839}};
    // long double b[2][2] = {{0.122335874893,0.14983058346}, {0.789604763456,0.174930570394769075584}};
    // for (int i = 0; i < 10; i++)
    // {
    //     mulmat(2,2,2,a,b);
    // }
    
    
    // for(int i = 0; i < 2; i++) {
    //     for(int j = 0; j < 2; j++) {
    //         printf("%.40Lf ", a[i][j]);
    //     }
    //     printf("\n");
    // }

    long double a = 1.12345678901234567890123456789012345678901234567890123456789012345678901234567890L;
    double b = 1.12345678901234567890123456789012345678901234567890123456789012345678901234567890;
    float c = 1.12345678901234567890123456789012345678901234567890123456789012345678901234567890;
    printf("%.60Lf\n", a);
    printf("%.60f\n", b);
    printf("%.60f\n", c);
    return 0;
}
