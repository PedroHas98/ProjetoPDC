#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define NX 10 // Dimensão em x
#define NY 10 // Dimensão em y
#define NZ 10 // Dimensão em z

#define DELTA_X 0.1 // Tamanho da grade em x
#define DELTA_Y 0.1 // Tamanho da grade em y
#define DELTA_Z 0.1 // Tamanho da grade em z

#define DELTA_T 0.001 // Passo de tempo
#define NT 10 // Número total de passos de tempo

#define C 1500 // Velocidade de propagação da onda no meio

#define PI 3.14159265358979323846 // Valor de pi


float s(float t, float xs, float ys, float zs) {
    // Fonte de onda sísmica
    float f = 25.0; // Frequência de pico da fonte
    float A = pow(2.0 * PI * f, 2.0); //Amplitude da onda 
    float t0 = 1.0 / f; //calcula o tempo de atraso t0 da onda sísmica, que é o tempo necessário para a onda chegar do ponto de origem até o ponto de registro
    float t_delay = t - t0; // calculado o tempo de atraso da onda sísmica
    float r = sqrt(pow(xs, 2.0) + pow(ys, 2.0) + pow(zs, 2.0)); //calcula a distância r entre a fonte de onda sísmica
    return A * pow(t_delay, 2.0) * exp(-A * pow(t_delay, 2.0)) * sin(2.0 * PI * f * (t - r / C - t0));
    /*o valor da função é calculado usando a amplitude da onda sísmica, o tempo de atraso, 
    a distância e a frequência. A função sin é usada para modelar a forma da onda sísmica, 
    enquanto as outras operações são usadas para ajustar o tempo e a amplitude da onda sísmica. 
    O valor calculado é retornado pela função s.*/
}   


int main() {
    double laplacianoX, laplacianoY, laplacianoZ, laplacianoT; // Laplacianos de x, y e z
    int t, x, y, z, i, j, k; // Índices de tempo e espaço

    // Inicializa o campo de pressão com zero
    double u[NX][NY][NZ];
    double u_past[NX][NY][NZ];
    double u_future[NX][NY][NZ];
    double s[NX][NY][NZ];

    double alphaX = C * DELTA_T / DELTA_X;
    double alphaY = C * DELTA_T / DELTA_Y;
    double alphaZ = C * DELTA_T / DELTA_Z;

    double alphaX2 = alphaX * alphaX;
    double alphaY2 = alphaY * alphaY;
    double alphaZ2 = alphaZ * alphaZ;

    // inicialização das condições iniciais
    for (i = 0; i < NX; i++) {
        for (j = 0; j < NY; j++) {
            for (k = 0; k < NZ; k++) {
                u_past[i][j][k] = 0.0;
                u[i][j][k] = 0.0;
                u_future[i][j][k] = 0.0;
                s[i][j][k] = sin(3.14 * i * DELTA_X) * sin(3.14 * j * DELTA_Y) * sin(3.14 * k * DELTA_Z);
            }
        }
    }

    // Simulação da propagação da onda
    for (t = 1; t < NT - 1; t++) {
        printf("t = %d\n", t);
        for (x = 2; x < NX - 2; x++) {
            for (y = 2; y < NY - 2; y++) {
                for (z = 2; z < NZ - 2; z++) {
                    laplacianoX = (-1.0/12.0 * u[i-2][j][k] + 4.0/3.0 * u[i-1][j][k] - 5.0/2.0 * u[i][j][k] + 4.0/3.0 * u[i+1][j][k] - 1.0/12.0 * u[i+2][j][k]) / (DELTA_X * DELTA_X);
                    laplacianoY = (-1.0/12.0 * u[i][j-2][k] + 4.0/3.0 * u[i][j-1][k] - 5.0/2.0 * u[i][j][k] + 4.0/3.0 * u[i][j+1][k] - 1.0/12.0 * u[i][j+2][k]) / (DELTA_Y * DELTA_Y);
                    laplacianoZ = (-1.0/12.0 * u[i][j][k-2] + 4.0/3.0 * u[i][j][k-1] - 5.0/2.0 * u[i][j][k] + 4.0/3.0 * u[i][j][k+1] - 1.0/12.0 * u[i][j][k+2]) / (DELTA_Z * DELTA_Z);
                    laplacianoT = (u_past[i][j][k] - 2.0 * u[i][j][k] + u_future[i][j][k]) / (DELTA_T * DELTA_T);
                    printf("%.2f ", u[x][y][z]);
                }
                 printf("\n");
            }
             printf("\n");
        }
        u_future[i][j][k] = 2.0 * u[i][j][k] - u_past[i][j][k] + alphaX2 * laplacianoX + alphaY2 * laplacianoY + alphaZ2 * laplacianoZ + laplacianoT;
    }


return 0;
}
