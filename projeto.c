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
#define NT 100 // Número total de passos de tempo

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
    float u[NT][NX][NY][NZ]; // Campo de pressão
    float laplacianoX, laplacianoY, laplacianoZ; // Laplacianos de x, y e z
    int t, x, y, z; // Índices de tempo e espaço

    // Inicializa o campo de pressão com zero
    for (t = 0; t < NT; t++) {
        for (x = 0; x < NX; x++) {
            for (y = 0; y < NY; y++) {
                for (z = 0; z < NZ; z++) {
                    u[t][x][y][z] = 0.0;
                }
            }
        }
    }

    // Simulação da propagação da onda
    for (t = 1; t < NT - 1; t++) {
        for (x = 2; x < NX - 2; x++) {
            for (y = 2; y < NY - 2; y++) {
                for (z = 2; z < NZ - 2; z++) {
                    // Cálculo do laplaciano de x
                    laplacianoX = (u[t][x-2][y][z] - 4.0 * u[t][x-1][y][z] + 6.0 * u[t][x][y][z] - 4.0 * u[t][x+1][y][z] + u[t][x+2][y][z]) / pow(DELTA_X, 2.0);
                    // Cálculo do laplaciano de y
                    laplacianoY = (u[t][x][y-2][z] - 4.0 * u[t][x][y-1][z] + 6.0 * u[t][x][y][z] - 4.0 * u[t][x][y+1][z] + u[t][x][y+2][z]) / pow(DELTA_Y, 2.0);
                    // Cálculo do laplaciano de z
                    laplacianoZ = (u[t][x][y][z-2] - 4.0 * u[t][x][y][z-1] + 6.0 * u[t][x][y][z] - 4.0 * u[t][x][y][z+1] + u[t][x][y][z+2]) / pow(DELTA_Z, 2.0);
                    // Cálculo do campo de pressão no próximo passo de tempo
                    u[t+1][x][y][z] = 2.0 * u[t][x][y][z] - u[t-1][x][y][z] + pow(C * DELTA_T, 2.0) * (laplacianoX + laplacianoY + laplacianoZ + s(t * DELTA_T, x * DELTA_X, y * DELTA_Y, z * DELTA_Z));
                }
            }
        }
    }

    // Impressão do campo de pressão em um plano XY
    printf("Campo de pressão em um plano XY:\n");
    int i, j;
    for (i = 0; i < NX; i++) {
        for (j = 0; j < NY; j++) {
            printf("%.3f ", u[NT-1][i][j][NZ/2]);
        }
    printf("\n");
}

return 0;
}
