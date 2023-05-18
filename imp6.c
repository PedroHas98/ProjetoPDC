#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

#define NX 25 // Dimensão em x
#define NY 25 // Dimensão em y
#define NZ 25 // Dimensão em z

#define DELTA_X 10 // Tamanho da grade em x
#define DELTA_Y 10 // Tamanho da grade em y
#define DELTA_Z 10 // Tamanho da grade em z

#define DELTA_T 0.001 // Passo de tempo
#define NT 2000 // Número total de passos de tempo

#define C 1500 // Velocidade de propagação da onda no meio

#define PI 3.14159265358979323846 // Valor de pi



float ricker_wavelet(float t, float fM) {
    float term1 = 1 - 2 * pow(PI * fM, 2) * pow(t, 2);
    float term2 = exp(-pow(PI * fM, 2) * pow(t, 2));
    return term1 * term2;
}

float s(float t, float x_s, float y_s, float z_s) {
    // Fonte de onda sísmica
    float fM = 10.0; // Frequência de pico da wavelet

    return ricker_wavelet(t, fM);
}

int main() {
    double laplacianoX, laplacianoY, laplacianoZ, laplacianoT; // Laplacianos de x, y e z
    int t, x, y, z, i, j, k; // Índices de tempo e espaço
    double start, end, result;

    // Inicialização das matrizes
    double u[NX][NY][NZ];
    double u_past[NX][NY][NZ];
    double u_future[NX][NY][NZ];

    double alphaX = C * C * DELTA_T * DELTA_T / (DELTA_X * DELTA_X);
    double alphaY = C * C * DELTA_T * DELTA_T / (DELTA_Y * DELTA_Y);
    double alphaZ = C * C * DELTA_T * DELTA_T / (DELTA_Z * DELTA_Z);

    double alphaX2 = alphaX * alphaX;
    double alphaY2 = alphaY * alphaY;
    double alphaZ2 = alphaZ * alphaZ;

    // Inicialização das condições iniciais
    for (i = 0; i < NX; i++) {
        for (j = 0; j < NY; j++) {
            for (k = 0; k < NZ; k++) {
                u_past[i][j][k] = 0.0;
                u[i][j][k] = 0.0;
                u_future[i][j][k] = 0.0;
            }
        }
    }
    start = omp_get_wtime();

    // Simulação da propagação da onda
    for (t = 1; t < NT - 1; t++) {
         #pragma omp parallel for num_threads(8) private(i, j, k, laplacianoX, laplacianoY, laplacianoZ, laplacianoT) shared(u, u_past, u_future) schedule(guided)
        for (i = 2; i < NX - 2; i++) {
            for (j = 2; j < NY - 2; j++) {
                for (k = 2; k < NZ - 2; k++) {
                    float laplacianoX = (-1.0/12.0 * u[i-2][j][k] + 4.0/3.0 * u[i-1][j][k] - 5.0/2.0 * u[i][j][k] + 4.0/3.0 * u[i+1][j][k] - 1.0/12.0 * u[i+2][j][k]) / (DELTA_X * DELTA_X);
                    float laplacianoY = (-1.0/12.0 * u[i][j-2][k] + 4.0/3.0 * u[i][j-1][k] - 5.0/2.0 * u[i][j][k] + 4.0/3.0 * u[i][j+1][k] - 1.0/12.0 * u[i][j+2][k]) / (DELTA_Y * DELTA_Y);
                    float laplacianoZ = (-1.0/12.0 * u[i][j][k-2] + 4.0/3.0 * u[i][j][k-1] - 5.0/2.0 * u[i][j][k] + 4.0/3.0 * u[i][j][k+1] - 1.0/12.0 * u[i][j][k+2]) / (DELTA_Z * DELTA_Z);
                    float laplacianoT = (u_past[i][j][k] - 2.0 * u[i][j][k] + u_future[i][j][k]) / (DELTA_T * DELTA_T);
                    u_future[i][j][k] = 2.0 * u[i][j][k] - u_past[i][j][k] + (1.0 / (C * C)) * (alphaX2 * laplacianoX + alphaY2 * laplacianoY + alphaZ2 * laplacianoZ - laplacianoT) - s(t, DELTA_X, DELTA_Y, DELTA_Z);
                }
            }
        } 
        // Atualização das matrizes u_past, u e u_future
        #pragma omp parallel for num_threads(8) private(i, j, k) shared(u, u_past, u_future) schedule(guided)
        for (i = 2; i < NX - 2; i++) {
            for (j = 2; j < NY - 2; j++) {
                for (k = 2; k < NZ - 2; k++) {
                    u_past[i][j][k] = u[i][j][k];
                    u[i][j][k] = u_future[i][j][k];
                }
            }
        }
        /*FILE* file;
        char name[30];
        if ((t%100)==0) {
                    sprintf(name, "Imp2%d.bin", t);
                    file = fopen(name, "wb");
                    #pragma omp parallel for num_threads(8)
                    for (i = 0; i < NX; i++) {
                        for (j = 0; j < NY; j++) {
                            for (k = 0; k < NZ; k++) {
                                fwrite(&u[i][j][k], sizeof(double), 1, file);
                                
                            }
                        }
                    }
            fclose(file);
        }*/
    }
    end = omp_get_wtime();
    result = end - start;

    printf("Tempo de execução: %f segundos\n", result);

    return 0;
}

