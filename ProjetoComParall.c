#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

#define NX 25 // Dimensão em x
#define NY 25 // Dimensão em y
#define NZ 25 // Dimensão em z
#define DELTAX 10 // Tamanho da grade em x
#define DELTAY 10 // Tamanho da grade em y
#define DELTAZ 10 // Tamanho da grade em z
#define DELTAT 0.001 // Passo de tempo
#define NT 1000 // Número total de passos de tempo
#define C 1500 // Velocidade de propagação da onda no meio
#define PI 3.14159265358979323846 // Valor de pi

float ricker_wavelet(float t, float fM) {
    float term1 = 1 - 2 * pow(PI * fM, 2) * pow(t, 2);
    float term2 = exp(-pow(PI * fM, 2) * pow(t, 2));
    return term1 * term2;
}

float s(float t, float xs, float ys, float zs) {
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
    double upast[NX][NY][NZ];
    double ufuture[NX][NY][NZ];

    double alphaX = C * C * DELTAT * DELTAT / (DELTAX * DELTAX);
    double alphaY = C * C * DELTAT * DELTAT / (DELTAY * DELTAY);
    double alphaZ = C * C * DELTAT * DELTAT / (DELTAZ * DELTAZ);

    // Inicialização das condições iniciais
    for (i = 0; i < NX; i++) {
        for (j = 0; j < NY; j++) {
            for (k = 0; k < NZ; k++) {
                upast[i][j][k] = 0.0;
                u[i][j][k] = 0.0;
                ufuture[i][j][k] = 0.0;
            }
        }
    }

    // Loop para executar o programa para diferentes valores de OMP_SCHEDULE
    const char *omp_schedules[] = {"static,8", "dynamic,8", "guided,8"};
    int num_schedules = sizeof(omp_schedules) / sizeof(omp_schedules[0]);
    
    int num_thread[] = {1, 2, 4, 8};
    int numT = sizeof(num_thread) / sizeof(num_thread[0]);
    
    for (int g = 0; g < numT; g++) {
        for (int schedule_idx = 0; schedule_idx < num_schedules; schedule_idx++) {
            const char *omp_schedule = omp_schedules[schedule_idx];
            setenv("OMP_SCHEDULE", omp_schedule, 1);
            printf("Executing with OMP_SCHEDULE=%s\n", omp_schedule);
            int thread = num_thread[g];
            
            #pragma omp num_threads(thread) parallel for private(t, i, j, k, laplacianoX, laplacianoY, laplacianoZ) shared(u, upast, ufuture, DELTAX, DELTAY, DELTAZ, alphaZ, alphaY, alphaX) schedule(runtime) 
                #pragma omp barrier
                start = omp_get_wtime();  
                // Simulação da propagação da onda
                    for (t = 1; t < NT - 1; t++) {
                    
                        for (i = 2; i < NX - 2; i++) {
                            for (j = 2; j < NY - 2; j++) {
                                for (k = 2; k < NZ - 2; k++) {
                                    laplacianoX = (-1.0 / 12.0 * u[i - 2][j][k] + 4.0 / 3.0 * u[i - 1][j][k] - 5.0 / 2.0 * u[i][j][k] + 4.0 / 3.0 * u[i + 1][j][k] - 1.0 / 12.0 * u[i + 2][j][k]) / (DELTAX * DELTAX);
                                    laplacianoY = (-1.0 / 12.0 * u[i][j - 2][k] + 4.0 / 3.0 * u[i][j - 1][k] - 5.0 / 2.0 * u[i][j][k] + 4.0 / 3.0 * u[i][j + 1][k] - 1.0 / 12.0 * u[i][j + 2][k]) / (DELTAY * DELTAY);
                                    laplacianoZ = (-1.0 / 12.0 * u[i][j][k - 2] + 4.0 / 3.0 * u[i][j][k - 1] - 5.0 / 2.0 * u[i][j][k] + 4.0 / 3.0 * u[i][j][k + 1] - 1.0 / 12.0 * u[i][j][k + 2]) / (DELTAZ * DELTAZ);
                                    ufuture[i][j][k] = 2.0 * u[i][j][k] - upast[i][j][k] + (1.0 / (C * C)) * (alphaX * laplacianoX + alphaY * laplacianoY + alphaZ * laplacianoZ) - s(t, DELTAX, DELTAY, DELTAZ);
                                }
                            }
                        }

                        // Atualização das matrizes upast, u e ufuture
                        for (i = 2; i < NX - 2; i++) {
                            for (j = 2; j < NY - 2; j++) {
                                for (k = 2; k < NZ - 2; k++) {
                                    upast[i][j][k] = u[i][j][k];
                                    u[i][j][k] = ufuture[i][j][k];
                                }
                            }
                        }
                    
                    }
                    end = omp_get_wtime();
                    result = end - start;
                
            printf("Tempo de execução: %f segundos|| NumT: %i\n", result, thread);
            printf("\n");
        }
    }

    return 0;
}
