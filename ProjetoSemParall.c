#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <locale.h>


#define NX 25 // Dimensão em x
#define NY 25 // Dimensão em y
#define NZ 25 // Dimensão em z

#define DELTA_X 10 // Tamanho da grade em x
#define DELTA_Y 10 // Tamanho da grade em y
#define DELTA_Z 10 // Tamanho da grade em z

#define DELTA_T 0.001 // Passo de tempo
#define NT 1000 // Número total de passos de tempo

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

float calcularMedia(double vetor[], int tamanho) {
    double soma = 0;
    int i;
    for (i = 0; i < tamanho; i++) {
        soma += vetor[i];
    }

    return soma / tamanho;
}

int main() {
    setlocale(LC_ALL, "pt_BR");
    double laplacianoX, laplacianoY, laplacianoZ, laplacianoT; // Laplacianos de x, y e z
    int t, x, y, z, i, j, k; // Índices de tempo e espaço
    clock_t start, end;
    double result;

    // Inicialização das matrizes
    double u[NX][NY][NZ];
    double u_past[NX][NY][NZ];
    double u_future[NX][NY][NZ];

    double alphaX = C * C * DELTA_T * DELTA_T / (DELTA_X * DELTA_X);
    double alphaY = C * C * DELTA_T * DELTA_T / (DELTA_Y * DELTA_Y);
    double alphaZ = C * C * DELTA_T * DELTA_T / (DELTA_Z * DELTA_Z);

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
    int g;
    double cal_mediana [5];
    for(g=0; g<5; g++){

    start = clock();

    // Simulação da propagação da onda
    for (t = 1; t < NT - 1; t++) {
        for (i = 2; i < NX - 2; i++) {
            for (j = 2; j < NY - 2; j++) {
                for (k = 2; k < NZ - 2; k++) {
                    float laplacianoX = (-1.0/12.0 * u[i-2][j][k] + 4.0/3.0 * u[i-1][j][k] - 5.0/2.0 * u[i][j][k] + 4.0/3.0 * u[i+1][j][k] - 1.0/12.0 * u[i+2][j][k]) / (DELTA_X * DELTA_X);
                    float laplacianoY = (-1.0/12.0 * u[i][j-2][k] + 4.0/3.0 * u[i][j-1][k] - 5.0/2.0 * u[i][j][k] + 4.0/3.0 * u[i][j+1][k] - 1.0/12.0 * u[i][j+2][k]) / (DELTA_Y * DELTA_Y);
                    float laplacianoZ = (-1.0/12.0 * u[i][j][k-2] + 4.0/3.0 * u[i][j][k-1] - 5.0/2.0 * u[i][j][k] + 4.0/3.0 * u[i][j][k+1] - 1.0/12.0 * u[i][j][k+2]) / (DELTA_Z * DELTA_Z);
                    u_future[i][j][k] = 2.0 * u[i][j][k] - u_past[i][j][k] + (1.0 / (C * C)) * (alphaX * laplacianoX + alphaY * laplacianoY + alphaZ * laplacianoZ - laplacianoT) - s(t, DELTA_X, DELTA_Y, DELTA_Z);

                }
            }
        }
            
        // Atualização das matrizes u_past, u e u_future
        for (i = 2; i < NX - 2; i++) {
            for (j = 2; j < NY - 2; j++) {
                for (k = 2; k < NZ - 2; k++) {
                    u_past[i][j][k] = u[i][j][k];
                    u[i][j][k] = u_future[i][j][k];
                }
            }
        }

    }
    end = clock();
    
    result = ((double)(end - start)) / CLOCKS_PER_SEC;
    cal_mediana[g] = result;
    printf("Execucao %d Tempo: %f segundos\n",g+1, result);
    
    }
    printf("Media de tempo: %f",calcularMedia(cal_mediana,5));
    return 0;
}
