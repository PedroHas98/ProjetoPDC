#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define NX 10
#define NY 10
#define NZ 10
#define C 1500.0
#define DT 0.001
#define DX 0.01
#define DY 0.01
#define DZ 0.01

void calculate_pressure_field(double u[][NY][NZ], double s[][NY][NZ]) {
    int t, x, y, z;
    double laplacian;
    double u_old[NX][NY][NZ];
    double u_new[NX][NY][NZ];

    // Inicialização do campo de pressão
    for (x = 0; x < NX; x++) {
        for (y = 0; y < NY; y++) {
            for (z = 0; z < NZ; z++) {
                u[x][y][z] = 0.0;
            }
        }
    }

    // Condições de contorno refletoras em x, y e z
    for (t = 0; t < 2; t++) {
        for (y = 0; y < NY; y++) {
            for (z = 0; z < NZ; z++) {
                u[t][y][z] = 0.0;
                u[NX - 1 - t][y][z] = 0.0;
            }
        }

        for (x = 0; x < NX; x++) {
            for (z = 0; z < NZ; z++) {
                u[x][t][z] = 0.0;
                u[x][NY - 1 - t][z] = 0.0;
            }
        }

        for (x = 0; x < NX; x++) {
            for (y = 0; y < NY; y++) {
                u[x][y][t] = 0.0;
                u[x][y][NZ - 1 - t] = 0.0;
            }
        }
    }

    // Cálculo do campo de pressão em cada passo de tempo
    for (t = 2; t < 5000; t++) {
        // Cópia do campo de pressão no passo anterior
        for (x = 0; x < NX; x++) {
            for (y = 0; y < NY; y++) {
                for (z = 0; z < NZ; z++) {
                    u_old[x][y][z] = u[x][y][z];
                }
            }
        }

        // Cálculo do campo de pressão no novo passo de tempo
        for (x = 2; x < NX - 2; x++) {
            for (y = 2; y < NY - 2; y++) {
                for (z = 2; z < NZ - 2; z++) {
                    laplacian = (-u_old[x][y-2][z] / 12.0 + 4.0 * u_old[x][y-1][z] / 3.0 - 5.0 * u_old[x][y][z] / 2.0 + 4.0 * u_old[x][y+1][z] / 3.0 - u_old[x][y+2][z] / 12.0) / pow(DX, 2)
                                + (-u_old[x-2][y][z] / 12.0 + 4.0 * u_old[x-1][y][z] / 3.0 - 5.0 * u_old[x][y][z] / 2.0 + 4.0 * u_old[x+1][y][z] / 3.0 - u_old[x+2][y][z] / 12.0) / pow(DY, 2)
                                + (-u_old[x][y][z-2] / 12.0 + 4.0 * u_old[x][y][z-1] / 3.0 - 5.0 * u_old[x][y][z] / 2.0 + 4.0 * u_old[x][y][z+1] / 3.0 - u_old[x][y][z+2] / 12.0) / pow(DZ, 2);
                                u_new[x][y][z] = 2.0 * u[x][y][z] - u_old[x][y][z] + pow(C * DT, 2) * laplacian + DT * DT * s[x][y][z];
                    // Exibição da propagação da onda
                    if (t >= 1000 && (x + y + z) % 10 == 0) {
                        printf("t = %d, x = %d, y = %d, z = %d, u = %f\n", t, x, y, z, u_new[x][y][z]);
                    }
                }
            }
        }

        // Atualização do campo de pressão
        for (x = 0; x < NX; x++) {
            for (y = 0; y < NY; y++) {
                for (z = 0; z < NZ; z++) {
                    u[x][y][z] = u_new[x][y][z];
                }
            }
        }
    }
}

int main() {
    double u[NX][NY][NZ];
    double s[NX][NY][NZ];
    // Fonte sonora
    int x, y, z, t0 = 1000, sigma = 10;
    double amplitude = 1.0;
    for (t0 = 1000; t0 < 1500; t0++) {
        for (x = 0; x < NX; x++) {
            for (y = 0; y < NY; y++) {
                for (z = 0; z < NZ; z++) {
                    s[x][y][z] += amplitude * exp(-(pow(x - NX / 2, 2) + pow(y - NY / 2, 2) + pow(z - NZ / 2, 2)) / (2.0 * pow(sigma, 2))) * pow(sin(2.0 * 3.14 * (t0 - 1000) * DT * 10.0), 2);
                }
            }
        }
    }

    // Cálculo do campo de pressão
    calculate_pressure_field(u, s);

return 0;

}

