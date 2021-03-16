#include <stdio.h>
#include <stdlib.h>

#define NX 64
#define NY 64
#define N (NX * NY)
#define NBN 4
#define NSTEPS 100
#define XLOOP(i) ((i + NX) % NX)
#define YLOOP(i) ((i + NY) % NY)
double matrix[NSTEPS][NX][NY];
int nblist[N][NBN];

void Init() {
    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NY; j++) {
            matrix[0][i][j] = ((double)rand() / RAND_MAX) * 0.1;
            int idx = i + j * NX;
            int lidx = XLOOP(i - 1) + j * NX, ridx = XLOOP(i + 1) + j * NX,
                uidx = i + YLOOP(j - 1) * NX, didx = i + YLOOP(j + 1) * NX;
            nblist[idx][0] = lidx; nblist[idx][1] = ridx;
            nblist[idx][2] = uidx; nblist[idx][3] = didx;
        }
    }
}

void Show(int t) {
    printf("------------------------------------------------\n");
    printf("The Step is %d\n", t);
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            printf("%lf ", matrix[t][i][j]);
            if (j == 4 - 1) printf("\n");
        }
    }
}

void Compute(int i, int t) {
    int ix = i % NX, iy = i / NX;
    double sum = matrix[t][ix][iy];
    for (int j = 0; j < NBN; j++) {
        int nbj = nblist[i][j];
        int jx = j % NX, jy = j / NX;
        sum += matrix[t][jx][jy];
    } 
    matrix[t + 1][ix][iy] = sum / 5.0;
}

int main() {
    Init(); 
    for (int t = 0; t < NSTEPS; t++) {
        Show(t);
        for (int i = 0; i < N; i++)
            Compute(i, t);
    }
    return 0;
}
