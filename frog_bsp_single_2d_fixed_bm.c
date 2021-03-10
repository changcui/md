#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include <string.h>
// Simulation Size
#define NX 4
#define NY 4
#define N (NX * NY)
#define NPROCS N
#define NSTEPS 10000
#define TAG 0
#define MASTER 0

// Simulation Params
#define KB 0.0000115
#define KA 0.00596
#define KC 0.7722087
#define XBALANCE 0.396047
#define YBALANCE 0.396047
#define XCYCLE (NX*XBALANCE)
#define YCYCLE (NY*YBALANCE)
#define MIND 0.353
#define T 0.007
#define S 0
#define A 1
#define V 2
#define X 0
#define Y 1
#define NBN 4

// Compute Marcos
#define L(i) (((i % NX) - 1 + NX) % NX + i / NX * NX)
#define R(i) (((i % NX) + 1 + NX) % NX + i / NX * NX)
#define U(i) (i % NX + ((i / NX - 1 + NY) % NY) * NX)
#define D(i) (i % NX + ((i / NX + 1 + NY) % NY) * NX)
#define copy(x, y) memcpy(x, y, sizeof(x))
#define zero(x) memset(x, 0, sizeof(x))
#define si s[i]
#define sit s[i][t]
#define sjt s[j][t]

typedef double Vec[2];
typedef Vec State[3];
State s[N][NSTEPS+1];
int nbl[N][NBN];

void Add(Vec v1, Vec v2, Vec res) {res[X] = v1[X] + v2[X]; res[Y] = v1[Y] + v2[Y];}
void Sub(Vec v1, Vec v2, Vec res) {res[X] = v1[X] - v2[X]; res[Y] = v1[Y] - v2[Y];}
void Mul(Vec v1, double c, Vec res) {res[X] = v1[X] * c; res[Y] = v1[Y] * c;}
double Dot(Vec v1, Vec v2) {return v1[X] * v2[X] + v1[Y] * v2[Y];}
void _slp(double &s, double cycle) {if(s < 0) s += cycle; if(s > cycle) s -= cycle;}
void SLoop(Vec v) {_slp(v[X], XCYCLE); _slp(v[Y], YCYCLE);}
void _dlp(double &d, double cycle) {if(d > cycle/2) d -= cycle; else if(d < -cycle/2) d += cycle; }
void DLoop(Vec v) {_dlp(v[X], XCYCLE); _dlp(v[Y], YCYCLE);}
double Distance(int i, int j, int t) {
    Vec dij;
    Sub(sjt[S], sit[S], dij);
    DLoop(dij);
    return sqrt(Dot(dij, dij));
}

double Vdw(double d) {
    if (!d) return 0;
    double d2 = d * d, d4 = d2 * d2, d6 = d2 * d4, d7 = d * d6, d13 = d6 * d7;
    return 6.0 * KA / d7 - 12.0 * KB / d13;
}

double _ed(double d) {
    if (!d) return 0;
    double d2 = d * d, d4 = d2 * d2, d6 = d2 * d4, d12 = d6 * d6;
    return -KA / d6 + KB / d12 + KC;
}

double Edi(int i, int t) {
    double e = 0.0;
    for (int j = 0; j < NBN; j++) {
        int idx = nbl[i][j];
        double dist = Distance(i, idx, t);
        e += _ed(dist);
    }
    return e / 2.0;
}

double Evi(int i, int t) {
    return 0.5 * Dot(sit[V], sit[V]);
}

void Balance(Vec v1, Vec v2) {
    v1[X] = (v1[X] - v2[X]) / 2;
    v2[X] = -v1[X];
    v1[Y] = (v1[Y] - v2[Y]) / 2;
    v2[Y] = -v2[Y];
}

void Aij(int i, int j, int t, Vec a) {
    int idx = nbl[i][j];
    Vec dij;
    Sub(s[idx][t][S], sit[S], dij);
    DLoop(dij);
    double ds = sqrt(Dot(dij, dij)), f = Vdw(ds);
    Mul(dij, f/ds, a);
}

void Init() {
    for (int i = 0; i < N; i++) {
        int ix = i % NX, iy = i / NX;
        si[0][S][X] = ix * XBALANCE;
        if (ix % 2 != 0) si[0][S][X] = s[(ix-1)+iy*NX][0][S][X] + MIND + 0.005*(ix%10-1);
        si[0][S][Y] = iy * (MIND + 0.0215);
        if (ix % 2 != 0) si[0][S][Y] = s[(ix-1)+iy*NX][0][S][Y] + 0.5*MIND + 0.005*(ix%10-1);
        nbl[i][0] = L(i); nbl[i][1] = R(i); nbl[i][2] = U(i); nbl[i][3] = D(i);
    } 
    for (int i = 0; i < N; i++) {
        zero(si[0][A]);
        zero(si[0][V]);
        for (int j = 0; j < NBN; j++) {
            Vec a;
            Aij(i, j, 0, a);
            Add(a, si[0][A], si[0][A]);
        } 
    }
}

void Frog(int i, int t) {
    Vec tmp;
    Mul(si[t][V], T, tmp);
    Add(tmp, si[t][S], si[t+1][S]);
    Mul(si[t][A], T*T/2, tmp);
    Add(tmp, si[t+1][S], si[t+1][S]);
    SLoop(si[t+1][S]);

    zero(si[t+1][A]);
    for (int j = 0; j < NBN; j++) {
        Vec a;
        Aij(i, j, t, a);
        Add(a, si[t+1][A], si[t+1][A]);
    }
    Add(sit[A], si[t+1][A], si[t+1][V]);
    Mul(si[t+1][V], T/2, si[t+1][V]);
    Add(sit[V], si[t+1][V], si[t+1][V]);
}

void Comm(int i, int t) {
    MPI_Request reqs[4];
    int lpid = L(i), rpid = R(i), upid = U(i), dpid = D(i);
    MPI_Isend(&s[i][t+1], sizeof(State), MPI_DOUBLE, lpid, TAG, MPI_COMM_WORLD, &reqs[0]);
    MPI_Isend(&s[i][t+1], sizeof(State), MPI_DOUBLE, rpid, TAG, MPI_COMM_WORLD, &reqs[1]);
    MPI_Isend(&s[i][t+1], sizeof(State), MPI_DOUBLE, upid, TAG, MPI_COMM_WORLD, &reqs[2]);
    MPI_Isend(&s[i][t+1], sizeof(State), MPI_DOUBLE, dpid, TAG, MPI_COMM_WORLD, &reqs[3]);
    MPI_Recv(&s[lpid][t+1], sizeof(State), MPI_DOUBLE, lpid, TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&s[rpid][t+1], sizeof(State), MPI_DOUBLE, rpid, TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&s[upid][t+1], sizeof(State), MPI_DOUBLE, upid, TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&s[dpid][t+1], sizeof(State), MPI_DOUBLE, dpid, TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Waitall(4, reqs, MPI_STATUS_IGNORE);
}

int main() {
    MPI_Init(NULL, NULL);
    int my_pid, nprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_pid);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    Init();
    double start_time = MPI_Wtime();
    for (int t = 0; t < NSTEPS; t++) {
        Frog(my_pid, t);
        Comm(my_pid, t);
    }
    double end_time = MPI_Wtime();
    if (my_pid == 0) {
        printf("Success!\n");
        printf("Exec time is %lfs\n", end_time - start_time);
    }
    MPI_Finalize();
    return 0;
}
