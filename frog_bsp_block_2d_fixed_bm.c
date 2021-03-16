#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
const int TAG = 0;
const int MASTER = 0;
const int NITERS = 2;
#define XLOOP(i) ((i + NPX) % NPX)
#define YLOOP(i) ((i + NPY) % NPY)

#define KB   0.0000115
#define KA   0.00596
#define XBALANCE    0.396047
#define YBALANCE    0.396047//10
#define MIND    0.353
#define NX 64
#define NY 64
#define N (NX*NY)
#define NPX 4
#define NPY 4
#define NPROCS (NPX*NPY)
#define NLX (NX/NPX)
#define NLY (NY/NPY)
#define NATOMS (NLX*NLY)
#define XCYCLE  (NX*XBALANCE)
#define YCYCLE  (NY*YBALANCE)
#define SIGN(x) (x>0?1:-1)
#define MAX(x, y) (x>y?x:y)
#define MIN(x, y) (x<y?x:y)
#define L(i) (((i%NX)-1+NX)%NX + (i/NX)*NX)
#define R(i) (((i%NX)+1+NX)%NX + (i/NX)*NX)
#define U(i) ((i%NX) + (((i/NX)+1+NY)%NY)*NX)
#define D(i) ((i%NX) + (((i/NX)-1+NY)%NY)*NX)
#define copy(x, y) memcpy(x, y, sizeof(x))
#define NSTEPS 100
#define T 0.007//0.011-frog
#define S 0
#define A 1
#define V 2
#define DV 3
#define X 0
#define Y 1
#define NBN 4
#define NBATOMS (NATOMS+2*(NLX+NLY))
#define si s[i]
#define sit s[i][t]
typedef double Vec[2];
typedef Vec State[4];
State s[N][NSTEPS+1], b[N][NSTEPS+1];
Vec aij[N][NSTEPS+1][NBN], aji[N][NSTEPS+1][NBN];
Vec a0[N][NSTEPS+1], v0[N][NSTEPS+1]; int jofi[N][NBN];
int nbl[N][NBN], nbj_idx[N][NBN], atomid[NPROCS][NBATOMS];

double deij[N][NSTEPS+1][NBN]={0};
Vec ndv[N][NSTEPS+1][NBN]={0}, mdv[N][NSTEPS+1][NBN]={0};
Vec dvij[N][NSTEPS+1][NBN]={0}, dvji[N][NSTEPS+1][NBN]={0};

void zero(Vec v){v[X]=0; v[Y]=0;}
void Add(Vec v1, Vec v2, Vec res){ res[X] = v1[X]+v2[X]; res[Y] = v1[Y]+v2[Y];}
void Sub(Vec v1, Vec v2, Vec res){ res[X] = v1[X]-v2[X]; res[Y] = v1[Y]-v2[Y];}
void Mul(Vec v1, double c, Vec res){ res[X] = v1[X]*c; res[Y] = v1[Y]*c;}
double Dot(Vec v1, Vec v2){ return v1[X]*v2[X] + v1[Y]*v2[Y];}
void slp(double &s, double cycle) {if(s<0) s+=cycle; if(s>cycle) s-=cycle;}
void SLoop(Vec v){slp(v[X], XCYCLE); slp(v[Y], YCYCLE);}
void dlp(double &d, double cycle) {if(d>cycle/2) d-=cycle; else if(d<-cycle/2) d+=cycle; }
void DLoop(Vec v){dlp(v[X], XCYCLE); dlp(v[Y], YCYCLE);}
double Dist_ij(int i, int j, int t){ Vec dij; Sub(b[j][t][S],sit[S],dij); DLoop(dij); return sqrt(Dot(dij,dij));}
double Vdw(double d){if(!d) return 0; double d2=d*d,d4=d2*d2,d6=d2*d4,d7=d*d6,d13=d6*d7; return 6.0*KA/d7-12.0*KB/d13;}
double Ed(double d) {if(!d) return 0; double d2=d*d,d4=d2*d2,d6=d2*d4,d12=d6*d6; return -KA/d6+KB/d12+0.7722087;}
double Edi(int i, int t){ double e=0; for(int j=0;j<NBN;j++){int nbj=nbl[i][j]; e += Ed(Dist_ij(i,nbj,t));} return e/2;}
double Evi(int i, int t){ return 0.5*Dot(sit[V], sit[V]);}
double Esum(int t) {double Er=0; for (int i=0; i<N; i++) Er+=Edi(i,t)+Evi(i,t); return Er;}
double M_sum(int t){double mom=0; for(int i=0; i<N; i++) mom += sit[V][X] + sit[V][Y]; return mom;}
void Balance_A(Vec v1, Vec v2) {v1[X]=(v1[X]-v2[X])/2; v2[X]=-v1[X]; v1[Y]=(v1[Y]-v2[Y])/2; v2[Y]=-v1[Y];}
int inblock(int j_idx) {return (j_idx<NATOMS);}
void AIJ(Vec nbs, Vec selfs, Vec a){
        Vec dij; Sub(nbs, selfs, dij); DLoop(dij); double ds=sqrt(Dot(dij, dij)), f=Vdw(ds); Mul(dij, f/ds, a);}
void compute_ss(int t, int t0, Vec s, Vec a, Vec v, int i, int ifin, Vec st1, Vec st2, Vec at) {
    if (!ifin) {
        Vec ss, tmp; double dt = (t-t0)*T;
        Mul(v, dt, tmp); Add(tmp, s, ss); Mul(a, dt*dt/2, tmp); Add(tmp, ss, ss); SLoop(ss);
        AIJ(ss, st2, at);
    }
    else AIJ(st1, st2, at);
}
void new_Aji(int i, int j, int t, int t0){
    int nbj=nbl[i][j];
    compute_ss(t, t0, s[i][t0][S], a0[i][t0], v0[i][t0], i, inblock(nbj_idx[i][j]), sit[S], b[nbj][t][S], aji[i][t][j]);
}
void new_Aij(int i, int j, int t, int t0){
    int nbj=nbl[i][j];
    Vec snbj, tmp; Mul(s[nbj][t-1][V],T,tmp); Add(tmp,s[nbj][t-1][S],snbj); Mul(s[nbj][t-1][A],T*T/2,tmp); Add(tmp,snbj,snbj); SLoop(snbj);
    compute_ss(t, t0, b[nbj][t0][S], b[nbj][t0][A], b[nbj][t0][V], i, inblock(nbj_idx[i][j]), snbj, sit[S], aij[i][t][j]);
}
void init(){
    for(int i=0;i<N;i++) { int ix = i%NX, iy = i/NX;
        si[0][S][X] = ix*XBALANCE; if(ix%2!=0) si[0][S][X]=s[(ix-1)+iy*NX][0][S][X]+MIND+0.005*(ix%10-1);
        si[0][S][Y] = iy*(MIND+0.0215); if(ix%2!=0) si[0][S][Y]=s[(ix-1)+iy*NX][0][S][Y]+0.5*MIND+0.005*(ix%10-1);
        nbl[i][0]=L(i); nbl[i][1]=R(i); nbl[i][2]=U(i); nbl[i][3]=D(i);
    }
    for(int i=0;i<N;i++) for(int j=0;j<NBN;j++) {int nbj=nbl[i][j]; for(int k=0;k<NBN;k++) if(nbl[nbj][k]==i){jofi[i][j]=k; break;}}
    for(int i=0;i<N;i++) {
        zero(si[0][A]);
        for(int j=0;j<NBN;j++) {AIJ(s[nbl[i][j]][0][S], si[0][S], aij[i][0][j]); Add(aij[i][0][j], si[0][A], si[0][A]);}
    }
    for(int i=0;i<N;i++) {zero(si[0][V]); copy(a0[i][0], si[0][A]);  copy(v0[i][0], si[0][V]);}
    for(int i=0;i<N;i++) for(int j=0;j<NBN;j++) copy(aji[i][0][j], aij[nbl[i][j]][0][jofi[i][j]]);
    for(int i=0;i<N;i++) for(int t=1;t<99;t++) copy(sit, si[0]);
    for(int pid=0;pid<NPROCS;pid++) {int num=0; int px=pid/NPY, py=pid%NPY;
        for(int i=0;i<N;i++) {int ix=(i/NLX)%NPX,iy=i/(NX*NLY); if(ix==px&&iy==py) atomid[pid][num++]=i;}
        for(int n=0;n<NATOMS;n++){int i=atomid[pid][n];
            for(int j=0;j<NBN;j++){int nbj=nbl[i][j],k;
                for(k=0;k<NATOMS;k++){if(atomid[pid][k]==nbj) {nbj_idx[i][j]=k;break;}}
                if(k==NATOMS) {nbj_idx[i][j]=num; atomid[pid][num++]=nbj;}}}
    }
}
void comp_adjt(double de, Vec aij, Vec vi, Vec vj, Vec res) {
    if(!(aij[X]||aij[Y])) zero(res);
    Vec dvij; Sub(vi, vj, dvij);
    double BB = Dot(aij, dvij), AA = Dot(aij, aij), bac = BB*BB-4*AA*de;
    if(bac<0) bac=0;
    double t = (-BB + SIGN(BB)*sqrt(bac))/(2*AA);
    Mul(aij, t, res);
}
void frog(int i, int t, int t0){ Vec tmp;
    Mul(si[t][V],T,tmp); Add(tmp,si[t][S],si[t+1][S]); Mul(si[t][A],T*T/2,tmp); Add(tmp,si[t+1][S],si[t+1][S]); SLoop(si[t+1][S]);
    for(int j=0;j<NBN;j++) new_Aij(i,j, t+1,t0);
    zero(si[t+1][A]); for(int j=0;j<NBN;j++) Add(aij[i][t+1][j], si[t+1][A], si[t+1][A]);
    Add(sit[A], si[t+1][A], si[t+1][V]); Mul(si[t+1][V], T/2, si[t+1][V]); Add(sit[V], si[t+1][V], si[t+1][V]);
    copy(a0[i][t+1], si[t+1][A]); copy(v0[i][t+1], si[t+1][V]);
}
State bufs[9][NATOMS][2];
int main() {
    MPI_Init(NULL, NULL);
    int my_pid, nprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_pid);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    int ix = my_pid % NPX, iy = my_pid / NPX;
    init(); copy(b, s);double E0=Esum(0); if(my_pid==0)printf("pid=%d, E_origin=%.16f\n", my_pid, E0);
    double e=0, esum=0, m=0, de=0, desum=0, msum=0;
    int t0=0, t1; double start_time = MPI_Wtime();
    for (int npr = 0; npr < 100; npr++) {
        t1=t0+1;
        for(int n=0;n<NATOMS;n++) frog(atomid[my_pid][n], t0, t0);
        for(int n=0;n<NATOMS;n++) copy(bufs[4][n][0], s[atomid[my_pid][n]][t0]);
        for(int n=0;n<NATOMS;n++) copy(bufs[4][n][1], s[atomid[my_pid][n]][t1]);

        MPI_Request reqs[4];
        int lpid = iy * NPX + XLOOP(ix - 1);
        int rpid = iy * NPX + XLOOP(ix + 1);

        MPI_Isend(&bufs[4], 2*NATOMS*sizeof(State), MPI_CHAR, lpid, TAG, MPI_COMM_WORLD, &reqs[0]);
        MPI_Isend(&bufs[4], 2*NATOMS*sizeof(State), MPI_CHAR, rpid, TAG, MPI_COMM_WORLD, &reqs[1]);
        MPI_Recv( &bufs[3], 2*NATOMS*sizeof(State), MPI_CHAR, lpid, TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv( &bufs[5], 2*NATOMS*sizeof(State), MPI_CHAR, rpid, TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        int upid = YLOOP(iy - 1) * NPX + ix;
        int dpid = YLOOP(iy + 1) * NPX + ix;
        MPI_Isend(&bufs[3], 3*2*NATOMS*sizeof(State), MPI_CHAR, upid, TAG, MPI_COMM_WORLD, &reqs[2]);
        MPI_Isend(&bufs[3], 3*2*NATOMS*sizeof(State), MPI_CHAR, dpid, TAG, MPI_COMM_WORLD, &reqs[3]);
        MPI_Recv( &bufs[0], 3*2*NATOMS*sizeof(State), MPI_CHAR, upid, TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv( &bufs[6], 3*2*NATOMS*sizeof(State), MPI_CHAR, dpid, TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Waitall(4, reqs, MPI_STATUS_IGNORE);

        for(int x=-1;x<=1;x++)for(int y=-1;y<=1;y++){
            int nb_pid = YLOOP(iy+y)*NPX + XLOOP(ix+x);
            int bufoff = 4+y*3+x;
            for(int n=0;n<NATOMS;n++) copy(b[atomid[nb_pid][n]][t0], bufs[bufoff][n][0]);
            for(int n=0;n<NATOMS;n++) copy(b[atomid[nb_pid][n]][t1], bufs[bufoff][n][1]);
        }
        for(int n=0;n<NATOMS;n++) {
            int i=atomid[my_pid][n];
            Vec dvit0;  Mul(si[t0][DV], 0.5, dvit0); Sub(si[t0][V], dvit0, dvit0);
            for(int j=0;j<NBN;j++) {
                int nbj = nbl[i][j];
                 Vec dvjt0;  Mul(b[nbj][t0][DV], 0.5, dvjt0); Sub(b[nbj][t0][V], dvjt0, dvjt0);
                 deij[i][t0][j] += Dot(dvit0, dvij[i][t0][j]) + Dot(dvjt0, dvji[i][t0][j]);
            }
            zero(si[t1][DV]); zero(si[t1][A]);
            Vec vit05; Add(si[t0][V], si[t1][V], vit05);
            for(int j=0;j<NBN;j++) { int nbj = nbl[i][j];
                new_Aji(i, j, t1, t0);
                Vec ait05; Add(aij[i][t0][j], aij[i][t1][j], ait05);
                Vec ajt05; Add(aji[i][t0][j], aji[i][t1][j], ajt05);
                Vec vjt05; Add(b[nbj][t0][V], b[nbj][t1][V], vjt05);
                double  dev = (Dot(ait05, vit05) + Dot(ajt05, vjt05))*T/4,
                        des = Ed(Dist_ij(i,nbj,t1)) - Ed(Dist_ij(i,nbj,t0));
                deij[i][t1][j] = deij[i][t0][j] + dev + des;   //comp deij
                comp_adjt(deij[i][t1][j], aij[i][t1][j], si[t1][V], b[nbl[i][j]][t1][V], ndv[i][t1][j]);
                // ndv 能量修正
                
                Vec at0; Add(aij[i][t0][j], aji[i][t0][j], at0);
                Vec at1; Add(aij[i][t1][j], aji[i][t1][j], at1);
                Add(at0, at1, mdv[i][t1][j]); Mul(mdv[i][t1][j], -T/4, mdv[i][t1][j]);
                // mdv 动量修正
                Balance_A(aij[i][t1][j], aji[i][t1][j]);
                Add(aij[i][t1][j],si[t1][A],si[t1][A]);
                Add(mdv[i][t1][j], ndv[i][t1][j], dvij[i][t1][j]);
                Sub(mdv[i][t1][j], ndv[i][t1][j], dvji[i][t1][j]);
                Add(dvij[i][t1][j], si[t1][DV], si[t1][DV]); //comp ndv
            }
            Add(si[t1][DV], si[t1][V], si[t1][V]);
        }
        e=0; esum=0; for(int n=0;n<NATOMS;n++) {int i=atomid[my_pid][n]; e += Edi(i,t0)+Evi(i,t0);}
        MPI_Reduce(&e, &esum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        de=0; desum=0; for(int n=0;n<NATOMS;n++) for(int j=0;j<NBN;j++) de += deij[atomid[my_pid][n]][t0][j];
        MPI_Reduce(&de, &desum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        m=0; msum=0; for(int n=0;n<NATOMS;n++) {int i=atomid[my_pid][n]; for(int j=0;j<NBN;j++) m += si[t0][V][X]+si[t0][V][Y];}
        MPI_Reduce(&m, &msum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if(my_pid==0) printf("E_%d=%.16f dE=%.2e, M=%.2e\n", t1, esum, desum*0.5+E0-esum, msum);
        t0=t1;
    }
    double exec_time = MPI_Wtime() - start_time;
    if (my_pid == MASTER) {
        printf("BSP Iters: %d\n", NITERS);
        printf("Time Sum: %.4lfs\n", exec_time);
        printf("Speed: %.4lf iters/s\n", NITERS/exec_time);
    }
    MPI_Finalize();
    return 0;
}
