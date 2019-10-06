#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <sys/time.h>

double get_sec()
{
    struct timeval time;
    gettimeofday(&time, NULL);
    return (time.tv_sec + 1e-6 * time.tv_usec);
}


void dgemm0(double *a, double *b, double *c, int n)
{
    double t0,t2;
    double optimized_time;
//    clock_t run_start,run_finish;
//    run_start = clock();
    t0 = get_sec();
    for(int i =0;i<n;i++){
        for (int j = 0; j < n ; j++){
            for(int k=0;k<n;k++){
                c[i*n+j] += a[i*n+k] * b[k*n+j];
            }
        }
    }
    t2 = get_sec();
    optimized_time = t2-t0;
    printf("%s - %4d: elapsed time is %8.5f second(s).\n", "dgemm0", n, optimized_time);
//    printf("\033[0;32m dgemm0 optimized kernel takes %8.5f seconds, performance is %5.2f GFLOPs. \033[0m\n", \
//            optimized_time, 2. * 1e-9 * n*n*n / optimized_time);
//    run_finish = clock();
//    double run_time = (double)(run_finish - run_start)/CLOCKS_PER_SEC*1000;
//    printf("dgemm0 run time is %fms\n",run_time);
}

void dgemm1(double *a, double *b, double *c, int n){

    double t0,t2;
    double optimized_time;
    t0 = get_sec();
    for (int i=0; i<n; i++)
        for (int j=0; j<n; j++) {
            register double r = c[i*n+j] ;
            for (int k=0; k<n; k++)
                r += a[i*n+k] * b[k*n+j];
            c[i*n+j] = r;
        }
    t2 = get_sec();
    optimized_time = t2-t0;
    printf("%s - %4d: elapsed time is %8.5f second(s).\n", "dgemm1", n, optimized_time);
//    printf("\033[0;32m dgemm1 optimized kernel takes %8.5f seconds, performance is %5.2f GFLOPs. \033[0m\n", \
//            optimized_time, 2. * 1e-9 * n*n*n / optimized_time);

}

void dgemm2(double *a, double *b, double *c, int n){
    int i,j,k;
    double t0,t2;
    double optimized_time;
    t0 = get_sec();
    for (i = 0; i < n; i+=2)
        for (j = 0; j < n; j+=2)
            for (k = 0; k < n; k+=2){
                c[i*n + j]         = a[i*n + k]*b[k*n + j] + a[i*n + k+1]*b[(k+1)*n + j]
                                     + c[i*n + j];
                c[(i+1)*n + j]     = a[(i+1)*n + k]*b[k*n + j] + a[(i+1)*n + k+1]*b[(k+1)*n + j]
                                     + c[(i+1)*n + j];
                c[i*n + (j+1)]     = a[i*n + k]*b[k*n + (j+1)] + a[i*n + k+1]*b[(k+1)*n + (j+1)]
                                     + c[i*n + (j+1)];
                c[(i+1)*n + (j+1)] = a[(i+1)*n + k]*b[k*n + (j+1)]
                                     + a[(i+1)*n + k+1]*b[(k+1)*n + (j+1)] + c[(i+1)*n + (j+1)];
            }
    t2 = get_sec();
    optimized_time = t2-t0;
    printf("%s - %4d: elapsed time is %8.5f second(s).\n", "dgemm2", n, optimized_time);

}

//double allocandRandomSample(double *a, double *b, double *c, int n){
//
//
//    a = (double *) calloc(sizeof(double), n);
//    b = (double *) calloc(sizeof(double), n);
//    c = (double *) calloc(sizeof(double), n);
//
//    srand(time(NULL));
//    for (int i = 0; i < n; i++) {
//        a[i] = (double)(rand() % 100 + 1);
//        b[i] = (double)(rand() % 100 + 1);
//    }
//    return *a, *b, *c;
//
//}

int main(int argc, char  *argv[])
{
    double *a,*b,*c;

    int nVal[] = {64,128,256,512,1024};
    int numAlgorithm = 3;
    int length = sizeof(nVal) / sizeof(nVal[0]);
    int fun;
    printf("\n*********** Register Reuse ***********\n");
    for( fun = 0; fun <numAlgorithm; fun++ ) {
        for (int i = 0; i < length; i++) {
            int n = nVal[i] * nVal[i];

//            allocandRandomSample(a,b,c,n);
//
            a = (double *) calloc(sizeof(double), n);
            b = (double *) calloc(sizeof(double), n);
            c = (double *) calloc(sizeof(double), n);

            srand(time(NULL));
            for (int i = 0; i < n; i++) {
                a[i] = rand() % 100 + 1;
                b[i] = rand() % 100 + 1;
            }
//            for (int i = 0; i < n; i++) {
//              printf("%f", a[i]);
//           }
            if(fun == 0)
                dgemm0(a, b, c, sqrt(n));
            else if(fun ==1)
                dgemm1(a, b, c, sqrt(n));
            else if(fun == 2)
                dgemm2(a, b, c, sqrt(n));
        }
    }
    return 0;
}
