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
    printf("\033[0;32m dgemm0 optimized kernel takes %8.5f seconds, performance is %5.2f GFLOPs. \033[0m\n", \
            optimized_time, 2. * 1e-9 * n*n*n / optimized_time);
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
    printf("\033[0;32m dgemm1 optimized kernel takes %8.5f seconds, performance is %5.2f GFLOPs. \033[0m\n", \
            optimized_time, 2. * 1e-9 * n*n*n / optimized_time);


}


int main(int argc, char  *argv[])
{
    double *a,*b,*c;

    int nVal[] = {64,128,256,512};
    int length = sizeof(nVal) / sizeof(nVal[0]);

    for(int i =0;i<length;i++) {
        printf("\033[0;31m======Begin with %d======\n", nVal[i]);
        printf("\n");
        int n = nVal[i] * nVal[i];
        a = (double *) calloc(sizeof(double), n);
        b = (double *) calloc(sizeof(double), n);
        c = (double *) calloc(sizeof(double), n);


        srand(time(NULL));
        for (int i = 0; i < n; i++) {
            a[i] = rand() % 100 + 1;
            b[i] = rand() % 100 + 1;
        }

        dgemm0(a, b, c, sqrt(n));
        dgemm1(a, b, c, sqrt(n));
        printf("\n");

//        printf("======================\n")
    }

    return 0;
}
