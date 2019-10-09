#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <sys/time.h>

void calLastRowandColumn(double *a,double *b,double *c,int n);
double get_sec();




void dgemm0(double *a, double *b, double *c, int n)
{
    int i,j,k;
    double t0,t2;
    double optimized_time;
//    clock_t run_start,run_finish;
//    run_start = clock();
    t0 = get_sec();
    for(i =0;i<n;i++){
        for (j = 0; j < n ; j++){
            for(k=0;k<n;k++){
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
    int i,j,k;
    double t0,t2;
    double optimized_time;
    t0 = get_sec();
    for (i=0; i<n; i++)
        for (j=0; j<n; j++) {
            register double r = c[i*n+j] ;
            for (k=0; k<n; k++)
                r += a[i*n+k] * b[k*n+j];
            c[i*n+j] = r;
        }
    t2 = get_sec();
    optimized_time = t2-t0;
    printf("%s - %4d: elapsed time is %8.5f second(s).\n", "dgemm1", n, optimized_time);

//    printf("\033[0;32m dgemm1 optimized kernel takes %8.5f seconds, performance is %5.2f GFLOPs. \033[0m\n", \
//            optimized_time, 2. * 1e-9 * n*n*n / optimized_time);

}

void dgemm2(double *a, double *b, double *c2, int n){
    int i,j,k;
    double t0,t2;
    double optimized_time;
    t0 = get_sec();
    for (i = 0; i < n; i+=2) {
        for (j = 0; j < n; j += 2) {
            register int t = i*n+j;
            register int tt = t+n;
            register double c00 = c2[t];
            register double c01 = c2[t+1];
            register double c10 = c2[tt];
            register double c11 = c2[tt+1];

            for (k = 0; k < n; k += 2) {
                /* 2 by 2 mini matrix multiplication using registers*/
                register int ta = i * n + k;
                register int tta = ta + n;

                register int tb = k * n + j;
                register int ttb = tb + n;

                register double a00 = a[ta];
                register double a01 = a[ta + 1];
                register double a10 = a[tta];
                register double a11 = a[tta + 1];

                register double b00 = b[tb];
                register double b01 = b[tb + 1];
                register double b10 = b[ttb];
                register double b11 = b[ttb + 1];

                c00 += a00 * b00 + a01 * b10;
                c01 += a00 * b01 + a01 * b11;
                c10 += a10 * b00 + a11 * b10;
                c11 += a10 * b01 + a11 * b11;
            }

            c2[t] = c00;
            c2[t+1] = c01;
            c2[tt] = c10;
            c2[tt+1] = c11;
        }
    }
    t2 = get_sec();
    optimized_time = t2-t0;
    printf("%s - %4d: elapsed time is %8.5f second(s).\n", "dgemm2", n, optimized_time);


}

void dgemm3(double *a, double *b, double *c, int n){
    int i,j,k;
    int boundary = n-(n%3);
    double t0,t2;
    double optimized_time;
    t0 = get_sec();
    for (i = 0; i <boundary;i+=3) {
        for(j=0; j<boundary;j+=3){
            register int t = i*n+j; //row 0
            register int tt = t+n;  //row 1
            register int ttt = tt+n; //row 2
             //define 9 register for c result.
            register double c00 = c[t];
            register double c01 = c[t+1];
            register double c02 = c[t+2];

            register double c10 = c[tt];
            register double c11 = c[tt+1];
            register double c12 = c[tt+2];

            register double c20 = c[ttt];
            register double c21 = c[ttt+1];
            register double c22 = c[ttt+2];

            register int ta ;
            register int tta ;
            register int ttta ;

            register int tb;
            register int ttb;
            register int tttb;
            //define 6 register for a and b.
            register double a00;
            register double a10;
            register double a20;

            register double b00;
            register double b01;
            register double b02;

            for(k = 0; k<boundary; k+=3){
                ta = i*n+k;  // a row 0
                tta = ta+n;  // a row 1
                ttta = tta+n; // a row 2

                tb = k*n+j;   // b column 0
                ttb = tb +n;   // b column 1
                tttb = ttb +n;   // b column 2

                //start: calculate value of one block in one inner loop
                a00 = a[ta];  //row1
                a10 = a[tta]; // row2
                a20 = a[ttta];//row3

                b00 = b[tb];  // column 1
                b01 = b[tb+1];// column 2
                b02 = b[tb+2];//column 3

                //first
                c00 += a00 * b00;
                c01 += a00 * b01;
                c02 += a00 * b02;

                c10 += a10 * b00;
                c11 += a10 * b01;
                c12 += a10 * b02;

                c20 += a20 * b00;
                c21 += a20 * b01;
                c22 += a20 * b02;

                //second
                a00 = a[ta+1];
                a10 = a[tta+1];
                a20 = a[ttta+1];

                b00 = b[ttb];
                b01 = b[ttb+1];
                b02 = b[ttb+2];

                c00 += a00 * b00;
                c01 += a00 * b01;
                c02 += a00 * b02;

                c10 += a10 * b00;
                c11 += a10 * b01;
                c12 += a10 * b02;

                c20 += a20 * b00;
                c21 += a20 * b01;
                c22 += a20 * b02;

                //third
                a00 = a[ta+2];
                a10 = a[tta+2];
                a20 = a[ttta+2];

                b00 = b[tttb];
                b01 = b[tttb+1];
                b02 = b[tttb+2];

                c00 += a00 * b00;
                c01 += a00 * b01;
                c02 += a00 * b02;

                c10 += a10 * b00;
                c11 += a10 * b01;
                c12 += a10 * b02;

                c20 += a20 * b00;
                c21 += a20 * b01;
                c22 += a20 * b02;

            }
            //start to calculate the vaule from boundary to n
            for(;k<n;k++){
                ta = i*n+k;  // a row 0
                tta = ta+n;  // a row 1
                ttta = tta+n; // a row 2

                tb = k*n+j;   // b column 0

                a00 = a[ta];  //row1
                a10 = a[tta]; // row2
                a20 = a[ttta];//row3

                b00 = b[tb];  // column 1
                b01 = b[tb+1];// column 2
                b02 = b[tb+2];//column 3

                c00 += a00 * b00;
                c01 += a00 * b01;
                c02 += a00 * b02;

                c10 += a10 * b00;
                c11 += a10 * b01;
                c12 += a10 * b02;

                c20 += a20 * b00;
                c21 += a20 * b01;
                c22 += a20 * b02;
            }
            //end

            c[t] = c00;
            c[t+1] = c01;
            c[t+2] = c02;

            c[tt] = c10;
            c[tt+1] = c11;
            c[tt+2] = c12;

            c[ttt] = c20;
            c[ttt+1] = c21;
            c[ttt+2] = c22;
            //end 
        }  
    }
    calLastRowandColumn(a,b,c,n);

    t2 = get_sec();
    optimized_time = t2-t0;
    printf("%s - %4d: elapsed time is %8.5f second(s).\n", "dgemm3", n, optimized_time);
}


int main(int argc, char  *argv[])
{
    double *a,*b,*c;

    int nVal[] = {64,128,256,512,1024,2048};
   // int nVal[] = {64};
    int numAlgorithm = 4;
    int length = sizeof(nVal) / sizeof(nVal[0]);
    int fun,i;
    printf("\n*********** Register Reuse ***********\n");
    for( fun = 0; fun <numAlgorithm; fun++ ) {
        for (i = 0; i < length; i++) {
            int n = nVal[i] * nVal[i];

//            allocandRandomSample(a,b,c,n);

            a = (double *) calloc(sizeof(double), n);
            b = (double *) calloc(sizeof(double), n);
            c = (double *) calloc(sizeof(double), n);

            srand(time(NULL));
            int j;
            for (j = 0; j < n; j++) {
                a[j] = rand() % 100 + 1;
                b[j] = rand() % 100 + 1;
                c[j] = 0.0;
            }

            if(fun == 0)
                dgemm0(a, b, c, sqrt(n));
            else if(fun ==1)
                dgemm1(a, b, c, sqrt(n));
            else if(fun == 2)
                dgemm2(a, b, c, sqrt(n));
            else if (fun == 3)
                dgemm3(a, b, c, sqrt(n));
        }
   }
    return 0;
}


double get_sec()
{
    struct timeval time;
    gettimeofday(&time, NULL);
    return (time.tv_sec + 1e-6 * time.tv_usec);
}

void calLastRowandColumn(double *a,double *b,double *c,int n){
    int remainAmount = n % 3;
    if(remainAmount > 0){
        int i,j,k;
        register double result;
        for(i = n - remainAmount;i<n;i++) {
            for (j = 0; j < n; j++) {
                result = c[i*n+j];
                for (k = 0; k < n; k++) {
                    result += a[i * n + k] * b[k * n + j];
                }
                c[i * n + j] = result;
            }
        }

        for(i = 0;i<n-remainAmount;i++){
            for(j = n-remainAmount; j<n;j++){
                result = c[i*n+j];
                for(k=0;k<n;k++){
                    result += a[i * n + k] * b[k * n + j];
                }
                c[i * n + j] = result;
            }
        }

    }
}