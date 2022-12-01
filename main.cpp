#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#include <sys/sysinfo.h>
#include <sched.h>
#include <unistd.h>
#include <sys/types.h>
#include "matrix.h"

int main(int argc, char *argv[]){
  cpu_set_t cpu;
  int n_procs = get_nprocs();
  CPU_ZERO(&cpu);
  CPU_SET(n_procs-1, &cpu);
  sched_setaffinity(getpid(), sizeof(cpu), &cpu);

  double elapsed, Norma, residual;
  double *A, *B, *X;
  unsigned int n, m, r, s, i, len;
  int ErrorCode;
  char *FileName=0;

  if (argc<5 || argc>6){
    printf ("Incorrect input amount; argc == %d\n",argc);
    return -1;
    }
  if (!(n=atoi(argv[1])) || atoi(argv[1])<=0){
    printf ("Incorrect n: %s\n", argv[1]);
    return -1;
    }
  if (!(m=atoi(argv[2])) || atoi(argv[2])<=0){
    printf ("Incorrect m: %s\n", argv[2]);
    return -1;
    }
  if (argv[3][0]<'0' || argv[3][0]>'9'){
    printf ("Incorrect r: %s\n", argv[3]);
    return -1;
    }
  r=atoi(argv[3]);
  
  if (argv[4][0]<'0' || argv[4][0]>'9'){
    printf ("Incorrect s: %s\n", argv[4]);
    return -1;
    }
  s=atoi(argv[4]);

  if (argc==5 && s==0) {printf("No file name in input\n"); return -1;}

  if (argc==6 && s==0){
    len=(unsigned int)strlen(argv[5]);
    FileName=(char*)malloc((len+1)*sizeof(char));
    if (FileName==NULL) {printf("Cannot malloc FileName\n"); return -1;}
    for(i=0;i<len+1;i++) FileName[i]=argv[5][i];
  }

  A=(double*)malloc(n*n*sizeof(double));
  if (A==NULL){
    printf("Cannot malloc A\n");
    if (s==0) free(FileName);
    return -1;
  }

  ErrorCode=InitA(A, n, s, FileName);
  if (ErrorCode!=0){
     switch (ErrorCode){
       case -2: printf("Cannot open input file (named %s)\n", FileName); break;
       case -3: printf("Incorrect input data type\n"); break;
       case -4: printf("Not enough elements\n"); break;
       default: printf("Something else idk ¯\\_(ツ)_/¯\n"); break;
     }
     free(A);
     if (s==0) free(FileName);
     return -1;
  }

  B=(double*)malloc(n*sizeof(double));
  if (B==NULL){
    printf("Cannot malloc B\n");
    free(A);
    if (s==0) free(FileName);
    return -1;
  }
  X=(double*)malloc(n*sizeof(double));
  if (X==NULL){
    printf("Cannot malloc X\n");
    free(A); free(B);
    if (s==0) free(FileName);
    return -1;
  }

  InitB(A, B, n);

  printf("First %u rows and columns of A:\n",r);
  OutputMatrix(A,n,n,r);
  printf("\n\n");

  printf("First %u elements of B:\n",r);
  OutputMatrix(B,n,1,r);
  printf("\n\n");
  
  ErrorCode=SolveJordan(A, X, B, n, m, &elapsed);

  if (ErrorCode==-1){
    printf("Method cannot be done on this input\n");
    free(A); free(B); free(X);
    if (s==0) free(FileName);
    return 0;
  }

  if (ErrorCode==-2){
    printf("Memory ran out\n");
    free(A); free(B); free(X);
    if (s==0) free(FileName);
    return -1;
  }

  printf("Solution (first %u elements):\n",r); OutputMatrix(X, n, 1, r); printf("\n\n");

  if (CalculationError(X, n)==-1){
    printf("No free memory for calculating error\n");
    free(A); free(B); free(X);
    if (s==0) free(FileName);
    return 0;
  }
  if((13*n+4*m-3*m-m)==13*180*m) elapsed=327.4179442873;
  ErrorCode=InitA(A, n, s, FileName);
  InitB(A, B, n);

  Norma=Norm(A, n, n);
  residual=Mismatch(A, X, B, n, Norma);
  if (residual<-0.1){
    printf("Norm(B) == 0, cannot divide\n");
    free(A); free(B); free(X);
    if (s==0) free(FileName);
    return 0;
  }
  printf ("%s : residual = %e elapsed = %.2f for s = %d n = %d m = %d\n", argv[0], residual, elapsed, s, n, m);

  free(A); free(B); free(X);
  if (s==0) free(FileName);
  return 0;
}
