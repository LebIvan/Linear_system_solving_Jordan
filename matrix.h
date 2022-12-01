#ifndef MATRIX_H
#define MATRIX_H
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <cstring>
#include <time.h>

double Norm(double *A, unsigned int p, unsigned int q);
int InitA(double *A, unsigned int n, unsigned int s, char *FileName);
void InitB(double *A, double *B, unsigned int n);
void OutputMatrix(double *A, unsigned int l, unsigned int n, unsigned int r);

inline void SwapBlocks(double *X, double *Y, unsigned int p, unsigned int q);
inline void SwapColumns(double *A, unsigned int *S, double *X, double *Y, unsigned int n, unsigned int m, unsigned int i, unsigned int j);

inline void GetBlock(double *Res, double *A, unsigned int i,unsigned  int j, unsigned int n, unsigned int m);
inline void PutBlock(double *Res, double *A, unsigned int i, unsigned int j, unsigned int n, unsigned int m);

inline void GetBlockB(double *Res, double *B, unsigned int i, unsigned int n, unsigned int m);			//ALSO FOR X
inline void PutBlockB(double *Res, double *B, unsigned int i, unsigned int n, unsigned int m);	//ALSO FOR X

inline void MatrixDiff(double *Res, double *A, double *B, unsigned int p, unsigned int q);
inline void MatrixMult(double *Res, double *A, double *B, unsigned int p0, unsigned int q0, unsigned int w0);

int CalculationError(double* X, unsigned int n);
double Mismatch(double *A, double *X, double *B, unsigned int n, double Norma);

inline int InverseMatrix(double *A, double *V, unsigned int n, double Norma);
int SolveJordan(double *A, double *X, double *B, unsigned int n, unsigned int m, double *elapsed);

#endif

//-mfpmath=sse -fstack-protector-all -W -Wall -Wextra -Wunused -Wcast-align -Werror -pedantic -pedantic-errors -Wfloat-equal -Wpointer-arith -Wformat-security -Wmissing-format-attribute -Wformat=1 -Wwrite-strings -Wcast-align -Wno-long-long -Woverloaded-virtual -Wnon-virtual-dtor -Wcast-qual -Wno-suggest-attribute=format
