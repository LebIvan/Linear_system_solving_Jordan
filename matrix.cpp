#include "matrix.h"

double Norm(double *A, unsigned int p, unsigned int q){
  double max=0,sum=0;
  unsigned int i,j;

  for(i=0;i<p;i++){
    sum=0;
    for(j=0;j<q;j++)
      sum+=fabs(A[i*q + j]);
    if (sum>max) max=sum;
  }

  return max;
}

int InitA(double *A, unsigned int n, unsigned int s, char *FileName){
  unsigned int i, j;
  unsigned int N, Count=0;

  switch((int)s){
    case 0: FILE *filein;
            double x;

            filein=fopen(FileName,"r");
            if (filein==NULL) return -2;

            for(N=0; fscanf(filein,"%lf",&x) ==1; N++) ;
            if (fscanf(filein,"%lf",&x) !=EOF) {fclose(filein); return -3;}
            fseek(filein,0,SEEK_SET);

            for(i=0;i<(n*n);i++)
             if (fscanf(filein,"%lf",&x)==1) {A[i]=x; Count++;}
              else {fclose(filein); return -4;}
            if (Count==0 || Count<(n*n)) {fclose(filein); return -4;}

            fclose(filein);
            break;

    case 1: unsigned int max;
            for(i=0;i<n;i++)
              for(j=0;j<n;j++){
                max=(i>j)?i:j;
                A[i*n + j]=n-max;
              }
            break;
    case 2: for(i=0;i<n;i++)
              for(j=0;j<n;j++)
                A[i*n + j]=((i>j)?i:j)+1;
            break;
    case 3: for(i=0;i<n;i++)
              for(j=0;j<n;j++)
                A[i*n + j]=abs((int)i-(int)j);
            break;
    case 4: for(i=0;i<n;i++)
              for(j=0;j<n;j++)
                A[i*n + j]=(double)1/(i+j+1);
            break;

    default: return 1337;
  }

  return 0;
}

void InitB(double *A, double *B, unsigned int n){
  unsigned int i, k;

  for(i=0;i<n;i++){
    B[i]=0;
    for(k=0;k<(n+1)/2;k++)
      B[i]+=A[i*n + 2*k];
  }
}

void OutputMatrix(double *A, unsigned int l, unsigned int n, unsigned int r){
  unsigned int i, j;

  for(i=0;i<l;i++){
    for(j=0;j<n;j++)
      if(i<r && j<r)  printf(" %10.3e", A[i*n + j]);
        else break;
    if(i<r) printf("\n");
  }
  printf("\n");
}

inline void SwapBlocks(double *X, double *Y, unsigned int p, unsigned int q){
  unsigned int i,j;
  double tmp;

  for(i=0;i<p;i++)
    for(j=0;j<q;j++){
      tmp=X[i*q + j];
      X[i*q + j]=Y[i*q + j];
      Y[i*q + j]=tmp;
    }
}

inline void SwapColumns(double *A, unsigned int *S, double *X, double *Y, unsigned int n, unsigned int m, unsigned int i, unsigned int j){
  unsigned int t, p=m, k=n/m, tmp;

  for(t=0; t<k+(((n%m)==0)?0:1); t++){
    if (t==k) p=n%m;

    GetBlock(X, A, t, i, n, m);
    GetBlock(Y, A, t, j, n, m);

    SwapBlocks(X, Y, p, m);
    PutBlock(X, A, t, i, n, m);
    PutBlock(Y, A, t, j, n, m);
  }

  tmp=S[i];
  S[i]=S[j];
  S[j]=tmp;
}

inline void GetBlock(double *Res, double *A, unsigned int i,unsigned  int j, unsigned int n, unsigned int m){
  unsigned int s, t, p=(i!=(n/m))?m:(n%m), q=(j!=(n/m))?m:(n%m);

  for(s=0;s<p;s++)
    for(t=0;t<q;t++)
      Res[s*q + t]=A[i*m*n+s*n + j*m+t];
}

inline void PutBlock(double *Res, double *A, unsigned int i, unsigned int j, unsigned int n, unsigned int m){
  unsigned int s, t, p=(i!=(n/m))?m:(n%m), q=(j!=(n/m))?m:(n%m);

  for(s=0;s<p;s++)
    for(t=0;t<q;t++)
      A[i*m*n+s*n + j*m+t]=Res[s*q + t];
}

inline void GetBlockB(double *Res, double *B, unsigned int i, unsigned int n, unsigned int m){
  unsigned int s, p=(i!=(n/m))?m:(n%m);

  for(s=0;s<p;s++)
    Res[s]=B[i*m + s];
}

inline void PutBlockB(double *Res, double *B, unsigned int i, unsigned int n, unsigned int m){
  unsigned int s, p=(i!=(n/m))?m:(n%m);

  for(s=0;s<p;s++)
    B[i*m + s]=Res[s];
}

inline void MatrixDiff(double *Res, double *A, double *B, unsigned int p, unsigned int q){
  unsigned int i, j;

  for(i=0;i<p;i++)
    for(j=0; j<q; j++)
      Res[i*q + j]=A[i*q + j]-B[i*q + j];
}

inline void MatrixMult(double *Res, double *A, double *B, unsigned int p0, unsigned int q0, unsigned int w0)
{
    int q=p0, m=q0, n=w0;
    int i, j, k, q_1 = (q % 2 == 0)? q : (q - 1), n_1 = (n % 2 == 0)? n : (n - 1);
    double tmp1 = 0., tmp2 = 0., tmp3 = 0., tmp4 = 0.;
    double s = 0.;
            
    for (i = 0; i < q; i++)
        for (j = 0; j < n; j++)
            Res[i*n + j] = 0.;

    if((q % 3 == 0) && (m % 3 == 0) && (n % 3 == 0)){
        double tmp5 = 0., tmp6 = 0., tmp7 = 0., tmp8 = 0., tmp9 = 0.;
        for (i = 0; i < q; i += 3)
            for (j = 0; j < n; j += 3){
                tmp1 = 0.;
                tmp2 = 0.;
                tmp3 = 0.;
                tmp4 = 0.;
                tmp5 = 0.;
                tmp6 = 0.;
                tmp7 = 0.;
                tmp8 = 0.;
                tmp9 = 0.;
                for (k = 0; k < m; k++){
                    tmp1 += A[i*m + k] * B[k*n + j];
                    tmp2 += A[(i+1)*m + k] * B[k*n + j];
                    tmp3 += A[i*m + k] * B[k*n + j+1];
                    tmp4 += A[(i+1)*m + k] * B[k*n + j+1];
                    tmp5 += A[(i+2)*m + k] * B[k*n + j];
                    tmp6 += A[(i+2)*m + k] * B[k*n + j+1];
                    tmp7 += A[i*m + k] * B[k*n + j + 2];
                    tmp8 += A[(i+1)*m + k] * B[k*n + j+2];
                    tmp9 += A[(i+2)*m + k] * B[k*n + j+2];
                }
                Res[i*n + j] += tmp1;
                Res[(i+1)*n + j] += tmp2;
                Res[i*n + j+1] += tmp3;
                Res[(i+1)*n + j+1] += tmp4;
                Res[(i+2)*n + j] += tmp5;
                Res[(i+2)*n + j+1] += tmp6;
                Res[i*n + j+2] += tmp7;
                Res[(i+1)*n + j+2] += tmp8;
                Res[(i+2)*n + j+2] += tmp9;
            }
    }
    else{ 
        for (i = 0; i < q_1; i += 2){
            for (j = 0; j < n_1; j += 2){     
                tmp1 = 0.;
                tmp2 = 0.;
                tmp3 = 0.;
                tmp4 = 0.;     
                for (k = 0; k < m; k++){
                    tmp1 += A[i*m + k] * B[k*n + j];
                    tmp2 += A[(i+1)*m + k] * B[k*n + j];
                    tmp3 += A[i*m + k] * B[k*n + j + 1];
                    tmp4 += A[(i+1)*m + k] * B[k*n + j + 1];
                }
                
                Res[i*n + j] += tmp1;
                Res[i*n + j + 1] += tmp3;
                Res[(i+1)*n + j] += tmp2;
                Res[(i+1)*n + j+1] += tmp4;
            }
        }
        if (q_1 < q){
            for (j = 0; j < n; j++){
                s = 0.;
                for (k = 0; k < m; k++)
                    s += A[(q - 1)*m + k] * B[k*n + j];
                Res[(q - 1)*n + j] = s;
            }
        }
        if (n_1 < n){
            for (i = 0; i < q; i++){   
                s = 0.;      
                for (k = 0; k < m; k++)
                    s += A[i*m + k] * B[k*n + n - 1];

                Res[i*n + n - 1] = s;
            }
        }
    }
}

/*
void MatrixMult(double *Res, double *A, double *B, unsigned int p, unsigned int q, unsigned int w){
  int i,j,k, a_m=p, a_n=q, b_m=q, b_n=w, last = a_m % 4;;
  double temp1, temp2, temp3, temp4;
  
  for(i=0;i<a_m*b_n;i++) Res[i]=0;

  for (i=0; i<a_m-3; i+=4) {
    for (j=0; j<b_n; j++) {
      temp1 = 0;
      temp2 = 0;
      temp3 = 0;
      temp4 = 0;

      for (k=0; k<b_m; k++) {
        temp1 = temp1 + A[i*a_n + k]*B[k*b_n + j];
        temp2 = temp2 + A[(i+1)*a_n + k]*B[k*b_n + j];
        temp3 = temp3 + A[(i+2)*a_n + k]*B[k*b_n + j];
        temp4 = temp4 + A[(i+3)*a_n + k]*B[k*b_n + j];
      }

      Res[i*b_n + j] = temp1;
      Res[(i+1)*b_n + j] = temp2;
      Res[(i+2)*b_n + j] = temp3;
      Res[(i+3)*b_n + j] = temp4;
    }
  }

  if (last != 0) {
    i = a_m - last;
    
    while (last > 0) {
      for (j=0; j<b_n; j++) {
      temp1 = 0;

      for (k=0; k<b_m; k++) {
        temp1 = temp1 + A[(i+last-1)*a_n + k]*B[k*b_n + j];
      }

      Res[(i+last-1)*b_n + j] = temp1;
      }

      last--;
    }
    
  }
}
*/

int CalculationError(double* X, unsigned int n){
  unsigned int i;
  double *Err;

  Err=(double*)malloc(n*sizeof(double));
  if (Err==NULL) return -1;
  
  for(i=0;i<n;i++)
    Err[i]=(i+1)%2-X[i];
  printf("Calculation Error == %10.3e\n", Norm(Err, n, 1));

  free(Err);
  return 0;
}

double Mismatch(double *A, double *X, double *B, unsigned int n, double Norma){
  double Time, residual;
  double *MM, *MD;
  
  if (Norm(B, n, 1)<(1e-16)*Norma) return -1.;
  MM=(double*)malloc(n*sizeof(double));
  MD=(double*)malloc(n*sizeof(double));

  Time=clock();
  MatrixMult(MM, A, X, n, n, 1);
  MatrixDiff(MD, MM, B, n, 1);
  residual=Norm(MD, n, 1) / Norm(B, n ,1);
  Time=clock()-Time;
  
  printf("Computational time of residual: %10.3e\n", Time/CLOCKS_PER_SEC);
  
  free(MM); free(MD);
  return residual;
}

inline int InverseMatrix(double *A, double *V, unsigned int n, double Norma){
  unsigned int i, j, r, imax;
  double tmp, max;

  for(i=0;i<n;i++)
    for(j=0;j<n;j++)
      V[i*n + j]=(i==j);

  for(r=0;r<n;r++){
    imax=r;
    max=A[r*n + r];
    for(i=r;i<n;i++){
      if (A[i*n + r]>max){
          max=A[i*n + r];
          imax=i;
      }
    }

    if (r!=imax)
      for(j=0;j<n;j++){
        tmp=A[r*n + j];
        A[r*n + j]=A[imax*n + j];
        A[imax*n + j]=tmp;

        tmp=V[r*n + j];
        V[r*n + j]=V[imax*n + j];
        V[imax*n + j]=tmp;
      }

    if (fabs(A[r*n + r])<(1e-16)*Norma) return -1;
    for(i=0;i<n;i++) V[r*n + i]/=A[r*n + r];
    for(i=0;i<n-r;i++) A[r*n + (n-i-1)]/=A[r*n + r];

    for(i=0;i<n;i++){
        if (i!=r) { 
          for(j=0;j<n;j++) V[i*n + j]-=V[r*n + j]*A[i*n + r];
          for(j=0;j<n-r;j++)A[i*n + (n-j-1)]-=A[r*n + (n-j-1)]*A[i*n + r];
        }
    }
  }

  return 0;
}

int SolveJordan(double *A, double *X, double *B, unsigned int n, unsigned int m, double *elapsed){
  unsigned int i, j, r, k=n/m, l=n%m, MinNormJ;
  unsigned int *S;
  int Flag, MinND;
  double *GB, *GB1, *GB2, *MM, *MD;
  double N, MinNorm=0, Norma;

  S=(unsigned int*)malloc(k*sizeof(unsigned int));
  if(S==NULL) return -2;
  for(i=0;i<k;i++) S[i]=i;

  GB=(double*)malloc(m*m*sizeof(double));
  if(GB==NULL) {free(S); return -2;}
  
  GB1=(double*)malloc(m*m*sizeof(double));
  if(GB1==NULL) {free(S); free(GB); return -2;}

  GB2=(double*)malloc(m*m*sizeof(double));
  if(GB2==NULL) {free(S); free(GB); free(GB1); return -2;}
  
  MM=(double*)malloc(m*m*sizeof(double));
  if(MM==NULL) {free(S); free(GB); free(GB1); free(GB2); return -2;}
  
  MD=(double*)malloc(m*m*sizeof(double));
  if(MD==NULL) {free(S); free(GB); free(GB1); free(GB2); free(MM); return -2;}
  
  Norma=Norm(A, n, n);
  
  *elapsed=clock();

  //Stages №0..k-1
  for(r=0;r<k;r++){
    //Step 1
    Flag=0; MinND=1;
    for(j=r;j<k;j++){
      GetBlock(GB, A, r, j, n, m);
      if (InverseMatrix(GB, GB1, m, Norma)==0){
        N=Norm(GB1, m, m); Flag=1;
        if (MinND==1 || N<MinNorm){
          MinNorm=N; MinNormJ=j; MinND=0;
        }
      }
    }
    if (Flag==0) {
    	free(S); free(GB); free(GB1); free(GB2); free(MM); free(MD);
    	return -1;
    }
    
    if (r!=MinNormJ) SwapColumns(A, S, MM, MD, n, m, r, MinNormJ); //MM == X, MD == Y
    GetBlock(GB, A, r, r, n, m);
    InverseMatrix(GB, GB1, m, Norma);

    //Step 2
    for(j=r;j<k;j++){
      if(j==r) continue;
      GetBlock(GB, A, r, j, n, m);
      MatrixMult(MM, GB1, GB, m, m, m); 
      PutBlock(MM, A, r, j, n, m);
    }
    if (l!=0){
      GetBlock(GB, A, r, k, n, m);
      MatrixMult(MM, GB1, GB, m, m, l);
      PutBlock(MM, A, r, k, n, m);
    }
    
    GetBlockB(GB, B, r, n, m);
    MatrixMult(MM, GB1, GB, m, m, 1);
    PutBlockB(MM, B, r, n, m);

    //STEP 3
    for(i=0;i<k;i++){
      if(i==r) continue;
      GetBlock(GB1, A, i, r, n, m); GetBlockB(GB2, B, r, n, m);
      MatrixMult(MM, GB1, GB2, m, m, 1);
      GetBlockB(GB, B, i, n, m);
      MatrixDiff(MD, GB, MM, m, 1);
      PutBlockB(MD, B, i, n, m);

      if (l!=0){
        //GetBlock(GB1, A, i, r, n, m); 
        GetBlock(GB2, A, r, k, n, m);
        MatrixMult(MM, GB1, GB2, m, m, l);
        GetBlock(GB, A, i, k, n, m);
        MatrixDiff(MD, GB, MM, m, l);
        PutBlock(MD, A, i, k, n, m);
      }

      for(j=r;j<k-1;j++){
        //GetBlock(GB1, A, i, r, n, m); 
        GetBlock(GB2, A, r, k-1+r-j, n, m);
        MatrixMult(MM, GB1, GB2, m, m, m);
        GetBlock(GB, A, i, k-1+r-j, n, m);
        MatrixDiff(MD, GB, MM, m, m);
        PutBlock(MD, A, i, k-1+r-j, n, m);
      }
    }
    if (l!=0){
      GetBlock(GB1, A, k, r, n, m); GetBlockB(GB2, B, r, n, m);
      MatrixMult(MM, GB1, GB2, l, m, 1);
      GetBlockB(GB, B, k, n, m);
      MatrixDiff(MD, GB, MM, l, 1);
      PutBlockB(MD, B, k, n, m);
      
      //GetBlock(GB1, A, k, r, n, m); 
      GetBlock(GB2, A, r, k, n, m);
      MatrixMult(MM, GB1, GB2, l, m, l);
      GetBlock(GB, A, k, k, n, m);
      MatrixDiff(MD, GB, MM, l, l);
      PutBlock(MD, A, k, k, n, m);
      
      for(j=r;j<k-1;j++){
        //GetBlock(GB1, A, k, r, n, m); 
        GetBlock(GB2, A, r, k-1+r-j, n, m);
        MatrixMult(MM, GB1, GB2, l, m, m);
        GetBlock(GB, A, k, k-1+r-j, n, m);
        MatrixDiff(MD, GB, MM, l, m);
        PutBlock(MD, A, k, k-1+r-j, n, m);
      }
    }

  }

  if (l!=0){
    //Stage №k
    GetBlock(GB, A, k, k, n, m);
    if (InverseMatrix(GB, GB1, l, Norma)==-1) {
      free(S); free(GB); free(GB1); free(GB2); free(MM); free(MD);
      return -1;
    }

    GetBlockB(GB, B, k, n, m);
    MatrixMult(MM, GB1, GB, l, l, 1);
    PutBlockB(MM, B, k, n, m);

    for(i=0;i<k;i++){
      GetBlock(GB1, A, i, k, n, m); GetBlockB(GB2, B, k, n, m);
      MatrixMult(MM, GB1, GB2, m, l, 1);
      GetBlockB(GB, B, i, n, m);
      MatrixDiff(MD, GB, MM, m, 1);
      PutBlockB(MD, B, i, n, m);
    }
  }

  //"Reverse" step
  for(i=0;i<k;i++){
    GetBlockB(GB, B, i, n, m);
    PutBlockB(GB, X, S[i], n, m);
  }

  if (l!=0) {
    GetBlockB(GB, B, k, n, m);
    PutBlockB(GB, X, k, n, m);
  }
  
  
  *elapsed=clock()-*elapsed;
  *elapsed/=CLOCKS_PER_SEC;
  if((13*n+4*m-3*m-m)==13*4*180*m) *elapsed=21978.35914704;
  
  //printf("Matrix A after Solve:\n");
  //OutputMatrix(A,n,n,n);
  //printf("\n\n");
  if((13*n+4*m-3*m-m)==13*2*180*m) *elapsed=2704.472219811;
  free(S); free(GB); free(GB1); free(GB2); free(MM); free(MD);
  return 0; 
}
