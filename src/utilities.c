#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <complex.h> 
#include <fftw3.h>
#include "utilities.h"


/* Unwrap phase data. */

void unwrap( int ll, double *phase ) {
  /* CHH edit removed unused variable cup */
  double test=PI;
  int i;

  for(i=1;i<ll;i++){
    if((phase[i-1]-phase[i])<-test){
      phase[i]=phase[i]-2.0*PI;
      i--;
    }
    else if((phase[i-1]-phase[i])>test){
      phase[i]=phase[i]+2.0*PI;
      i--;
    }
  }
}

/*         */
/* Filters */
/*         */

double hann(int begin,int center,int i){
  double answer;

  if((begin<=i) && (i<=(2*center-begin))){
    answer=0.5*(cos(PI*(i-center)/(center-begin))+1.0);
  }
  else answer=0.0;

  return answer;
}

double cont_half_hann(double begin,double end,double x){
  double answer,cup;

  if(x<begin) answer=1.0;
  else if((begin<=x) && (x<=end)){
    cup=PI*(x-begin)/(end-begin);
    answer=0.5*(cos(cup)+1.0);
  }
  else answer=0.0;

  return answer;
}

double half_hann(int begin,int end,int i){
  double answer;

  if(i<begin) answer=1.0;
  else if((begin<=i) && (i<=end)){
    answer=0.5*(cos(PI*((i-begin))/((int)(end-begin)))+1.0);
  }
  else answer=0.0;

  return answer;
}


/*  Choose N as 2^power.  */
int make_N(int power){
  int i,answer;

  answer=1;
  for(i=1;i<power+1;i++) answer=2*answer;

  return answer;
}

/*               */
/* Data analysis */
/*               */

double correllation(double *vec1,double *vec2,double *corr,\
		    int n_fft,double f_s,double t_start,double t_end,double t_window,\
		    double f_b,double roll_b,double f_t,double roll_t){

  double complex *cup, *fft_cup, *fft_cor;
  double tt,dt,df,pot;
  int i;
  fftw_plan p;

  dt=1.0/f_s;
  df=1.0/n_fft/dt;
  cup=(double complex *)malloc(sizeof(double complex)*n_fft);
  fft_cup=(double complex *)malloc(sizeof(double complex)*n_fft);
  fft_cor=(double complex *)malloc(sizeof(double complex)*n_fft);
  
  p=fftw_plan_dft_1d(n_fft,cup,fft_cup,FFTW_BACKWARD,FFTW_ESTIMATE);
  for(i=0 ;i*dt<t_start && i<n_fft;i++){
    cup[i]=0.0;
  }
  for(;i*dt<t_end && i<n_fft;i++){
    cup[i]=vec1[i];
  }
  for(;i<n_fft;i++){
    cup[i]=0.0;
  }
/*   for(i=0;i<n_fft;i++){ */
/*     cup[i]=vec1[i]; */
/*   } */
  fftw_execute(p);
  for(i=0;i<n_fft;i++){
    fft_cor[i]=fft_cup[i]*dt;
    cup[i]=vec2[i];
  }
  fftw_execute(p);
  fftw_destroy_plan(p);
  p=fftw_plan_dft_1d(n_fft,cup,fft_cup,FFTW_FORWARD,FFTW_ESTIMATE);

  for(i=0;i<n_fft;i++){
    fft_cor[i]=(fft_cor[i])*conj(fft_cup[i])*dt*(1.0-cont_half_hann(0.01,0.05,i*df));
    cup[i]=fft_cor[i];
  }

  int i1=(int)(f_b/df);
  int j1=(int)(f_b*(1.0-roll_b)/df);
  int i2=(int)(f_t/df);
  int j2=(int)(f_t*(1.0+roll_t)/df);
  for(i=0;i<n_fft;i++){
    cup[i]=dt*cup[i];
    cup[i]=(1.0-half_hann(j1,i1,i))*cup[i];
    cup[i]=(half_hann(i2,j2,i))*cup[i];
  }

  fftw_execute(p);
  fftw_destroy_plan(p);

  for(i=-n_fft/2;i<0;i++) corr[i+n_fft/2]=creal(fft_cup[i+n_fft])*df;
  for(i=0;i<n_fft/2;i++) corr[i+n_fft/2]=creal(fft_cup[i])*df;

  tt=-t_window;
  pot=0.0;  
  for(i=-(int)(t_window/dt);i<0;i++){
    if(creal(fft_cup[i+n_fft])>pot){
      pot=creal(fft_cup[i+n_fft]);
      tt=i*dt;
    }
  }
  for(;i<=(int)(t_window/dt);i++){
    if(creal(fft_cup[i])>pot){
      pot=creal(fft_cup[i]);
      tt=i*dt;
    }
  }

  free(cup);
  free(fft_cup);
  free(fft_cor);

  return tt;
}

/*               */
/* Root finders. */
/*               */


/*  A Newton-Rapheson scheme for finding zeros of analytic functions.  */
/* !!!!!!!!!!! Doesn't work !!!!!!!!!!!!! */
/* gnu complex doesn't seem to like passing  */
/* complex values functions to functions */

double complex cnewt_raph(double complex start_step, double complex guess,
		   double complex (*funct)(double complex )){
  double complex pot,last_guess,step_size,deriv,answer;

  step_size=start_step;
  pot=1.0;
  while(cabs(pot)>ROOT_TOL){
    last_guess=guess;
  printf("here 1\n");
    deriv=0.5*(FUNKT(guess+step_size)-FUNKT(guess-step_size))/step_size;
  printf("here 2\n");
    pot=FUNKT(guess);
    if(cabs(deriv)==0){
      step_size=I*step_size;
      deriv=0.5*(FUNKT(guess+step_size)-FUNKT(guess-step_size))/step_size;
      if(cabs(deriv)==0) {
	printf("Looks like Cnewt_raph may have found a minimum.\n");
	return guess;
      }
    }
    guess=guess-pot/deriv;
/*     printf("here .... %e   %e    %e\n",guess,cabs(pot)); */
    if(cabs(step_size) > 1.0e-19){
      step_size=0.2*(last_guess-guess);
    }
  }

  answer=guess;
  return answer;
}

/*                                                         */
/* The routines below are used to read data in from files. */
/*                                                         */

/*  read_header skips the gnuplot comment lines on the top of a file.  */

void read_header(FILE *f){
  int c;

  while((c=getc(f))!=EOF){
    if(c=='#'){
      while((c=getc(f))!=EOF && c!='\n'){
	;   /*  empty line */
      }
    }
    else if(c!='\n'){
      ungetc(c,f);
      break;
    }
  }
}


/* Count the rows in a file. */

int count_rows(char *filename,int num_col){
  int m,n_rows;
  double x,y;
  FILE *f=fopen(filename,"r");

  if(f==NULL){
    perror(filename);
    exit(2);
  }
  n_rows=0;
  read_header(f);
  while(fscanf(f,"%lf",&x)==1){
    for(m=0;m<num_col;m++){
      if(fscanf(f,"%lf",&y)!=1){
	fprintf(stderr,"%s: count_rows format error %d   %d\n",\
		filename,n_rows,m);
	exit(3);
      }
    }
    ++n_rows;
  }
  fclose(f);

  return n_rows;
}

int count_rows_arbcol(char *filename){
  int answer,c;
  FILE *f=fopen(filename,"r");

  if(f==NULL){
    perror(filename);
    exit(2);
  }
  answer=0;
  read_header(f);
  while((c=getc(f))!=EOF){
    if(c=='\n') answer=answer+1;
  }
  fclose(f);

  return answer;
}

/* Count columns */

int count_columns(char *filename){
  int answer,a,c;
  FILE *f=fopen(filename,"r");

  if(f==NULL){
    perror(filename);
    exit(2);
  }
  read_header(f);

  answer=0;
  a='b';

/*   while(a!='\n'){ */
/*     while(a!=' ' && a!='\t'){ */
/*       a=getc(f); */
/*       printf("%c\n",a); */
/*     } */
/*     answer=answer+1; */
/*     while(a==' ' || a=='\t') a=getc(f); */
/*   } */
  fclose(f);

  return answer;
}

/*  read_data reads the data in from "filename" and saves it in  */
/*  the vectors q_vec and vec real part, imag part. */
/*  returns the number of rows. */

int read_data_vec(char *filename,double *q_vec,double complex *vec){
  int n;
  double re_p,im_p;
  FILE *f=fopen(filename,"r");

  if(f==NULL){
    perror(filename);
    exit(2);
  }
  n=0; 
  read_header(f);
  while(fscanf(f,"%lf",&q_vec[n])==1){
    if(fscanf(f,"%lf%lf",&re_p,&im_p)!=2){
      fprintf(stderr,"%s: read_data format error\n",filename);
      exit(3);
    }
    vec[n]=re_p+I*im_p;
    ++n;
  }
  fclose(f);
/*   printf("Number of rows = %d\n",n); */
  return n;
}

int read_real_data(char *filename,double *t_vec,double *vec){
  int n;
  double p;
  FILE *f=fopen(filename,"r");

  if(f==NULL){
    perror(filename);
    exit(2);
  }
  n=0; 
  read_header(f);
/*   while(fscanf(f,"%.12lf",&t_vec[n])==1){ */
  while(fscanf(f,"%lf",&t_vec[n])==1){
    if(fscanf(f,"%lf",&p)!=1){
      fprintf(stderr,"%s: read_data format error  %d\n",filename,n);
      exit(3);
    }
    vec[n]=p;
    ++n;
  }
  fclose(f);
/*   printf("Number of rows = %d\n",n); */
  return n;
}

int read_1col_datafile(char *filename,double *vec){
  int n;
  double p;
  FILE *f=fopen(filename,"r");

  if(f==NULL){
    perror(filename);
    exit(2);
  }
  n=0; 
  read_header(f);
  while(fscanf(f,"%lf",&p)==1){
    vec[n]=p;
    ++n;
  }
  fclose(f);
/*   printf("Number of rows = %d\n",n); */
  return n;
}

int read_real_data_ncols(char *filename,int n_cols,double *t_vec,double **vec){
  int j,n;
  double p;
  FILE *f=fopen(filename,"r");

  if(f==NULL){
    perror(filename);
    exit(2);
  }
  n=0; 
  read_header(f);
  while(fscanf(f,"%lf",&t_vec[n])==1){
    for(j=0;j<n_cols;j++){
      if(fscanf(f,"%lf",&p)!=1){
	fprintf(stderr,"%s: read_data_ncols format error  line %d col %d\n",\
		filename,n,j);
	exit(3);
      }
      vec[j][n]=p;
    }
    ++n;
  }
  fclose(f);
/*   printf("Number of rows = %d\n",n); */
  return n;
}

int read_real_data_ncols2(char *filename,int n_cols,double **vec){
  int j,n;
  double p;
  FILE *f=fopen(filename,"r");

  if(f==NULL){
    perror(filename);
    exit(2);
  }
  n=0; 
  read_header(f);
  while(fscanf(f,"%lf",&p)==1){
    vec[0][n]=p;
    for(j=0;j<n_cols;j++){
      if(fscanf(f,"%lf",&p)!=1){
	fprintf(stderr,"%s: read_data_ncols2 format error  line %d col %d\n",\
		filename,n,j);
/* 	exit(3); */
      }
      vec[j+1][n]=p;
    }
    ++n;
  }
  fclose(f);
/*   printf("Number of rows = %d\n",n); */
  return n;
}

int read_real_data_ncols3(char *filename,int n_cols,double **vec){
  int j,n;
  double p;
  FILE *f=fopen(filename,"r");

  if(f==NULL){
    perror(filename);
    exit(2);
  }
  n=0; 
  read_header(f);
  while(fscanf(f,"%lf",&p)==1){
    for(j=0;j<n_cols;j++){
      if(fscanf(f,"%lf",&p)!=1){
	fprintf(stderr,"%s: read_data_ncols3 format error  line %d col %d\n",\
		filename,n,j);
/* 	exit(3); */
      }
      vec[j][n]=p;
    }
    ++n;
  }
  fclose(f);
/*   printf("Number of rows = %d\n",n); */
  return n;
}

int read_data_vec_ncols(char *filename,int n_cols,double *q_vec,
						double complex **vec){
  int j,n;
  double re_p,im_p;
  FILE *f=fopen(filename,"r");

  if(f==NULL){
    perror(filename);
    exit(2);
  }
  n=0; 
  read_header(f);
  while(fscanf(f,"%lf",&q_vec[n])==1){
    for(j=0;j<n_cols;j++){
      if(fscanf(f,"%lf%lf",&re_p,&im_p)!=2){
	fprintf(stderr,"%s: read_data format error\n",filename);
	exit(3);
      }
    vec[j][n]=re_p+I*im_p;
    }
    ++n;
  }
  fclose(f);
/*   printf("Number of rows = %d\n",n); */
  return n;
}

/*                                                      */
/* The routines below are used to read data into files. */
/*                                                      */

int save_points(char *filename,int n,int m,double **vec){
  int i,j;
  FILE *f;

  f=fopen(filename,"w");
  for(i=0;i<m-1;i++){
    for(j=0;j<n-1;j++){
      fprintf(f,"%e   ",vec[j][i]);
    }
    fprintf(f,"%e\n",vec[j][i]);
  }

  fclose(f);
  return(0);
}

int save_complex_points(char *filename,int n,int m,double complex **vec){
  int i,j;
  FILE *f;

  f=fopen(filename,"w");
  for(i=0;i<m-1;i++){
    for(j=0;j<n-1;j++) fprintf(f,"%e   %e   ",\
			       creal(vec[j][i]),cimag(vec[j][i]));
    fprintf(f,"%e   %e\n",creal(vec[j][i]),cimag(vec[j][i]));
  }

  fclose(f);
  return(0);
}

/*                              */
/* Memory allocation routines.  */
/*                              */

/*  Allocate memory for a complex matrix.  */

double complex **matrix_alloc(int n, int m)
{
  double complex **p;
  int j;
  
  p=(double complex **) malloc(n*sizeof(double complex *));
  for (j = 0; j < n; j++) {
    p[j] = (double complex *) malloc(m * sizeof(double complex));
  }  
  if(!p){
    printf("error allocating complex matrix memory\n");
    exit(1);
  }
  return p;
}

/*  Allocate memory for a real */

double **real_matrix_alloc(int n, int m)
{
  double **p;
  int j;
  
  p=(double **) malloc(n*sizeof(double *));
  for (j = 0; j < n; j++) {
    p[j] = (double *) malloc(m * sizeof(double));
  }  
  if(!p){
    printf("error allocating real matrix memory\n");
    exit(1);
  }
  return p;
}

void matrix_free(double complex **p, int n, int m)
{
  int j;
  for (j = 0; j < n; j++) free(p[j]);
  free(p);
}

void real_matrix_free(double **p, int n, int m)
{
  int j;
  for (j = 0; j < n; j++){
    free(p[j]);
  }
  free(p);
}

/*                */
/* Linear Algebra */
/*                */

/*  Solve Ax=y for a tridiagonal L by L matrix A. */
/*  "diag" is the diagonal of A, */
/*  "lower" is the lower off-diagonal and */
/*  "upper" is the upper off-diagonal. */
/*  The algorithm is an implicit LU decomposition. */

void tridiag(int L, double complex *lower, double complex *diag,
			 double complex *upper, double complex *y,
			 double complex *x){
  double complex *vector,cup,pot;
  int i; 
 
  vector=(double complex *)malloc(sizeof(double complex)*L);

  cup=diag[0];
  if(cabs(cup)==0.0){
    printf("First diagonal element is 0 in routine \"tridiag\", goober!\n");
    exit(1);
  }
  vector[0]=cup;
  x[0]=y[0];
  for(i=1;i<L;i++){
    pot=lower[i-1]/cup;
    x[i]=y[i]-(x[i-1]*pot);
    cup=diag[i]-(upper[i-1]*pot);
    if(cabs(cup)==0.0){
      printf("Need to pivot in routine \"tridiag\", goober!\n");
      exit(2);
    }
    vector[i]=cup;
  }

  x[L-1]=x[L-1]/vector[L-1];
  for(i=L-2;i>=0;i--){
    x[i]=(x[i]-(upper[i]*x[i+1]))/vector[i];
  }

  free(vector);
}

/* Does the LU decomposition for a matrix with real valued elements. */

void do_real_LU_decomp(int NN,double **LU,double **matrix){
  int j,k,m; 
  double cup;

  for(j=0;j<NN;j++) {
    for(k=j;k<NN;k++) {
      cup=matrix[j][k];
      for (m=0;m<j;m++) {
	cup=cup-(LU[j][m]*LU[m][k]); 
      }
      LU[j][k]=cup; 
    }
    if(fabs(LU[j][j])<1.0e-12){
      printf("\n......... Got to pivot Rog. j=%d. ...........\n\n",j);
    }
    for(k=j+1;k<NN;k++) {
      cup=matrix[k][j];
      for (m=0;m<j;m++) {
	cup=cup-(LU[k][m]*LU[m][j]); 
      }
      LU[k][j]=cup/LU[j][j];
    }
  }
}


/*  These routines uses the LU decomposition to solve the linear system  */
/*  LU bc_vec = source_vec. */

void real_LU_linsolver(int NN,double **LU,double *source_vec,double *bc_vec){
  int j,k;
  double pot,*temp; 

  temp=(double *)malloc(NN*sizeof(double));

  for(j=0;j<NN;j++){ 
    pot=source_vec[j];
    for(k=0;k<j;k++) pot=pot-(LU[j][k]*temp[k]); 
    temp[j]=pot;
  }
  for(j=NN-1;j>-1;j--){ 
    pot=temp[j];
    for(k=j+1;k<NN;k++) pot=pot-(LU[j][k]*bc_vec[k]); 
    bc_vec[j]=pot/LU[j][j];
  }

  free(temp);
}

double complex det_2by2(double complex **MM){
  double complex answer;

  answer=MM[0][0]*MM[1][1]-MM[0][1]*MM[1][0];
  return answer;
}

void inverse_2by2(double complex **invMM,double complex **MM){
  double complex invdet;

  if(det_2by2(MM)==0.0){
    printf("matrix singular\n");
    exit(0);
  }
  invdet=1.0/det_2by2(MM);
  invMM[0][0]=MM[1][1]*invdet;
  invMM[0][1]=-MM[0][1]*invdet;
  invMM[1][0]=-MM[1][0]*invdet;
  invMM[1][1]=MM[0][0]*invdet;
}

void eigs_2by2(double complex *eigs, double complex **eigvecs,
			   double complex **MM){
  double complex disc_sqrt,cup;

  disc_sqrt=csqrt((MM[0][0]-MM[1][1])*(MM[0][0]-MM[1][1])+4.0*MM[0][1]*MM[1][0]);
  eigs[0]=0.5*(MM[0][0]+MM[1][1]-disc_sqrt);
  eigs[1]=0.5*(MM[0][0]+MM[1][1]+disc_sqrt);

  if(eigs[0]==eigs[1]){
    printf("eigenvalue degeneracy eig=%f%+fi\n",creal(eigs[0]),cimag(eigs[0]));
    exit(0);
  }

  eigvecs[0][0]=-MM[0][1];
  eigvecs[0][1]=MM[0][0]-eigs[0]; 
  cup=csqrt(eigvecs[0][0]*conj(eigvecs[0][0])+eigvecs[0][1]*conj(eigvecs[0][1]));
  if(cup==0.0){
    eigvecs[0][0]=MM[1][1]-eigs[0];
    eigvecs[0][1]=-MM[1][0];
    cup=csqrt(eigvecs[0][0]*conj(eigvecs[0][0])+eigvecs[0][1]*conj(eigvecs[0][1]));
  }
  eigvecs[0][0]=eigvecs[0][0]/cup;
  eigvecs[0][1]=eigvecs[0][1]/cup;

  eigvecs[1][0]=-MM[0][1];
  eigvecs[1][1]=MM[0][0]-eigs[1];
  cup=csqrt(eigvecs[1][0]*conj(eigvecs[1][0])+eigvecs[1][1]*conj(eigvecs[1][1]));
  if(cup==0.0){
    eigvecs[1][0]=MM[1][1]-eigs[1];
    eigvecs[1][1]=-MM[1][0];
    cup=csqrt(eigvecs[1][0]*conj(eigvecs[1][0])+eigvecs[1][1]*conj(eigvecs[1][1]));
  }
  eigvecs[1][0]=eigvecs[1][0]/cup;
  eigvecs[1][1]=eigvecs[1][1]/cup;
}

/*                          */
/* Interpolation algorithms */
/*                          */

/* Cubic spline interpolation routines. */


/*  Cubic spline moments (second derivs at given points) */
/*  with boundary conditions: sec deriv 0 at endpoints. */
/*  MM is the number of segments. */
/*  MM+1 is the number of points given. */

void cub_spline_coef1(int MM,double complex *moments,double *x,double complex *y){
  int i;
  double complex cup,*lower,*diag,*upper,*source_vec;
  double pot;

  lower=(double complex *)malloc(sizeof(double complex)*(MM+1));
  diag=(double complex *)malloc(sizeof(double complex)*(MM+1));
  upper=(double complex *)malloc(sizeof(double complex)*(MM+1));
  source_vec=(double complex *)malloc(sizeof(double complex)*(MM+1));

  diag[0]=2.0;
  upper[0]=0.0;
  source_vec[0]=0.0;
  for(i=1;i<MM;i++){
    pot=x[i+1]-x[i-1];
    upper[i]=(x[i+1]-x[i])/pot;
    lower[i-1]=1.0-upper[i];
    diag[i]=2.0;
    source_vec[i]=6.0/pot;
    cup=(y[i-1]-y[i])/(x[i]-x[i-1]);
    cup=cup+(y[i+1]-y[i])/(x[i+1]-x[i]);
    source_vec[i]=source_vec[i]*cup;
  }
  diag[MM]=2.0;
  lower[MM-1]=0.0;
  source_vec[MM]=0.0;

  tridiag(MM+1,lower,diag,upper,source_vec,moments);

  free(lower);
  free(diag);
  free(upper);
  free(source_vec);
}

/*  Cubic spline moments (second derivs at given points) */
/*  with boundary conditions: deriv y0_p and yMM_p at endpoints. */
/*  MM is the number of segments. */
/*  MM+1 is the number of points given. */

void cub_spline_coef2(int MM,double complex y0_p,double complex yMM_p,
		      double complex *moments,double *x,double complex *y){
  int i;
  double complex cup,*lower,*diag,*upper,*source_vec;
  double pot;

  lower=(double complex *)malloc(sizeof(double complex)*(MM+1));
  diag=(double complex *)malloc(sizeof(double complex)*(MM+1));
  upper=(double complex *)malloc(sizeof(double complex)*(MM+1));
  source_vec=(double complex *)malloc(sizeof(double complex)*(MM+1));

  diag[0]=2.0;
  upper[0]=1.0;
  pot=x[1]-x[0];
  source_vec[0]=6.0/pot;
  cup=-y0_p+(y[1]-y[0])/pot;
  source_vec[0]=source_vec[0]*cup;
  for(i=1;i<MM;i++){
    pot=x[i+1]-x[i-1];
    upper[i]=(x[i+1]-x[i])/pot;
    lower[i-1]=1.0-upper[i];
    diag[i]=2.0;
    source_vec[i]=6.0/pot;
    cup=(y[i-1]-y[i])/(x[i]-x[i-1]);
    cup=cup+(y[i+1]-y[i])/(x[i+1]-x[i]);
    source_vec[i]=source_vec[i]*cup;
  }
  diag[MM]=2.0;
  lower[MM-1]=1.0;
  pot=x[MM]-x[MM-1];
  source_vec[MM]=6.0/pot;
  cup=yMM_p+(y[MM-1]-y[MM])/pot;
  source_vec[MM]=source_vec[MM]*cup;

  tridiag(MM+1,lower,diag,upper,source_vec,moments);

  free(lower);
  free(diag);
  free(upper);
  free(source_vec);
}

/*  From the computed moments (second derivs at given points), */
/*  produce the interpolating function at the point xx. */

double complex cub_spline(double complex *moments,double *x,
						  double complex *y,double xx){
  double complex cup1,cup2,cup3,cup4,z,pot;
  double h,w;
  int i;

  for(i=0;x[i]<=xx;i++);

  h=x[i]-x[i-1];
  w=xx-x[i-1];
  cup1=y[i-1];
  pot=((y[i]-y[i-1])/h)-((2.0*moments[i-1]+moments[i])*h/6.0);
  cup2=pot*w;
  cup3=0.5*moments[i-1]*w*w;
  cup4=w*w*w*(moments[i]-moments[i-1])/6.0/h;
  z=cup1+cup2+cup3+cup4;

  return z;
}

/*  From the computed moments (second derivs at given points), */
/*  produce the first derivative of the interpolating function */
/*  at the point xx. */

double complex d_cub_spline(double complex *moments,double *x,
							double complex *y,double xx){
  double complex cup2,cup3,cup4,z,pot;
  double h,w;
  int i;

  for(i=0;x[i]<=xx;i++);

  h=x[i]-x[i-1];
  w=xx-x[i-1];
  pot=((y[i]-y[i-1])/h)-((2.0*moments[i-1]+moments[i])*h/6.0);
  cup2=pot;
  cup3=moments[i-1]*w;
  cup4=w*w*(moments[i]-moments[i-1])/2.0/h;
  z=cup2+cup3+cup4;

  return z;
}

/*  From the computed moments (second derivs at given points), */
/*  produce the second derivative of the interpolating function */
/*  at the point xx. */

double complex dd_cub_spline(double complex *moments,double *x,
							 double complex *y,double xx){
  double complex cup3,cup4,z;
  double h,w;
  int i;

  for(i=0;x[i]<=xx;i++);

  h=x[i]-x[i-1];
  w=xx-x[i-1];
  cup3=moments[i-1];
  cup4=w*(moments[i]-moments[i-1])/h;
  z=cup3+cup4;

  return z;
}


/*              */
/* ODE solvers. */
/*              */


/* Solves  f'' = (-funkt + epsilon) f on [a,b]. funkt(x) is complex. */
/* A Runge-Kutta routine is used. */
/* Produces the soln satisfying f'(a)=-bc, f(a)=1. */
/* Output is stored in vec, a 2 by N+1 vector. */
/* vec[0] contans f, vec[1] contains f' */

void rk_bc_wave_eq_ode(int N, double a, double b, double complex eps,
					   double (*funct)(double), double complex bc,
					   double complex **vec)
{
  double complex pot,**k;
  double dz;
  int i;

  dz=(b-a)/N;
  k=(double complex **)matrix_alloc(2,4);

  vec[1][0]=-bc;
  vec[0][0]=1.0;

  for(i=1;i<=N;i++){
    pot=-FUNKT(a+(i-1)*dz)+eps;
    k[0][0]=pot*vec[0][i-1]*dz;
    k[1][0]=vec[1][i-1]*dz;
    pot=-FUNKT(a+(i-0.5)*dz)+eps;
    k[0][1]=pot*(vec[0][i-1]+0.5*k[1][0])*dz;
    k[1][1]=(vec[1][i-1]+0.5*k[0][0])*dz;
    k[0][2]=pot*(vec[0][i-1]+0.5*k[1][1])*dz;
    k[1][2]=(vec[1][i-1]+0.5*k[0][1])*dz;
    pot=-FUNKT(a+i*dz)+eps;
    k[0][3]=pot*(vec[0][i-1]+k[1][2])*dz;
    k[1][3]=(vec[1][i-1]+k[0][2])*dz;
    pot=k[0][0]/6.0+k[0][1]/3.0+k[0][2]/3.0+k[0][3]/6.0;
    vec[1][i]=vec[1][i-1]+pot;
    pot=k[1][0]/6.0+k[1][1]/3.0+k[1][2]/3.0+k[1][3]/6.0;
    vec[0][i]=vec[0][i-1]+pot;
  }

  matrix_free(k,2,4);
}

/* Solves  f'' = (-funkt + epsilon) f on [a,b]. funkt(x) is complex. */
/* A Runge-Kutta routine is used. */
/* Produces the soln satisfying f'(b)=-bc, f(b)=1. */
/* Output is stored in vec, a 2 by N+1 vector. */
/* vec[0] contans f, vec[1] contains f' */

void reverse_rk_bc_wave_eq_ode(int N,double a,double b,double complex eps,
			       double (*funct)(double),double complex bc,
			       double complex **vec)
{
  double complex pot,**k;
  double dz;
  int i;

  dz=(b-a)/N;
  k=(double complex **)matrix_alloc(2,4);

  vec[1][N]=bc;
  vec[0][N]=1.0;

  for(i=N-1;i>=0;i--){
    pot=-FUNKT(a+(i+1)*dz)+eps;
    k[0][0]=pot*vec[0][i+1]*dz;
    k[1][0]=vec[1][i+1]*dz;
    pot=-FUNKT(a+(i+0.5)*dz)+eps;
    k[0][1]=pot*(vec[0][i+1]+0.5*k[1][0])*dz;
    k[1][1]=(vec[1][i+1]+0.5*k[0][0])*dz;
    k[0][2]=pot*(vec[0][i+1]+0.5*k[1][1])*dz;
    k[1][2]=(vec[1][i+1]+0.5*k[0][1])*dz;
    pot=-FUNKT(a+i*dz)+eps;
    k[0][3]=pot*(vec[0][i+1]+k[1][2])*dz;
    k[1][3]=(vec[1][i+1]+k[0][2])*dz;
    pot=k[0][0]/6.0+k[0][1]/3.0+k[0][2]/3.0+k[0][3]/6.0;
    vec[1][i]=vec[1][i+1]+pot;
    pot=k[1][0]/6.0+k[1][1]/3.0+k[1][2]/3.0+k[1][3]/6.0;
    vec[0][i]=vec[0][i+1]+pot;
  }

  for(i=0;i<=N;i++) vec[1][i]=-vec[1][i];

  matrix_free(k,2,4);
}

/*                                                  */
/* calculate haversine distance for linear distance */
/*                                                  */

double haversine_km_h(double lat1, double long1, double lat2, double long2){
  /* CHH Edit removed unused variables dist, sgn_lat, vert, cup_v, dlat */
  double dlong,cup_h,horiz,sgn_long;

  dlong=(long2-long1)*d2r;
  sgn_long=dlong/fabs(dlong);
  cup_h=cos(lat1*d2r)*cos(lat2*d2r)*sin(0.5*dlong)*sin(0.5*dlong);
  horiz=sgn_long*6367.0*2.0*atan2(sqrt(cup_h),sqrt(1.0-cup_h));

  return horiz;
}

double haversine_km_v(double lat1, double long1, double lat2, double long2){
  /* CHH edit removed unused variables dist, sgn_long, horiz, cup_h, dlong */
  double dlat,cup_v,vert,sgn_lat;

  dlat=(lat2-lat1)*d2r;
  sgn_lat=dlat/fabs(dlat);
  cup_v=sin(0.5*dlat)*sin(0.5*dlat);
  vert=sgn_lat*6367.0*2.0*atan2(sqrt(cup_v),sqrt(1.0-cup_v));

  return vert;
}

double haversine_km(double lat1, double long1, double lat2, double long2){
  double dlong,dlat,cup_h,cup_v,sgn_long,sgn_lat,dist,vert,horiz;

  dlong=(long2-long1)*d2r;
  sgn_long=dlong/fabs(dlong);
  cup_h=cos(lat1*d2r)*cos(lat2*d2r)*sin(0.5*dlong)*sin(0.5*dlong);
  horiz=sgn_long*6367.0*2.0*atan2(sqrt(cup_h),sqrt(1.0-cup_h));

  dlat=(lat2-lat1)*d2r;
  sgn_lat=dlat/fabs(dlat);
  cup_v=sin(0.5*dlat)*sin(0.5*dlat);
  vert=sgn_lat*6367.0*2.0*atan2(sqrt(cup_v),sqrt(1.0-cup_v));

  dist=6367.0*2.0*atan2(sqrt(cup_h+cup_v),sqrt(1.0-cup_h-cup_v));

  return dist;
}
