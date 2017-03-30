#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>

#define PI M_PI
#define RR 50.0

double intersection_toy_model_curve(double phi);
int write_sphere_grid(char *filename);
int write_theta_phi_grid(char *filename);
int write_theta_phi_vertices(char *filename);

int N_phi,N_phi0,N_the;

int main(int argc,char **argv){
  int status;

  N_the=25;
  N_phi0=30;
  if(N_phi0<=N_the){
    printf("azimuthal array smaller than elevation array\n");
    exit(0);
  }

  status=write_sphere_grid("sphere_section_grid.dat");

  status=write_theta_phi_grid("theta_phi_grid.dat");
  status=write_theta_phi_vertices("theta_phi_grid1.dat");

  return(status);
}

int write_sphere_grid(char *filename){
  double phi,theta,phi_skip,the_skip,cup;
  int i,j,N_phi;
  FILE *f;

  the_skip=(double)(0.5*PI/N_the);

  f=fopen(filename,"w");
  cup=1.0;
  for(j=1;j<N_the;j++){
    N_phi=N_phi0+1-j;
    phi_skip=(double)(2.0*PI/N_phi);
    for(i=0;i<=N_phi;i++){
      phi=i*phi_skip;
      theta=0.5*PI-(cup*intersection_toy_model_curve(phi)+j*the_skip);
      fprintf(f,"%e   %e   %e\n",RR*sin(theta)*cos(phi),RR*sin(theta)*sin(phi),RR*cos(theta));
    }
    fprintf(f,"\n");
    cup=0.8*cup;
  }
  fprintf(f,"%e   %e   %e\n",0.0,0.0,RR);
  fclose(f);

  return(0);
}

int write_theta_phi_vertices(char *filename){
  double phi,theta,phi_skip,the_skip,cup;
  int i,j,N_phi;
  FILE *f;

  the_skip=(double)(0.5*PI/N_the);

  f=fopen(filename,"w");
  cup=1.0;
  for(j=0;j<N_the;j++){
    N_phi=N_phi0+1-j;
    phi_skip=(double)(2.0*PI/N_phi);
    for(i=0;i<=N_phi;i++){
      phi=i*phi_skip;
      theta=0.5*PI-(cup*intersection_toy_model_curve(phi)+j*the_skip);
      fprintf(f,"%e   %e\n",phi,theta);
    }
    fprintf(f,"\n");
    cup=0.8*cup;
  }
  fprintf(f,"%e   %e\n",0.0,0.0);
  fclose(f);

  return(0);
}

int write_theta_phi_grid(char *filename){
  double phi,theta,phi_skip,the_skip,cup;
  int i,j,N_phi;
  FILE *f;

  the_skip=(double)(0.5*PI/N_the);

  f=fopen(filename,"w");
  cup=1.0;
  for(j=0;j<N_the;j++){
    N_phi=N_phi0+1-j;
    phi_skip=(double)(2.0*PI/N_phi);
    for(i=0;i<N_phi;i++){
      phi=i*phi_skip;
      theta=0.5*PI-(cup*intersection_toy_model_curve(phi)+j*the_skip);
      fprintf(f,"%e   %e\n",phi,theta);
      phi=i*(double)(2.0*PI/(N_phi-1));
      theta=0.5*PI-(0.8*cup*intersection_toy_model_curve(phi)+(j+1)*the_skip);
      fprintf(f,"%e   %e\n",phi,theta);
    }
    fprintf(f,"\n");
    cup=0.8*cup;
  }
  fprintf(f,"%e   %e\n",0.0,0.0);
  fprintf(f,"%e   %e\n",2*PI,0.0);
  fclose(f);

  return(0);
}

double intersection_toy_model_curve(double phi){
  double answer;

  answer=0.2*exp(-(pow(0.6*(phi-1.2*PI),6)));

  return answer;
}
