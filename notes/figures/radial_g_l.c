/* cc -I/usr/local/include radial_g_l.c -lgsl */

#include <stdio.h>
#include <gsl/gsl_sf_bessel.h>

int main (int argc,char** argv){
  int i,l;
  double x,y,dx;
  FILE *f;

  dx=0.01;

  printf("%s\n",argv[1]);
  f=fopen(argv[1],"w");
  for(i=0;i*dx<=20;i++){
    x=i*dx;
    fprintf(f,"%f", x);
    for(l=0;l<5;l++){
      y = gsl_sf_bessel_jl(l,x);
      fprintf (f,"   %.18e",y);
    }
    fprintf(f,"\n");
  }
  fclose(f);

  return 0;
}
