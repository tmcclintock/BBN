#include "initialAbundances.h"

#include <math.h>
#include <stdio.h>
#include "../Constants/physicalConstants.h"

double getYpInit(double T9){
  return 1./(1+exp(-QPN/kb/T9));
}

double getYnInit(double T9){
  return 1./(1+exp(QPN/kb/T9));
}

double getYiInit(double T9, int i){
  if(i==0)return getYpInit(T9);
  if(i==1)return getYnInit(T9);
  return Ymin;
}

void initializeAbundances(double *Y, double T9){
  int i;
  Y[0]=getYnInit(T9);
  Y[1]=getYpInit(T9);
  Y[2]=Y[0]*Y[1]*rhob0*exp(Qvals[3]/T9)/(pow(T9,1.5)*4.71e9);
  for(i=3;i<totalnnuc;i++)
    Y[i]=Ymin;//getYiInit(T9,i);

  /* Test Print */
/*   printf("Initial Abundances:  "); */
/*   for(i=0;i<nnuc;i++) */
/*     printf("Y[%i]=%e\t",i,Y[i]); */
/*   printf("\n"); */
}
