#include "updateVariables.h"

#include "../Constants/programConstants.h"

//This is used both for updating real variables and temporary variables
void updateVariables(double *t,double *T9,double *h,double *phie,double *dvdt,double *Y,double *dYdt,double dt){
  int i;//iteration variable

  *t+=dt;
  *T9=*T9+dvdt[0]*dt;
  *h=*h+dvdt[1]*dt;
  *phie=*phie+dvdt[2]*dt;
  for(i=0;i<nnuc;i++)
    Y[i]=Y[i]+dYdt[i]*dt;
}
