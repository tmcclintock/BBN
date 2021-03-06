#include "outputAbundances.h"

#include <stdio.h>
#include <stdlib.h>
#include "../Constants/programConstants.h"

void outputAbundances(FILE *ap,double *Y,double t,double T9){
  int i;
  double y;
  fwrite(&t,sizeof(double),1,ap);//the time
  fwrite(&T9,sizeof(double),1,ap);//the temperature
  for(i=0;i<nnuc;i++){
    y=Y[i]/Y[1];
    if(i==5)y*=4;//express 4He as a mass fraction
    fwrite(&y,sizeof(double),1,ap);//abundances
  }
}
