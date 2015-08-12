#include "outputRates.h"

#include <stdio.h>
#include <stdlib.h>
#include "../Constants/programConstants.h"

void outputRates(FILE *rp,double *dYdt,double t){
  int i;
  fwrite(&t,sizeof(double),1,rp);//the time
  for(i=0;i<nnuc;i++)
    fwrite(&dYdt[i],sizeof(double),1,rp);//abundance rates
}
