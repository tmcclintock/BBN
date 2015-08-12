#include "outputDynamics.h"

#include <stdio.h>
#include <stdlib.h>
#include "../Constants/programConstants.h"

void outputDynamics(FILE *dp,double T9,double h,double phie,double t,double H){
  fwrite(&t,sizeof(double),1,dp);
  fwrite(&T9,sizeof(double),1,dp);
  fwrite(&h,sizeof(double),1,dp);
  fwrite(&phie,sizeof(double),1,dp);
  fwrite(&H,sizeof(double),1,dp);
}
