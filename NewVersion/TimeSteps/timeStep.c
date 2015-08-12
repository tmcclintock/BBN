#include "timeStep.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../Constants/physicalConstants.h"
#include "../Constants/programConstants.h"

double getTtimeStep(double T9, double Tdot){
  return fabs(T9/Tdot)*CT;
}
double getYtimeStep(double *Y, double *dYdt){
  double dt=CY*fabs(Y[0]/dYdt[0])*(1+pow(log10(Y[0])/log10(Ymin),2));
  double dt_temp=dt;
  int i,imin=0;
  for(i=0;i<nnuc;i++){
    if(dYdt[i]>0&&Y[i]>Ymin){
      dt_temp = CY*fabs(Y[i]/dYdt[i])*(1+pow(log10(Y[i])/log10(Ymin),2));
    }
    if(dt_temp<dt){dt=dt_temp;imin=i;}
  }
  //printf("imin=%i\tY=%e\tdYdt=%e\n",imin,Y[imin],dYdt[imin]);
  return dt;
}

double getTimeStep(double T9,double h,double phie, double *dvdt,double *Y,double *dYdt,double dtprev){
  double Tdot=dvdt[0];
  double dt=0;
  double dt_T = getTtimeStep(T9,Tdot);
  double dt_Y = getYtimeStep(Y,dYdt);
  if(dt_T!=dt_T||dt_Y!=dt_Y){
    printf("Error: NAN");
    return -1;
  }
  if(dt_T<dt_Y){
    //    printf("T step= %e\n",dt_T);
    dt=dt_T;}
  
  else {dt=dt_Y;
    //    printf("Y step = %e\n",dt_Y,dt_T);
  }
  if(dt>TIMECHANGE*dtprev){
    //printf("dtT=%e\tdtY=%e\tdtPREV=%e\n",dt_T,dt_Y,TIMECHANGE*dtprev);  
    return TIMECHANGE*dtprev;
  }//The timestep cannot grow by more than 150%
  //printf("dtT=%e\tdtY=%e\n",dt_T,dt_Y);
  return dt;
}
