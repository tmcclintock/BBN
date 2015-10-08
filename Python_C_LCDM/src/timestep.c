#include "timestep.h"

/* This calculates the timestep given the abundances and the
   abundance derivatives, as well as the dynamic variables and their
   derivatives.
 */

double get_dt_Y(double*dYdt,double*Y){
  int i, imin;//count variables
  double dt=fabs(Y[3]/dYdt[3])*(1+pow(log10(Y[3])/log10(Ymin),2))*CY;
  double dt_temp = dt;
  for(i=3;i<nnuc;i++){
    if(dYdt[i]>0 && Y[i]>Ymin) //might be fabs(dYdt[i])
      dt_temp = fabs(Y[i]/dYdt[i])*(1+pow(log10(Y[i])/log10(Ymin),2))*CY;
    if(dt_temp < dt){
      dt=dt_temp;
      imin=i;
    }
  }
  return dt;
}

int get_dt(double*dydt,double*y,double*dt_p,double dt){
  double T9=y[0];
  double T9dot=dydt[0];
  double dt_prev = dt;//for now
  double dt_T = fabs(T9/T9dot)*CT;
  double dt_Y = get_dt_Y(dydt,y);
  //Do some error checking
  if(dt_T!=dt_T){
    printf("Error: dt_T == nan.\n");
    fflush(stdout);
    return 1;
  }
  if(dt_Y!=dt_Y){
    printf("Error: dt_T == nan.\n");
    fflush(stdout);
    return 1;
  }
  //Take the minimal timestep
  if(dt_T<dt_Y)
    dt = dt_T;
  else
    dt = dt_Y;
  if(dt>TIMECHANGE*dt_prev)
    dt = TIMECHANGE*dt_prev;
  *dt_p = dt;
  return 0;
}
