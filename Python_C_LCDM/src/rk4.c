#include "rk4.h"

/* This is a general implementation of the 4th order Runga-Kutta algorithm.
 */

int rk4(double*y,double*yout,double*t_p,double*dt_p,double*dtprev_p,
	int inc, int ip){

  int n = 3+totalnnuc;//3 Dynamic variables plus all nuclei
  int fail = 0; ///Success checker

  int i;//iteration variable
  double t = *t_p;//Current time
  double dt = *dtprev_p;//dt holds the step size

  //Assemble Y0, the abundances at the start of the step
  double Y0[totalnnuc];
  for (i=0;i<totalnnuc;i++)
    Y0[i] = y[i+3];//first three spots are dynamic variables
  
  //The rest is standard RK4 with dynamic step size
  double k1[n],k2[n],k3[n],k4[n],temp[n];
  fail = get_dydt(k1,y,1,inc,ip,Y0,dt);
  if (fail) return 1;
  fail = get_dt(k1,y,dt_p,dt);
  if (fail) return 1;
  dt = *dt_p;
  for(i=0;i<n;i++){
    temp[i]=y[i]+k1[i]*dt/2.;
    if(temp[i]<Ymin)temp[i]=Ymin;
  }

  fail = get_dydt(k2,temp,2,inc,ip,Y0,dt);
  if (fail) return 1;
  for(i=0;i<n;i++){
    temp[i]=y[i]+k2[i]*dt/2.;
    if(temp[i]<Ymin)temp[i]=Ymin;
  }

  fail = get_dydt(k3,temp,3,inc,ip,Y0,dt);
  if (fail) return 1;
  for(i=0;i<n;i++){
    temp[i]=y[i]+k3[i]*dt/2.;
    if(temp[i]<Ymin)temp[i]=Ymin;
  }

  fail = get_dydt(k4,temp,4,inc,ip,Y0,dt);
  if (fail) return 1;
  for(i=0;i<n;i++){
    yout[i] = y[i]+dt*(k1[i]/6.+k2[i]/3.+k3[i]/3.+k4[i]/6.);
  }

  return 0;
}
