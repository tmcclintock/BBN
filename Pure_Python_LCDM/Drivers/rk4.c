#include "rk4.h"

/* This is a general implementation of the 4th order Runga-Kutta algorithm.
 */

int rk4(int n,double*y,double*yout,double*t_p,double*dt_p,double*dtprev_p){

  int fail = 0; ///Success checker

  int i;//iteration variable
  double t = *t_p;
  double dt = *dtprev_p;//
  double k1[n],k2[n],k3[n],k4[n],temp[n],y_old[n];

  for(i=0;i<n;i++)
    y_old[i]=y[i];

  fail = get_dydt(k1,y);
  if (fail) return 1;
  fail = get_dt(k1,y,dt_p,dt);
  if (fail) return 1;
  dt = *dt_p;
  for(i=0;i<n;i++){
    temp[i]=y[i]+k1[i]*dt/2.;
    //if(temp[i]<Ymin)temp[i]=Ymin;
  }

  fail = get_dydt(k2,temp);
  if (fail) return 1;
  for(i=0;i<n;i++){
    temp[i]=y[i]+k2[i]*dt/2.;
    //if(temp[i]<Ymin)temp[i]=Ymin;
  }

  fail = get_dydt(k3,temp);
  if (fail) return 1;
  for(i=0;i<n;i++){
    temp[i]=y[i]+k3[i]*dt/2.;
    //if(temp[i]<Ymin)temp[i]=Ymin;
  }

  fail = get_dydt(k4,temp);
  if (fail) return 1;
  for(i=0;i<n;i++){
    yout[i] = y[i]+dt*(k1[i]/6.+k2[i]/3.+k3[i]/3.+k4[i]/6.);
  }

  return 0;
}
