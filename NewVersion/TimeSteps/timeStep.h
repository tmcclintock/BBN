#ifndef TIMESTEPS_H
#define TIMESTEPS_H

double getTtimeStep(double T9,double Tdot);
double getYtimeStep(double *Y,double *dYdt);
double getTimeStep(double T9,double h,double phie,double *dvdt,double *Y,double *dYdt,double dtprev);
#endif
