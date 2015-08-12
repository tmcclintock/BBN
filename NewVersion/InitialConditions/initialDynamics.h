//location "../InitialConditions/"

#ifndef INITDYNAM_H
#define INITDYNAM_H

double getTInit();
double gethInit();
double getphieInit(double h,double T9,double Yp);
double getTimeInit(double T9);
void setrhob0(double T9,double h);

#endif
