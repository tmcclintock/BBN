//location "../InitialConditions/"

#ifndef INITABUN_H
#define INITABUN_H

double getYpInit(double T9);//initial proton abundance
double getYnInit(double T9);//initial neutron abundance
double getYiInit(double T9,int i);//initial abundances
void initializeAbundances(double *Y,double T9);//wrapper

#endif
