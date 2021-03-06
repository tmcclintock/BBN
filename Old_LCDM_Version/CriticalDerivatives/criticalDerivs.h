#ifndef CRITICALD_H
#define CRITICALD_H

/*Returns the expansion rate, overlaps in functionality with getDerivs, 
  but it is necessary*/
double getExpansionRate(double phie,double h,double T9,double *Y);

/*Returns an array pointing to each of the critical derivatives*/
double getDerivs(double phie, double h, double T9, double *Y, double *dYdt, double *v);

/*The critical derivatives*/
double dTdt(double H,double rho,double pres,double rhoe,double S,double dSdt,double N,double M,double phie,double h,double T9,double *Y,double *dYdt);
double dhdt(double Tdot,double H,double h,double T9);
double dphiedt(double phie,double Tdot,double H,double N,double M,double T9,double S,double dSdt);

/*Temperature derivative of r = log(a^3)*/
double drdT(double H,double rho,double pres,double rhoe,double S,double dSdt,double N,double M,double phie,double h,double T9,double *Y,double *dYdt);

#endif
