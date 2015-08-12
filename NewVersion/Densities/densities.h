//location: "../Densities/"
#ifndef DENSITIES_H
#define DENSITIES_H

double getrhoGamma(double T9);//photon energy density
double getrhoNu(double h,double T9);//neutrino energy density
double getrhoe(double phie,double T9);//electron + positron energy density
double getrhob(double h,double T9,double *Y);//baryon energy density
double getH(double rho);//expansion rate

/*Pressures*/
double getPGamma(double T9);
double getPNu(double T9);
double getPe(double phie,double T9);
double getPb(double h, double T9, double *Y);

#endif
