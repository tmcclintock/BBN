//location: "../Densities/"

#ifndef DENSDERIVS_H
#define DENSDERIVS_H

double drhoedt(double rhoe,double phie,double H,double h,double T9,double S,double dSdt,double N,double M);//lepton en. dens. dt
double drhobdt(double h, double T9, double *dYdt);//baryon en. dens. dt

double drhoGdT(double T9);//photon en dens. dT
double drhoNdT(double T9);//neutrino en dens dT
double drhoedT(double rhoe, double phie, double T9);//lepton en. dens dT
double drhobdT(double h, double T9, double *Y);//baryon en dens dT

double drhoedphie(double phie,double T9);//lepton en dens d(chem potential)


#endif
