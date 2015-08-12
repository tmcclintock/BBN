//location: "../PartialDerivatives"

#ifndef PARTIALS_H
#define PARTIALS_H

//Phi_e derivatives
double dphiedT(double phie,double N, double M, double T9, double S);
double dphiedr(double phie,double N, double M, double T9, double S);
double dphiedS(double phie,double N, double M, double T9, double S);

//M and N and their derivatives
double getM(double h, double T9, double S, double dSdt);
double getN(double phie, double T9);
double dNdphie(double phie,double T9);
double dNdT(double N,double phie,double T9);
double dMdT(double M,double T9);
double dMdr(double M);
double dMdS(double M,double S);
#endif
