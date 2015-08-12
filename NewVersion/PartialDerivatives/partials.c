#include "partials.h"

#include <math.h>
#include "../Bessel/bessel.h"
#include "../Constants/physicalConstants.h"
#include "../Constants/programConstants.h"

//Phi_e derivatives
double dphiedT(double phie,double N,double M,double T9,double S){
  return (dMdT(M,T9)-dNdT(N,phie,T9))/dNdphie(phie,T9);
}
double dphiedr(double phie,double N,double M,double T9,double S){
  return dMdr(M)/dNdphie(phie,T9);
}
double dphiedS(double phie,double N,double M,double T9,double S){
  return dMdS(M,S)/dNdphie(phie,T9);
}


//M and N and their derivatives
double getM(double h, double T9, double S, double dSdt){
  return PI*PI/2*NA*pow(hbar*c/kb,3)*h*S;
}
double getN(double phie, double T9){
  double N=0,z=mecc/kb/T9;//unitless
  int i;//iteration variable
  for(i=1;i<itNum;i++)
    N+=pow(-1.,i+1)*sinh(i*phie)*LBess(i*z);
  return N*z*z*z;
}
double dNdphie(double phie,double T9){
  double dN=0,z=mecc/kb/T9;//unitless
  int i;//iteration variable
  for(i=1;i<itNum;i++)
    dN+=pow(-1.,i+1)*i*cosh(i*phie)*LBess(i*z);
  return dN*z*z*z;
}
double dNdT(double N,double phie,double T9){
  double dN,z=mecc/kb/T9;//unitless
  double sum=0;//answer and summation portion
  int i;//iteration variable
  for(i=1;i<itNum;i++)
    sum+=pow(-1.,i+1)*sinh(i*phie)*(KN(3,i*z)+K1(i*z));
  dN = -kb/mecc*(2*z*N-z*z*z*z/2.*sum);
  return dN;
}
double dMdT(double M,double T9){
  return -3.*M/T9;
}
double dMdr(double M){
  return -M;
}
double dMdS(double M,double S){
  if(S==0)return 0;
  return M/S;
}
