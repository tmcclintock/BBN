#include "densDerivs.h"

#include <math.h>
#include <stdio.h>
#include "../Bessel/bessel.h"
#include "../Constants/physicalConstants.h"
#include "../Constants/programConstants.h"
#include "../PartialDerivatives/partials.h"
#include "densities.h"

double drhobdt(double h, double T9, double *dYdt){
  double sum=0;
  int i;
  for(i=0;i<nnuc;i++)
    sum+=(zeta*T9+deltaM[i]/M_u)*dYdt[i];
  return h*T9*T9*T9*sum;
}
double drhoedt(double rhoe,double phie,double H, double h, double T9, double S, double dSdt,double N,double M){//lepton dens. time deriv.
  return 3.*H*drhoedphie(phie,T9)*(dphiedr(phie,N,M,T9,S)+dphiedS(phie,N,M,T9,S)*dSdt/3./H);
}


double drhoGdT(double T9){//photon density temp. deriv.
  return 4.*getrhoGamma(T9)/T9;
}
double drhoNdT(double T9){//neutrino density temp. deriv
  return 7.*getrhoGamma(T9)/2./T9;
}
double drhoedT(double rhoe,double phie, double T9){//lepton density temp deriv.
  double z = mecc/kb/T9, sum=0;
  int i;
  for(i=1;i<itNum;i++){
    sum+=pow(-1.,i+1)*cosh(i*phie)*(3*KN(4,i*z)/4.+KN(2,i*z)+KN(0,i*z)/4.);
  }
  return kb*z/mecc*(rhoe+1./PI/PI*pow(mecc,4)/pow(hbar*c,3)/c/c*sum);
}
double drhobdT(double h, double T9, double *Y){//baryon dens. temp. deriv.
  double sum=0;
  int i;
  for(i=0;i<nnuc;i++)
    sum+=Y[i];
  return h*T9*T9*T9*zeta*sum;
}

double drhoedphie(double phie,double T9){//lepton dens. chem. pot. deriv.
  double z=mecc/kb/T9,sum=0;
  int i;
  for(i=1;i<itNum;i++)
    sum+=pow(-1.,i+1)*i*sinh(i*phie)*MBess(i*z);
  return 2./PI/PI*pow(mecc,4)/pow(hbar*c,3)/c/c*sum;
}
