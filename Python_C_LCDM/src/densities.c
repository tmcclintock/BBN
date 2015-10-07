#include "densities.h"

/* This contains the densities, pressures, and density derivatives
   for the constituents.
 */

double getrhoGamma(double T9){//photons
  return ar*T9*T9*T9*T9;
}

double getrhoNu(double h,double T9){//neutrinos
  return XNU*(7./8.)*ar*pow(T0,4)*pow(h*T9*T9*T9/rhob0,4./3.);
}

double getrhoe(double phie,double T9){//leptons
  double sum=0, z=mecc/kb/T9,re;//temp; unitless; final answer
  int i;
  re = 2./PI/PI*pow(mecc,4)/pow(hbar*c,3)/c/c;
  for(i=1;i<itNum;i++)
    sum+=pow(-1,i+1)*cosh(i*phie)*MBess(i*z);
  return re*sum;
}

double getrhob(double h,double T9,double*Y){//baryons
  int i;
  double sum=0, rb=h*T9*T9*T9;//sum and final answer
  for(i=0;i<nnuc;i++)
    sum+=(deltaM[i]/M_u+zeta*T9)*Y[i];
  return rb*(1+sum);
}

double getPGamma(double T9){//photons; per c^2
  return getrhoGamma(T9)/3.;
}

double getPNu(double T9){//neutrinos; per c^2
  return 7.*getPGamma(T9)/8.;
}

double getPe(double phie,double T9){//leptons; per c^2
  double sum=0, z=mecc/kb/T9, re;//temp; unitless; final answer
  int i;
  re=2./PI/PI*pow(mecc,4)/pow(hbar*c,3)/c/c;
  for(i=1;i<itNum;i++)
    sum+=pow(-1,i+1)/i*cosh(i*phie)*LBess(i*z);
  return re*sum/z;
}

double getPb(double h,double T9,double*Y){//baryons; per c^2
  double sum=0;
  int i;
  for(i=0;i<nnuc;i++)
    sum+=Y[i];
  return h*T9*T9*T9*2./3.*zeta*T9*sum;
}
