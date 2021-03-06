#include "densities.h"

#include <math.h>
#include <stdio.h>
#include "../Bessel/bessel.h"
#include "../Constants/physicalConstants.h"
#include "../Constants/programConstants.h"

double getrhoGamma(double T9){//photon energy density
  return ar*T9*T9*T9*T9;
}
double getrhoNu(double h,double T9){//neutrino energy density
  // This takes into account the fact that the temperature is
  // different between the baryons/radiation and the neutrinos
  return XNU*(7./8.)*(ar)*pow(T0,4.)*pow(h*T9*T9*T9/rhob0,4./3.);
}
double getrhoe(double phie,double T9){//lepton energy density
  double sum=0, z=mecc/kb/T9,re;//temp; unitless;final answer
  int i;//iteration variable
  re = 2./PI/PI*pow(mecc,4)/pow(hbar*c,3)/c/c;
  for(i=1;i<itNum;i++){
    sum+=pow(-1.,i+1)*cosh(i*phie)*MBess(i*z);
    //printf("i=%i\tcosh(i*phie)=%e\n",i,cosh(i*phie));//cosh(i*phie));
  }
  //printf("re=%e\tphie=%e\tre*sum=%e\n",re,phie,re*sum);
  
  return re*sum;
}
double getrhob(double h,double T9,double *Y){//baryon energy density
  int i;//iteration variable
  double sum=0;//summation
  double rb=h*T9*T9*T9;//final answer
  for(i=0;i<nnuc;i++)//add up the part from each nuclei
    sum+=(deltaM[i]/M_u+zeta*T9)*Y[i];
  return rb*(1.+sum);
}

double getH(double rho){//expansion rate
  return sqrt(8.*PI*G*rho/3.);
}

/*Pressures*/
double getPGamma(double T9){//per c^2
  return getrhoGamma(T9)/3.;
}
double getPNu(double T9){//per c^2
  return 7.*getPGamma(T9)/8.;
}
double getPe(double phie,double T9){//per c^2
  double sum=0, z=mecc/kb/T9,re;//temp; unitless;final answer
  int i;//iteration variable
  re = 2./PI/PI*pow(mecc,4)/pow(hbar*c,3)/c/c;
  for(i=1;i<itNum;i++)
    sum+=pow(-1.,i+1)/i*cosh(i*phie)*LBess(i*z);
  return re*sum/z;
}
double getPb(double h, double T9, double *Y){//per c^2
  double sum=0;
  int i;
  for(i=0;i<nnuc;i++)
    sum+=Y[i];
  return h*T9*T9*T9*2./3.*zeta*T9*sum;
}
