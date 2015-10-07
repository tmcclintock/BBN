#include "expansion_rate.h"

/* The expansion rate is calculated here
 */

double get_H_wrapper(double T9,double h,double phie,double*Y){
  double rhoG = getrhoGamma(T9), rhoNu = getrhoNu(h,T9),
    rhoe = getrhoe(phie,T9), rhob = getrhob(h,T9,Y);
  double rhoTotal = rhoG+rhoNu+rhoe+rhob;
  return getH(rhoTotal);
}

double get_H(double rho){//expansion rate
  return sqrt(8.*PI*G*rho/3.);
}

