#include "collect_derivs.h"

/* This combines the abundance derivativs and dynamic derivs into a single
   array so that it can be used in RK4.
 */

int get_dydt(double*dxdt,double*x){
  double v[3],dvdt[3];//dynamic variables
  double Y[totalnnuc],dYdt[totalnnuc];//abundances

  int fail = 0;//success checker
  fail = get_dynamics_derivs(dvdt,v,dYdt,Y);
  if (fail) return 1;
  fail = get_abundance_derivs(dYdt,Y,v);

  //Combine the derivatives into one array
  int i;
  for(i=0;i<3;i++)
    dxdt[i]=dvdt[i];
  for(i=3;i<3+totalnnuc;i++)
    dxdt[i]=dYdt[i-3];
  return 0;
}
