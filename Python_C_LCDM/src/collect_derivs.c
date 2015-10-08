#include "collect_derivs.h"

/* This combines the abundance derivativs and dynamic derivs into a single
   array so that it can be used in RK4.
 */

int get_dydt(double*dxdt,double*x,int ki,int inc,int ip,double*Y0,double dt){
  double v[3],dvdt[3];//dynamic variables
  double Y[totalnnuc],dYdt[totalnnuc];//abundances

  //Put the x values into Y and v
  v[0]=x[0],v[1]=x[1],v[2]=x[2];
  int i;
  for(i=0;i<totalnnuc;i++)
    Y[i]=x[i+3];

  int fail = 0;//success checker
  fail = get_abundance_derivs(dYdt,Y,v,ki,inc,ip,Y0,dt);
  if (fail) return 1;
  fail = get_dynamics_derivs(ki,dvdt,v,dYdt,Y);
  if (fail) return 1;

  //Combine the derivatives into one array
  for(i=0;i<3;i++)
    dxdt[i]=dvdt[i];
  for(i=3;i<3+totalnnuc;i++)
    dxdt[i]=dYdt[i-3];
  /*for(i=0;i<3+totalnnuc;i++){
    if(ki==1||ki==1){//||ki==3||ki==4){
    printf("%e\t%e\n",x[i],dxdt[i]);
    }
    }
    printf("\n");
  */
  return 0;
}
