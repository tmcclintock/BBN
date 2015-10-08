#include "rk4BBN.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../AbundanceDerivs/geAbundanceDerivs.h"
#include "../Constants/physicalConstants.h"
#include "../Constants/programConstants.h"
#include "../CriticalDerivatives/criticalDerivs.h"
//#include "../Densities/densDerivs.h"
#include "../Densities/densities.h"
#include "../TimeSteps/timeStep.h"
//#include "updateVariables.h"

double rk4BBN(double *T9,double *h,double *phie,double *t,double *Hout,double *Y,double *dYdt,double dtprev,int inc,int ip){
  
  int fail=0;//correctness variable; 0 means success so far
  int i;//iteration variable

  /*Declare expansion rate*/
  double H;

  /*Declare variable holding the timestep*/
  double dt=dtprev;

  /*Declare the arrays holding dTdt dhdt and dphiedt*/
  //  double dvdt[3];
  double k1_dvdt[3],k2_dvdt[3],k3_dvdt[3],k4_dvdt[3];

  /*Declare temporary array and variables for abundances and critical values*/
  double Y_temp[totalnnuc];//, dYdt_total[totalnnuc];
  double T9_temp,h_temp,phie_temp;

  /*Declare temporary array for abundance changes*/
  double k1_dYdt[totalnnuc],k2_dYdt[totalnnuc],k3_dYdt[totalnnuc],k4_dYdt[totalnnuc];

  /*Calculate the initial Expansion Rate*/
  H = getExpansionRate(*phie,*h,*T9,Y);

  /*Begin with Y_temp = Y*/
  for(i=0;i<nnuc;i++)
    Y_temp[i]=Y[i];
  /*Note: in getdYdt Y is relabeled as Y0 and Y
    is used for Y_temp to reflect the fact that
    Y is the abundances at the RK4 step we are
    considering*/

  /*Calculate the abundance derivatives at time t*/
  /*Calculate the critical derivatives at time t*/
  /*These arrays are k1/h in RK4 lingo*/
  fail=getdYdt(*h,*T9,H,1,inc,ip,Y,Y_temp,k1_dYdt,dt);
  getDerivs(*phie,*h,*T9,Y,k1_dYdt,k1_dvdt);
  if(fail!=0)return -1;

  /*Calculate the current timestep*/
  dt = getTimeStep(*T9,*h,*phie,k1_dvdt,Y,k1_dYdt,dtprev);
  if(dt==-1)return -1;

  /*As one can see, the k1 term depends on dtprev, while all
    other k terms depend on the current dt. This is due to 
    simple fact that this is a "Verlet" RK4 scheme.*/

  /*Execute the runga-kutta algorithm*/
  /*Calculating the k2/h values*/
  /*printf("T9=%e  dTdt=%e\n",*T9,k1_dvdt[0]);
  printf("h= %e  dhdt=%e\n",*h,k1_dvdt[1]);
  printf("pe=%e  dpdt=%e\n",*phie,k1_dvdt[2]);*/
  for(i=0;i<nnuc;i++){
    //printf("Y[%i] = %e  dYdt[%i] = %e\n",i,Y[i],i,k1_dYdt[i]);
    Y_temp[i]= Y[i] + dt*k1_dYdt[i]/2;
    if(Y_temp[i]<Ymin)Y_temp[i]=Ymin;
  }
  T9_temp =     *T9 + dt*k1_dvdt[0]/2;
  h_temp =       *h + dt*k1_dvdt[1]/2;
  phie_temp = *phie + dt*k1_dvdt[2]/2;
  H = getExpansionRate(phie_temp,h_temp,T9_temp,Y_temp);
  fail=getdYdt(h_temp,T9_temp,H,2,inc,ip,Y,Y_temp,k2_dYdt,dt);
  getDerivs(phie_temp,h_temp,T9_temp,Y_temp,k2_dYdt,k2_dvdt);
  if(fail!=0)return -1;

  /*Calculate the k3/h values*/
  for(i=0;i<nnuc;i++){
    Y_temp[i]= Y[i] + dt*k2_dYdt[i]/2;
    if(Y_temp[i]<Ymin)Y_temp[i]=Ymin;
  }
  T9_temp =     *T9 + dt*k2_dvdt[0]/2;
  h_temp =       *h + dt*k2_dvdt[1]/2;
  phie_temp = *phie + dt*k2_dvdt[2]/2;
  H = getExpansionRate(phie_temp,h_temp,T9_temp,Y_temp);
  fail=getdYdt(h_temp,T9_temp,H,3,inc,ip,Y,Y_temp,k3_dYdt,dt);
  getDerivs(phie_temp,h_temp,T9_temp,Y_temp,k3_dYdt,k3_dvdt);
  if(fail!=0)return -1;


  /*Calculate the k4/h values*/
  for(i=0;i<nnuc;i++){
    Y_temp[i]= Y[i] + dt*k3_dYdt[i];
    if(Y_temp[i]<Ymin)Y_temp[i]=Ymin;
  }
  T9_temp =     *T9 + dt*k3_dvdt[0];
  h_temp =       *h + dt*k3_dvdt[1];
  phie_temp = *phie + dt*k3_dvdt[2];
  H = getExpansionRate(phie_temp,h_temp,T9_temp,Y_temp);
  fail=getdYdt(h_temp,T9_temp,H,4,inc,ip,Y,Y_temp,k4_dYdt,dt);
  getDerivs(phie_temp,h_temp,T9_temp,Y_temp,k4_dYdt,k4_dvdt);
  if(fail!=0)return -1;

  /*Find final critical variables*/
  double Tdot    = k1_dvdt[0]/6+k2_dvdt[0]/3+k3_dvdt[0]/3+k4_dvdt[0]/6;
  double hdot    = k1_dvdt[1]/6+k2_dvdt[1]/3+k3_dvdt[1]/3+k4_dvdt[1]/6;
  double phiedot = k1_dvdt[2]/6+k2_dvdt[2]/3+k3_dvdt[2]/3+k4_dvdt[2]/6;
  T9_temp =     *T9 + dt*Tdot;
  h_temp =       *h + dt*hdot;
  phie_temp = *phie + dt*phiedot;
/*   if(*T9>99.99)printf("\nT9=%4.3e\th=%4.3e\tphie=%4.3e\tH=%e\tTdot=%4.3e\thdot=%4.3e\tphiedot=%4.3e\t",*T9,*h,*phie,H,Tdot,hdot,phiedot); */

  
  /*Find total abundance changes*/
  for(i=0;i<nnuc;i++){
    dYdt[i]=k1_dYdt[i]/6. + k2_dYdt[i]/3. + k3_dYdt[i]/3. + k4_dYdt[i]/6.;
    //if(fabs(dYdt[i])<fabs(1e-5*H))dYdt[i]=0;//smaller than the expansion rate
  }

  /*Gather all abundances*/
  for(i=0;i<nnuc;i++){
    Y_temp[i]=  Y[i] + dt*dYdt[i];
  }

  /*Test print*/
/*   printf("\nTest:\t"); */
/*   for(i=0;i<nnuc;i++){ */
/*     printf("Y[%i]=%e\t",i,Y[i]); */
/*     printf("dYdt[%i]=%e\t",i,dYdt[i]); */
/*   } */
/*   printf("T9=%e\th=%e\tphie=%e dt=%e",T9_temp,h_temp,phie_temp,dt); */
/*   printf("\n"); */

  /*Update all parameters completely*/
  *T9=T9_temp,*h=h_temp,*phie=phie_temp;
  for(i=0;i<nnuc;i++){
    Y[i]=  Y_temp[i];
    if(Y[i]<Ymin)Y[i]=Ymin;//nothing is less than the minimum abundance
  }
  /*Test print 2*/
/*   printf("\nSuccess:\t"); */
/*   for(i=0;i<nnuc;i++){ */
/*     printf("Y[%i]=%e\t",i,Y[i]); */
/*     printf("dYdt[%i]=%e\t",i,dYdt[i]); */
/*   } */
/*   printf("\n"); */

  /*Update the expansion rate*/
  *Hout = H;

  /*Update t*/
  *t+=dt;
  return dt;
}
