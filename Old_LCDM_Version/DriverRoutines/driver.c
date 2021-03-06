#include "driver.h"

#include <stdlib.h>
#include <stdio.h>
#include "../Constants/physicalConstants.h"
#include "../Constants/programConstants.h"
#include "../OutputRoutines/outputAbundances.h"
#include "../OutputRoutines/outputDynamics.h"
#include "../OutputRoutines/outputHeader.h"
#include "../OutputRoutines/outputRates.h"
#include "rk4BBN.h"

int driver(FILE *hp,FILE *ap,FILE *dp,FILE *rp,double T9,double h,double phie,double t,double *Y,double *dYdt){

  /*Iteration variable*/
  int i;
  
  /*Declare nsteps: the number of steps taken to reach TFinal*/
  int nsteps=0;
  int nprints=0;
  int inc=50;//how often we need to add a correction
  int ip=50;

  /*Declare Expansion rate*/
  double H;

  /*Declare timestep*/
  double dt, dtprev=1e-4;
 
  /*End at the final temperature or after some number of steps*/
  while(T9>T9Final&&nsteps<10){//begin while T<TFinal
    //while(nsteps<1){
    dt = rk4BBN(&T9,&h,&phie,&t,&H,Y,dYdt,dtprev,inc,ip);
    if(dt<0)return -1;//error
    if(nsteps%1==0){
      
      /*       printf("T9=%16.15e\t",T9); */
      //printf("t=%4.3e\tdt=%e\n",t,dt);
      printf("dt=%4.3e\teh=%4.3e\tphie=%4.3e\tnb=%4.3e\tT9=%4.3e\n",dt,h,phie,h*T9*T9*T9/M_ug,T9);
      outputAbundances(ap,Y,t,T9);
      outputDynamics(dp,T9,h,phie,t,H);
      outputRates(rp,dYdt,t);
      nprints++;
    }
    if(ip==inc) ip=0;
    ip++;
    nsteps++;
    dtprev=dt;
  }//end while T<TFinal

  /* Print all abundnaces (for debugging and such) */
  for(i=0;i<nnuc;i++){
    if(i==5){
      printf("X[%i]=%e\t",i,Y[i]*4);
      //	  printf("dYdt[%i]=%e\n",i,dYdt[i]);
    } else if(i==1){printf("Y[%i]=%4.3e\t",i,Y[i]);
    } else {printf("RH[%i]=%4.3e\t",i,Y[i]/Y[1]);
    }
  }
  printf("\n");

  /* Print meaningful abundances */
  printf("Y_n    = %5.4e\n",Y[0]/Y[1]);
  printf("Y'_H   = %5.4e\n",Y[1]);
  printf("D/H    = %5.4e\n",Y[2]/Y[1]);
  printf("He3/H  = %5.4e\n",(Y[3]+Y[4])/Y[1]);
  printf("X_He4  = %5.4e\n",Y[5]*A[5]);
  printf("Li7/H  = %5.4e\n",(Y[6]+Y[7])/Y[1]);
  printf("Li6/H  = %5.4e\n",Y[8]/Y[1]);
  printf("Li6/Li7= %5.4e\n",Y[8]/Y[7]);


  printf("\nFINAL T9=%e\th=%e\tphie=%e\tt=%e\tnsteps=%i\n",T9,h,phie,t,nsteps);

  outputHeader(hp,nnuc,nreac,nprints);
  return 0;
}

