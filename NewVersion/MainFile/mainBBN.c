/*Written by Tom McClintock
  Calculates the nuclear abundances from BBN
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "../Constants/physicalConstants.h"
#include "../Constants/programConstants.h"
#include "../DriverRoutines/driver.h"
#include "../InitialConditions/initializeAll.h"
#include "../InitialConditions/initialAbundances.h"
#include "../InitialConditions/initialDynamics.h"

/*These are global arrays and must be redeclared here*/
double deltaM[totalnnuc];
double Z[totalnnuc];
double Qvals[totalnreac];
double A[totalnnuc];
double nni[11];
double nnj[11];
double nnk[11];
double nnl[11];
double reactionDetails[totalnreac][8];

int main(){
  /*Create the files to hold the abundances,dynamics, rates, and header info*/
  FILE *headFP=fopen("BinaryDataFiles/header.dat","w");
  FILE *abunFP=fopen("BinaryDataFiles/abundances.dat","w");
  FILE *dynaFP=fopen("BinaryDataFiles/dynamics.dat","w");
  FILE *rateFP=fopen("BinaryDataFiles/rates.dat","w");

  /* Decalare the critical variables as well as the expansion rate H*/
  double T9,h,phie,t;

  /*Create the arrays for abundances and abundance changes*/
  double Y[totalnnuc], dYdt[totalnnuc];

  /*Initialize all variables*/
  initializeAll(&T9,&h,&phie,&t,Y);
  printf("Initial dynamics: Ti=%3.2e\thi=%3.2e\tphiei=%3.2e\trhob0=%e\tti=%e\n\n",T9,h,phie,rhob0,t);
  fflush(stdout);
  
  /*Initialize calculation constants*/
  initializeConstants(deltaM,Z,Qvals,A,nni,nnj,nnk,nnl,reactionDetails);  

  /*Call the driver wrapper*/
  int bbn = driver(headFP,abunFP,dynaFP,rateFP,T9,h,phie,t,Y,dYdt);
  
  if(bbn!=-1)
    printf("Simulation complete\n");
  else{
    printf("Error: NAN\n");
    return -1;
  }
  fflush(stdout);  

  return 0;
}
