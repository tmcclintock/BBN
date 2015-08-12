#include "initialDynamics.h"

#include <math.h>
#include <stdio.h>
#include "../Constants/physicalConstants.h"
#include "../Constants/programConstants.h"
#include "../Bessel/bessel.h"

double rhob0;

double getTInit(){
  return T0;
}
double gethInit(){//g/cm^3 apparently, fuck kawano's code
  return 33683*eta0*2.75;//1.889e31*eta0;//mol*MeV/cm/cm/cm/GK/GK/GK
}
double getphieInit(double h,double T9,double Yp){
  double sum=0,z=mecc/kb/T9;
  int i;
  for(i=1;i<itNum;i++)
    sum+=pow(-1.,i+1)*i*LBess(i*z);
  return PI*PI/2.*NA*pow(hbar*c/kb,3)*h*Yp/z/z/z/sum;//unitless somehow, but really fucking grams
}
double getTimeInit(double T9){//seconds
  return 10.4*10.4/T9/T9;
}

void setrhob0(double T9,double h){
  rhob0=h*T9*T9*T9;//sets the initial baryon mass density
}
