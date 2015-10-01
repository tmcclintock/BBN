#include "dynamics_derivs.h"

/* This calculates the derivatives of the dynamical variables T9, h, and phie.
 */

/*Prototypes for critical derivatives and 
  the temperature derivative of r = log(a^3)
*/
double dTdt(double H,double rho,double pres,double rhoe,
	    double S,double dSdt,double N,double M,
	    double phie,double h,double T9,
	    double *Y,double *dYdt);
double dhdt(double Tdot,double H,double h,double T9);
double dphiedt(double phie,double Tdot,double H,
	       double N,double M,double T9,double S,double dSdt);
double drdT(double H,double rho,double pres,double rhoe,
	    double S,double dSdt,double N,double M,
	    double phie,double h,double T9,
	    double *Y,double *dYdt);

int get_dynamics_derivs(double*dvdt,double*v,double*dYdt,double*Y){
  double T9=v[0],h=v[1],phie=v[2];

  //Calculate all densities and the expansion rate
  double rhoG = getrhoGamma(T9), rhoNu = getrhoNu(h,T9),
    rhoe = getrhoe(phie,T9), rhob = getrhob(h,T9,Y);
  double rhoTotal = rhoG+rhoNu+rhoe+rhob;
  double rho = rhoG+rhoe; //neutrinos have decoupled
  double H = getH(rhoTotal); //expansion rate
  
  //Calculate all pressures
  double PG = getPGamma(T9), PNu = getPNu(T9), Pe = getPe(phie,T9),
    Pb = getPb(h,T9,Y), pressure = PG+Pe+Pb;//neutrinos are decoupled

  //Calculate charge abundnace and time derivative of charge abundance
  int i;
  double S=0, dSdt=0, sumdYdt=0;
  for(i=0;i<nnuc;i++){
    S+=Z[i]*Y[i];
    dSdt+=Z[i]*dYdt[i];
    sumdYdt+=dYdt[i];
  }

  //Calculate the left and right side of charge balance equation
  double M = getM(h,T9,S,dSdt);//left side
  double N = getN(phie,T9);//right side

  //Calculate the critical derivatives
  double Tdot = dTdt(H,rho,pressure,rhoe,S,dSdt,N,M,phie,h,T9,Y,dYdt);
  double hdot = dhdt(Tdot,H,h,T9);
  double phiedot = dphiedt(phie,Tdot,H,N,M,T9,S,dSdt);
  dvdt[0]=Tdot;dvdt[1]=hdot;dvdt[2]=phiedot;
  return 0;
}

double dTdt(double H,double rho,double pres,double rhoe,
	    double S,double dSdt,double N,double M,
	    double phie,double h,double T9,
	    double *Y,double *dYdt){
  return 3*H/drdT(H,rho,pres,rhoe,S,dSdt,N,M,phie,h,T9,Y,dYdt);
}
double dhdt(double Tdot,double H,double h,double T9){
  return -3*h*(H+Tdot/T9);
}
double dphiedt(double phie,double Tdot,double H,double N,double M,
	       double T9,double S,double dSdt){
  return dphiedT(phie,N,M,T9,S)*Tdot+
    dphiedr(phie,N,M,T9,S)*3*H+
    dphiedS(phie,N,M,T9,S)*dSdt;
}

double drdT(double H,double rho,double pres,double rhoe,
	    double S,double dSdt,double N,double M,
	    double phie,double h,double T9,
	    double *Y,double *dYdt){
  return -1*(drhoGdT(T9)+drhoedT(rhoe,phie,T9)+drhobdT(h,T9,Y))/
    (rho+pres+1./3./H*(drhobdt(h,T9,dYdt)+
		       drhoedt(rhoe,phie,H,h,T9,S,dSdt,N,M)));
}
