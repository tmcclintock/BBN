#include"criticalDerivs.h"

#include<math.h>
#include<stdio.h>
#include "../Constants/physicalConstants.h"
#include "../Constants/programConstants.h"
#include "../Densities/densities.h"
#include "../Densities/densDerivs.h"
#include "../PartialDerivatives/partials.h"

double getExpansionRate(double phie,double h,double T9,double *Y){
  double rhoG = getrhoGamma(T9);
  double rhoNu = getrhoNu(h,T9);
  double rhoe = getrhoe(phie,T9);
  double rhob = getrhob(h,T9,Y);
  double rhoTotal = rhoG+rhoe+rhob+rhoNu;
  double H = getH(rhoTotal);//and the expansion rate
  return H;
}

double getDerivs(double phie, double h, double T9,double *Y, double *dYdt, double *v){
  /*Calculate all densities and prssures here*/
  double rhoG = getrhoGamma(T9);
  double rhoNu = getrhoNu(h,T9);
  double rhoe = getrhoe(phie,T9);
  double rhob = getrhob(h,T9,Y);
  double rhoTotal = rhoG+rhoe+rhob+rhoNu;
  double rho = rhoG+rhoe;//neutrinos have decoupled
  double H = getH(rhoTotal);//and the expansion rate


  double PG = getPGamma(T9);
  //double PNu = getPNu(T9);//unused
  double Pe = getPe(phie,T9);
  double Pb = getPb(h,T9,Y);
  double pres = PG+Pe+Pb;//neutrinos have decoupled

  /*Calculate charge abundance and time derivative of charge abundance here*/
  int i;
  double S=0, dSdt=0,sumdYdt=0;
  for(i=0;i<nnuc;i++){
    S+=Z[i]*Y[i];
    dSdt+=Z[i]*dYdt[i];
    sumdYdt+=dYdt[i];
/*     if(T9>99.99) */
/*       printf("\ndYdt[%i]=%4.3e\t",i,dYdt[i]); */
  }
  
  /*Calculate left and right side of charge balance equation (eq. 31)*/
  double M=getM(h,T9,S,dSdt);//left side
  double N=getN(phie,T9);//right side

  /*Calculate critical derivatives*/
  double Tdot = dTdt(H,rho,pres,rhoe,S,dSdt,N,M,phie,h,T9,Y,dYdt);
  double hdot = dhdt(Tdot,H,h,T9);
  double phiedot = dphiedt(phie,Tdot,H,N,M,T9,S,dSdt);
  
  if(T9>99.99){
/*      printf("\nT9=%4.3e\trhob=%4.3e\tphie=%4.3e\t",T9,rhob,phie);  */
/*   printf("T9=%e\th=%e\tphie=%e\n",T9,h,phie); */
/*     printf("rG=%4.3e\trNu=%4.3e\trE=%4.3e\trb=%4.3e\tH=%4.3e\tsdYdt=%4.3e\n",rhoG,rhoNu,rhoe,rhob,H,sumdYdt); */
/*   printf("PG=%e\tPE=%e\tPb=%e\n",PG,Pe,Pb); */
/*   printf("S=%e\tdSdt=%e\tN=%e\tM=%e\n",S,dSdt,N,M); */
/*   printf("Tdot=%e\thdot=%e\tphiedot=%e\n",Tdot,hdot,phiedot); */
  }

  v[0]=Tdot;v[1]=hdot;v[2]=phiedot;
  return H;
}

double dTdt(double H,double rho,double pres,double rhoe,double S,double dSdt,double N,double M,double phie,double h,double T9,double *Y,double *dYdt){
  //double drrrdT =drdT(H,rho,pres,rhoe,S,dSdt,N,M,phie,h,T9,Y,dYdt);
  //printf("3*H=%e\tdrdT=%e\n",3.*H,drrrdT);
  return 3*H/drdT(H,rho,pres,rhoe,S,dSdt,N,M,phie,h,T9,Y,dYdt);
}
double dhdt(double Tdot,double H,double h,double T9){
  //  double left = H;
  //  double right = Tdot/T9;
  //  printf("H=%e\tTdot/T9=%e\n",left,right);
  return -3*h*(H+Tdot/T9);
}
double dphiedt(double phie,double Tdot,double H,double N,double M,double T9,double S,double dSdt){
  return dphiedT(phie,N,M,T9,S)*Tdot+
    dphiedr(phie,N,M,T9,S)*3*H+
    dphiedS(phie,N,M,T9,S)*dSdt;
}

double drdT(double H,double rho,double pres,double rhoe,double S,double dSdt,double N,double M,double phie,double h,double T9,double *Y,double *dYdt){
  return -1*(drhoGdT(T9)+drhoedT(rhoe,phie,T9)+drhobdT(h,T9,Y))/
    (rho+pres+1./3./H*(drhobdt(h,T9,dYdt)+drhoedt(rhoe,phie,H,h,T9,S,dSdt,N,M)));
}
