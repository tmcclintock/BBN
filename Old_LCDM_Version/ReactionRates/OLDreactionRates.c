#include "reactionRates.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../Constants/physicalConstants.h"
#include "../Constants/programConstants.h"


void reactionRates(double T9, double *f, double *r){
  /*Calculate the inverse of the Temperature*/
  double z=mecc/kb/T9;

  /*Name the Boltzmann constant k, in units of MeV/GK*/
  //  double k=kbMeV;//this is used in some if statements
  
  /*Temporary and iteration variables*/
  double T9A=0;

  /*Calculate multiples of T9*/
  double T92=T9*T9,T93=T9*T9*T9;
  double T94=T9*T9*T9*T9;
  double T95=T9*T9*T9*T9*T9;
  double T913=pow(T9,1./3),T923=pow(T9,2./3);
  double T943=pow(T9,4./3),T953=pow(T9,5./3);
  //double T973=pow(T9,7./3),T983=pow(T9,8./3);
  double T912=sqrt(T9),T932=pow(T9,3./2),T952=pow(T9,5./2);
  //  double T96=pow(T9,6);

  /*Nuclear partition functions given in Angulo 1999 and Wagoner 1969*/
  double Gpart[totalnnuc];
  Gpart[0]=1; /*p*/
  Gpart[1]=1; /*n*/
  Gpart[2]=1; /*d*/
  Gpart[3]=1; /*t*/
  Gpart[4]=1; /*3He*/
  Gpart[5]=1; /*4He*/
  Gpart[6]=1+0.52*pow(T9,-0.0575)*exp(-5.59/T9+0.0127*T9); /*7Li*/
  Gpart[7]=1+0.516*pow(T9,-0.518)*exp(-5.02/T9+0.012*T9); /*7Be*/
  //    1+3.47*pow(T9,-.347)*exp(-25.9/T9+.0558*T9), /*6Li*/
  

  /* n -> p */
  int ie;
  double a[14]={0.,0.15735,0.46172e1,-0.40520e2,0.13875e3,-0.59898e2,0.66752e2,-0.16705e2,0.38071e1,-0.39140,0.23590e-1,-0.83696e-4,-0.42095e-4,0.17675e-5};
  double b[11]={0.,0.22211e2,-0.72798e2,0.11571e3,-0.11763e2,0.45521e2,-0.37973e1,0.41266e0,-0.26210e-1,0.87934e-3,-0.12016e-4};
  f[0]=1.;
  for(ie=1;ie<=13;ie++) f[0]+=a[ie]/pow(z,ie);
  f[0]*=exp(-0.33979/z)/tau; /* n->p */
  
  if(z<5.10998997931){
    r[0]=-0.62173;
    for(ie=1;ie<=10;ie++) r[0]+=b[ie]/pow(z,ie);
    r[0]*=exp(-2.8602*z)/tau; /* p->n */
  }
  else r[0]=0.; /* Serpico 2004 */
  /* Reaction number 0 */
  //////////////////////////////////////////////////////////////////////

  /* WEAK REACTIONS */
  f[1]=1.78141141239e-9;	/* H3 -> e- + v + He3 */
  f[2]=0.827;	 	/* Li8 -> e- + v + 2He4 */
  f[3]=34.3;	 	/* B12 -> e- + B + C12 */
  f[4]=3.834e-12;	 	/* C14 -> e- + v + N14 */
  f[5]=0.9;	 	/* B8 -> e+ + v + 2He4 */
  f[6]=5.668e-4;	 	/* C11 -> e+ + v + B11 */
  f[7]=63.01;	 	/* N12 -> e+ + v + C12 */
  f[8]=0.001159;	 	/* N13 -> e+ + v + C13 */
  f[9]=0.0098171;	 	/* O14 -> e+ + v + N14 */
  f[10]=0.0056704; 	/* O15 -> e+ + v + N15 */
  for(ie=1;ie<11;ie++){f[ie]=0;r[ie]=0;}
  //////////////////////////////////////////////////////////////////////


  /* p + n -> g + d */
  if(T9<=100 && T9 >=11.6)/* Smith and Kawano 1993 */
    f[11] = 47420*(1-0.850*T912+0.49*T9-0.0962*T932+.00847*T92-0.00028*T952);
  else 
/*     f[11]= 44216.0*(1+3.75191*T9+1.92934*T92+0.746503*T93+0.0197023*T94+3.00491e-6*T95)/ */
/*       (1+5.4678*T9+5.62395*T92+0.489312*T93+0.00747806*T94); /\* Ando 2006 *\/ */
    f[11] = 44060*(1+0.106597*T912-2.75037*T9+4.62949*T932-3.52204*T92
		   +1.34596*T952-0.209351*T93);/* Serpico 2004 */
  /* Reaction number 11 */
  //////////////////////////////////////////////////////////////////////

  /* d + p -> g + 3He */
/*   if(T9 > 10) */
/*     f[12] = 2.58e3/T923*exp(-3.721/T913)*(1+3.96*T9+0.116*T92);/\* Angulo 1999 *\/ */
/*   else */
/*     f[12] = exp(1.29043/T913)/T923* */
/*       (-15.7097+126.821*T913-206.509*T923-721.914*T9+2120.73*T943-369.613*T953 */
/*        +173.239*T92+127.838*T973+100.688*T983-77.3717*T93);/\* Serpico 2004 *\/ */
  if(T9>0.11)
    f[12] = 2.58e3/T923*exp(-3.721/T913)*(1+3.96*T9+0.116*T92);/* Angulo 1999 */
  else
    f[12]=1.81e3/T923*exp(-3.721/T913)*(1+14.3*T9-90.5*T92+395*T93);/* Angulo 1999 */
  /* Reaction number 12 */
  //////////////////////////////////////////////////////////////////////

  /* d + d -> n + 3He */
  if(T9<10&&T9>0.001)
    f[13]=4.67e8/T923*exp(-4.259/T913)*
      (1+1.079*T9-0.1124*T92+0.00568*T93);
  else f[13]=4.98099e7;;/* Angulo 1999 */
  /* Reaction Number 13 */
  //////////////////////////////////////////////////////////////////////


  /* d + d -> p + t */
  f[14]=4.66e8/T923*exp(-4.259/T913)*(1+0.759*T9-0.0612*T92+0.00278*T93);/* Angulo 1999 */
  /*Reaction Number 14 */
  //////////////////////////////////////////////////////////////////////


  /* 3He + n -> p + t */
  if(T9 > 3)
    f[15]=7.21e8*(1-0.508*T912+0.228*T9);/* Smith 1993 */
  else 
    f[15] = 6.505e8*(1-0.655*T9+0.445*T92-0.082*T93);/* Descouvemont 2004 */
  /* Reaction number 15 */
  //////////////////////////////////////////////////////////////////////

  /* t + d -> n + 4He */
  f[16]=8.29e10/T923*exp(-4.524/T913-T92/0.08/0.08)*
    (1+17.2*T9+175*T92)+ 8.12e8*pow(T9,-0.712)*exp(-0.506/T9);/* Angulo 1999 */
  /* Reaction Number 16 */
  //////////////////////////////////////////////////////////////////////

  /* 3He + d -> p + 4He */
  if(T9>2)
    f[17] = 5.021e10/T923*exp(-7.144/T913-T92/0.27/0.27)*
      (1 + 0.058*T913+0.603*T923+0.245*T9+6.97*T943+7.19*T953)
      +5.212e8/T912*exp(-1.762/T9);/* Smith 1993 */
  else 
    f[17] = 5.477e10*exp(-7.1820/T913)/T923*(1+4.367*T9-4.329*T92+1.115*T93);/* Descouvemont 2004 */
  /* Reaction number 17 */
  //////////////////////////////////////////////////////////////////////

  /* 3He + 4He -> g + 7Be */
  f[18] = 5.46e6/T923*exp(-12.827/T913)*(1-0.307*T9+0.0881*T92-0.0106*T93+0.000446*T94);/* Angulo 1999 */
  /* Reaction Number 18 */
  //////////////////////////////////////////////////////////////////////

  /* t + 4He -> g + 7Li */
  f[19] = 8.2e5/T923*exp(-8.081/T913)*(1-0.389*T9+0.134*T92-0.0181*T93+0.000923*T94);/* Angulo 1999 */
  /* Reaction Number 19 */
  //////////////////////////////////////////////////////////////////////

  /* 7Be + n -> p + 7Li */
  
  //if(T9 > 0.2){
  T9A=T9/(1+13.076*T9);
  f[20]=2.675e9*(1-0.56*T912+0.179*T9-0.0283*T932+0.00221*T92-6.85e-5*T952)+
    9.391e8*pow(T9A/T9,3./2)+4.467e7/T932*exp(-0.07486/T9);/* Smith 1993 */
/*   }else{ */
/*   f[20] = 4.609e9*(1-7.518*T9+53.093*T92-135.953*T93);/\* Descouvemont 2004 *\/ */
/*   } */
  /*   if(T9<2.5) f[20]=6.8423032e9+1.7674863e10*T9+2.6622006e9*T9*T9-3.3561608e8*T9*T9*T9-5.9309139e6*pow(T9,4.)-1.4987996e10*sqrt(T9)-1.0576906e10*pow(T9,3./2.)+2.7447598e8*pow(T9,5./2.)+7.6425157e7*pow(T9,7./2.)-2.282944e7*pow(T9,-3./2.)/exp(0.050351813/T9); */
/*   else f[20]=1.28039e9;//alterbbn code */
  /* Reaction Number 20 */
  //////////////////////////////////////////////////////////////////////

  /* 7Li + p -> 4He + 4He */
  f[21] = 7.2e8/T923*exp(-8.473/T913 - T92/6.5/6.5)*
    (1+1.05*T9-0.653*T92+0.185*T93-0.0212*T94+0.00093*T95)+
    9.85e6*pow(T9,0.576)*exp(-10.415/T9);/* Angulo 1999 */
    /* Reaction Number 21 */
  //////////////////////////////////////////////////////////////////////
  /*
    The reactions above this comment are the MOST important strong reactions.
    Below here are less important strong reactions.
  */

  /* H2 + n -> g + H3 */
  if(T9 > 10)
    f[22] = 66.2 + 1251*T9; /* Fowler et al. 1967 */
  else
    f[22] = 214*pow(T9,0.075) + 742*T9;/* Nagai 2006 */
  /* Reaction Number 22 */
  //////////////////////////////////////////////////////////////////////

  /* He3 + n -> g + He4 */
  f[23] = 66.2*(1+905*T9); /* Wagoner 1969 */
  /* Reaction Number 23 */
  //////////////////////////////////////////////////////////////////////

  /* Be7 + n -> a + He4 */
  f[24] = 2.5e4 * (1 + 3760*T9);/* Wagoner 1969 */
  /* Reaction Number 24 */
  //////////////////////////////////////////////////////////////////////

  /* H3 + p -> g + He4 */
  f[25] = 2.2e4/T923*exp(-3.869/T913)*
    (1+0.108*T913+1.68*T923+1.26*T9+0.551*T943+1.06*T953);/* Caughlan 1988*/
  /* Reaction Number 25 */
  //////////////////////////////////////////////////////////////////////


  /*Cut off unused rates*/
  int i;
  for(i=nreac;i<totalnreac;i++){
    f[i]=0;
  }

/*   printf("T9=%e\t",T9); */
/*   for(i=0;i<nreac;i++){ */
/*     i=i; */
/*     printf("f[%i]=%e\t",i,f[i]); */
/*   } */
/*   printf("\n"); */
}
