#include "reaction_rates.h"

void reaction_rates(double T9, double*f, double*r){
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
  
  if(z<5.10998997931){/* This is set by the freeze out temperature */
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
  if( T9 >=1.5)/* Smith and Kawano 1993 */
    f[11] = 47420*(1-0.850*T912+0.49*T9-0.0962*T932+.00847*T92-0.00028*T952);
  else 
/*     f[11]= 44216.0*(1+3.75191*T9+1.92934*T92+0.746503*T93+0.0197023*T94+3.00491e-6*T95)/ */
/*       (1+5.4678*T9+5.62395*T92+0.489312*T93+0.00747806*T94); /\* Ando 2006 *\/ */
    f[11] = 44060*(1+0.106597*T912-2.75037*T9+4.62949*T932-3.52204*T92
		   +1.34596*T952-0.209351*T93);/* Serpico 2004 */
/*   if(T9<1.5) f[11]=44060.*(1.-2.7503695564153143*T9-3.5220409333117897*T9*T9-0.2093513619089196*T9*T9*T9+0.10659679579058313*sqrt(T9)+ 4.62948586627009*pow(T9,3./2.) +1.3459574632779876*pow(T9,5./2.)); */
/*   else */
/*     f[11]=(1.-sqrt(T9)*0.8504+T9*0.4895-pow(T9,1.5)*0.09623+T9*0.008471*T9-T9*2.8e-4*pow(T9,1.5))*47420.;/\* abbn code *\/ */


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
/*   if (T9<4.) f[12]=((-15.709674-721.9142*T9+173.23945*T9*T9-77.371692*T9*T9*T9+126.82087e0*pow(T9,1./3.)-206.50853*pow(T9,2./3.)+2120.7339*pow(T9,4./3.)-369.61306*pow(T9,5./3.)+127.8378*pow(T9,7./3.)+100.68769*pow(T9,8./3.))*pow(T9,-2./3.))/exp(1.29042942e0*pow(T9,-1./3.)); */
/*   else f[12]=2049.72;/\* abbn code*\/ */

  /* Reaction number 12 */
  //////////////////////////////////////////////////////////////////////

  /* d + d -> n + 3He */
  if(T9<10&&T9>0.001)
    f[13]=4.67e8/T923*exp(-4.259/T913)*
      (1+1.079*T9-0.1124*T92+0.00568*T93);
  else f[13]=4.98099e7;;/* Angulo 1999 */
/*   if(T9<4.) f[13]=((-1.8436156e6 - 6.1150115e7*T9 - 2.7251853e7*T9*T9 - 2.2800422e6*T9*T9*T9 - 252433.58*pow(T9,4.) - 284357.41*pow(T9,10./3.) + 906146.25*pow(T9,11./3.) + 1.2270083e7*pow(T9,1./3.) - 1.3680884e7*pow(T9,2./3.) + 1.328894e8*pow(T9,4./3.) - 1.1916242e7*pow(T9,5./3.) + 8.3705218e6*pow(T9,7./3.) + 2.2357751e6*pow(T9,8./3.))*pow(T9,-2./3.))/exp(1.*pow(T9,-1./3.)); */
/*   else f[13]=4.98099e7;/\* abbn code *\/ */

  /* Reaction Number 13 */
  //////////////////////////////////////////////////////////////////////


  /* d + d -> p + t */
  f[14]=4.66e8/T923*exp(-4.259/T913)*(1+0.759*T9-0.0612*T92+0.00278*T93);/* Angulo 1999 */
/*   if(T9<4.) f[14]=((-5.8523126e6+2.3222535e8*T9-9.877862e6*T9*T9+5.2331507e7*pow(T9,1./3.)-1.7022642e8*pow(T9,2./3.)-1.1875268e8*pow(T9,4./3.)+5.2922232e7*pow(T9,5./3.))*pow(T9,-2./3.))/exp(1.0676573*pow(T9,-1./3.)); */
/*   else f[14]=4.021e7;/\* abbn code *\/ */

  /*Reaction Number 14 */
  //////////////////////////////////////////////////////////////////////


  /* 3He + n -> p + t */
  if(T9 > 3)
    f[15]=7.21e8*(1-0.508*T912+0.228*T9);/* Smith 1993 */
  else 
    f[15] = 6.505e8*(1-0.655*T9+0.445*T92-0.082*T93);/* Descouvemont 2004 */
/*   if(T9<2.5) f[15]=7.064935e8+6.733213571736319e8*T9+1.7181155480346258e9*T9*T9-4.5367658146835446e8*T9*T9*T9-1.2216728981712557e8*pow(T9,4.)-4.92736677238425e8*sqrt(T9)-1.3659670893994067e9*pow(T9,3./2.)-6.629932739639357e8*pow(T9,5./2.)+4.834951929033479e8*pow(T9,7./2.); */
/*   else f[15]=4.81732e8;/\* abbn code *\/ */

  /* Reaction number 15 */
  //////////////////////////////////////////////////////////////////////

  /* t + d -> n + 4He */
  f[16]=8.29e10/T923*exp(-4.524/T913-T92/0.08/0.08)*
    (1+17.2*T9+175*T92)+ 8.12e8*pow(T9,-0.712)*exp(-0.506/T9);/* Angulo 1999 */
/*   if(T9<2.5) f[16]=6.2265733e8/(exp(0.49711597/T9)*pow(T9,0.56785403)) + exp(-0.23309803*T9*T9 - 1.342742*pow(T9,-1./3.))*(-8.1144927e7 + 2.2315324e9*T9 - 2.9439669e9*T9*T9 + 1.8764462e9*T9*T9 - 6.0511612e8*pow(T9,4.) + 9.5196576e7*pow(T9,5.) - 5.2901086e6*pow(T9,6.))*pow(T9,-2./3.); */
/*   else f[16]=3.40249e8;/\* abbn code *\/ */

  /* Reaction Number 16 */
  //////////////////////////////////////////////////////////////////////

  /* 3He + d -> p + 4He */
  if(T9>2)
    f[17] = 5.021e10/T923*exp(-7.144/T913-T92/0.27/0.27)*
      (1 + 0.058*T913+0.603*T923+0.245*T9+6.97*T943+7.19*T953)
      +5.212e8/T912*exp(-1.762/T9);/* Smith 1993 */
  else 
    f[17] = 5.477e10*exp(-7.1820/T913)/T923*(1+4.367*T9-4.329*T92+1.115*T93);/* Descouvemont 2004 */
/*   if(T9<2.5) f[17]=3.1038385e8/(exp(1.6190981/T9)*pow(T9,0.12159455))+exp(-0.0062340825*T9*T9-1.4540617*pow(T9,-1./3.))*(-3.1335916e7-6.2051071e8*T9-1.8782248e9*T9*T9+6.5642773e8*T9*T9*T9+1.530887e8*pow(T9,4.)-4.9542138e8*pow(T9,10./3.)-1.770285e8*pow(T9,11./3.)+1.14185e8*pow(T9,1./3.)-2.516526e7*pow(T9,13./3.)+1.7500204e8*pow(T9,2./3.)-1.7513362e9*pow(T9,4./3.)+5.2792247e9*pow(T9,5./3.)-3.32382e9*pow(T9,7./3.)+2.0346284e9*pow(T9,8./3.))*pow(T9,-2./3.); */
/*   else f[17]=1.55167e8;/\* abbn code *\/ */

  /* Reaction number 17 */
  //////////////////////////////////////////////////////////////////////

  /* 3He + 4He -> g + 7Be */
  f[18] = 5.46e6/T923*exp(-12.827/T913)*(1-0.307*T9+0.0881*T92-0.0106*T93+0.000446*T94);/* Angulo 1999 */
/*   if(T9<2.5) f[18]=((0.000046165644-0.00046036111*T9-0.021600946*T9*T9+0.069627779*T9*T9*T9+7.346612*pow(T9,4.)-95.123199*pow(T9,5.)+391.13123*pow(T9,6.)-187.23717*pow(T9,7.)+86.111544*pow(T9,8.)-21.630169*pow(T9,9.)+3.6006922*pow(T9,10.)-0.34322836*pow(T9,11.)+0.018106742*pow(T9,12.)-0.00035681506*pow(T9,13.))*pow(T9,-1./2.))/(exp(0.48102949*T9)*pow(1.+1.17917554*T9,3.)); */
/*   else f[18]=149.06;/\* abbn code *\/ */

  /* Reaction Number 18 */
  //////////////////////////////////////////////////////////////////////

  /* t + 4He -> g + 7Li */
  f[19] = 8.2e5/T923*exp(-8.081/T913)*(1-0.389*T9+0.134*T92-0.0181*T93+0.000923*T94);/* Angulo 1999 */
/*   if (T9<2.5) f[19]=((0.094614248-4.9273133*T9+99.358965*T9*T9-989.81236*T9*T9*T9+4368.45*pow(T9,4.)+931.93597*pow(T9,5.)-391.07855*pow(T9,6.)+159.23101*pow(T9,7.)-34.407594*pow(T9,8.)+3.3919004*pow(T9,9.)+0.017556217*pow(T9,10.)-0.036253427*pow(T9,11.)+0.0031118827*pow(T9,12.)-0.00008714468*pow(T9,13.))*pow(T9,-1./2.))/(exp(8.4e-7*T9)*pow(1.+1.78616593*T9,3.)); */
/*   else f[19]=807.406;/\* abbn code *\/ */

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
/*   if (T9<2.5) */
/*     { */
/*       f[23]=((-8.9654123e7-2.5851582e8*T9-2.6831252e7*T9*T9+3.8691673e8*pow(T9,1./3.)+4.9721269e8*pow(T9,2./3.)+2.6444808e7*pow(T9,4./3.)-1.2946419e6*pow(T9,5./3.)-1.0941088e8*pow(T9,7./3.)+9.9899564e7*pow(T9,8./3.))*pow(T9,-2./3.))/exp(7.73389632*pow(T9,-1./3.)); */
/*       f[23]+=exp(-1.137519e0*T9*T9-8.6256687*pow(T9,-1./3.))*(3.0014189e7-1.8366119e8*T9+1.7688138e9*T9*T9-8.4772261e9*T9*T9*T9+2.0237351e10*pow(T9,4.)-1.9650068e10*pow(T9,5.)+7.9452762e8*pow(T9,6.)+1.3132468e10*pow(T9,7.)-8.209351e9*pow(T9,8.)-9.1099236e8*pow(T9,9.)+2.7814079e9*pow(T9,10.)-1.0785293e9*pow(T9,11.)+1.3993392e8*pow(T9,12.))*pow(T9,-2./3.); */
/*     } */
/*   else */
/*     { */
/*       f[23]=1.53403e6; */
/*       f[23]+=84516.7; */
/*     }/\* abbn code *\/ */

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

  /* He3 + He3 -> p + p + He4 */
  f[26] = 5.59e10/T923*exp(-12.277/T913)*
    (1-0.135*T9+0.0254*T92-0.00129*T93);/* Angulo 1999 */
  /* Reaction Number 26 */
  //////////////////////////////////////////////////////////////////////

  /* Li7 + d -> n + He4 + He4 */
  if(T9 < 2.5)
    f[27] = 1.66e11/T923*exp(-10.254/T913)
      +1.71e6/T932*exp(-3.246/T9)+1.49e10/T932*exp(-4.0894/T9)*
      (2.57e-2/T9 + 2.6314/T923-4.1929/T913 -2.1241+4.1136*T913);
  else
    f[27] = 1.37518e9;/* Boyd 1993 */
  /* Reaction Number 27 */
  //////////////////////////////////////////////////////////////////////

  /* Be7 + d -> p + He4 + He4 */
  f[28] = 1.07e12/T923*exp(-12.428/T913);/* Cauglin 1988*/
  //printf("f[28]=%e\n",f[28]/f[20]);  
  /* Reaction Number 28 */
  //////////////////////////////////////////////////////////////////////

  /* d + 4He -> p + Li6 */
  f[29] = 14.82/T923*exp(-7.435/T913)*
    (1+6.572*T9+0.076*T92+0.0248*T93)
    +82.8/T932*exp(-7.904/T9);/* Angulo 1999 */
  /* Reaction Number 29 */
  //////////////////////////////////////////////////////////////////////

  /* Li6 + n -> g + Li7 */
  f[30] = 5.1e3;/* Malaney 1989 */
  /* Reaction Number 30 */
  //////////////////////////////////////////////////////////////////////

  /* Li6 + n -> He4 + t */
  f[31] = 2.54e9/T932*exp(-2.39/T9)+(1.-pow(1./(T9*49.18+1.),3./2)*0.261)*1.68e8;
  /* Caughlan 1988; Given in the paper as He4 + t -> n + Li6 */
  /* Reaction Number 31 */
  //////////////////////////////////////////////////////////////////////

  /* Li6 + p -> g + Be7 */
  f[32] = 1.25e6/T923*exp(-8.415/T913)*
    (1-0.252*T9+0.0519*T92-0.00292*T93);/* Angulo 1999 */
  /* Reaction Number 32 */
  //////////////////////////////////////////////////////////////////////

  /* Li6 + p -> He4 + He3 */
  f[33] = 3.54e10/T923*exp(-8.415/T913)*
    (1-0.137*T9+0.0241*T92-0.00128*T93);/* Angulo 1999 */
  /* Reaction Number 33 */
  //////////////////////////////////////////////////////////////////////


  /*Cut off unused rates*/
  int i;
  for(i=nreac;i<totalnreac;i++){
    f[i]=0;
  }
  /*
   printf("T9=%e\n",T9);
   for(i=0;i<nreac;i++){
     printf("f[%i]=%e\n",i,f[i]);
   } 
   printf("\n"); 
  */
}
