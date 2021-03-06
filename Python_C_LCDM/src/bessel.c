#include "bessel.h"

/* This contains some bessel functions.
   See Abramowitz and Stegun p 378-379 and recurrance relations.
*/


double I0(double z){//0th Modified Bessel Function of the first kind
  z=fabs(z);
  double t=z/3.75;
  if(z<3.75)
    return 1+3.5156229*t*t+3.0899424*t*t*t*t+1.2067492*pow(t,6)+
      0.2659732*pow(t,8)+.0360768*pow(t,10)+0.0045813*pow(t,13); 
  else
    return exp(z)/sqrt(z)*(0.39894228+0.01328592/t+0.00225319/t/t-
			   0.00157565/t/t/t+0.00916281/pow(t,4)-
			   0.02057706/pow(t,5)+0.02635537/pow(t,6)-
			   0.01647633/pow(t,7)+0.00392377/pow(t,8));
}

double I1(double z){//1st Modified Bessel Function of the first kind
  z=fabs(z);
  double t=z/3.75;
  if(z<3.75)
    return z*(0.5+0.87890594*t*t+0.51498869*t*t*t*t+0.15084934*pow(t,6)+
		 0.02658733*pow(t,8)+0.00301532*pow(t,10)+
		 0.00032411*pow(t,12));
  else
    return exp(z)/sqrt(z)*(0.39894228-0.03988024/t-0.00362018/t/t+
			   0.00163801/t/t/t-0.01031555/pow(t,4)+
			   0.02282967/pow(t,5)-0.02895312/pow(t,6)+
			   0.01787654/pow(t,7)-0.00420059/pow(t,8)); 
}

double IN(int n,double z){
  double i0=I0(z),i1=I1(z),in=0;
  int i;
  if(n==0)return i0;
  if(n==1)return i1;
  for(i=1;i<n;i++){
    in=i0-(2./z)*i*i1;
    i0=i1;
    i1=in;
  }
  return in;
}

double K0(double z){//0th Modified Bessel Function of the second kind
  z=fabs(z/2);
  if(z<=1)
    return -log(z)*I0(2*z)-0.57721566+0.42278420*z*z+0.23069756*pow(z,4)+
      0.03488590*pow(z,6)+0.00262698*pow(z,8)+0.00010750*pow(z,10)+
      0.00000740*pow(z,12);
  else
    return exp(-2*z)/sqrt(2*z)*(1.25331414-0.07832358/z+0.02189568/z/z-
			    0.01062446/z/z/z+0.00587872/z/z/z/z-
			    0.00251540/pow(z,5)+0.00053208/pow(z,6));
}
double K1(double z){//1st Modified Bessel Function of the second kind
  z=fabs(z/2);
  if(z<=1)
    return log(z)*I1(2*z)+1./(2*z)*
      (1.+0.15443144*z*z-0.67278579*z*z*z*z-0.18156897*pow(z,6)-
       0.01919402*pow(z,8)-0.00110404*pow(z,10)-0.00004686*pow(z,12));
  else
    return exp(-2*z)/sqrt(2*z)*
      (1.25331414+0.23498619/z-0.0365562/z/z+0.01504268/z/z/z-
       0.00780353/z/z/z/z+0.00325614/pow(z,5)-0.00038245/pow(z,6));
}

double KN(int n, double x){//Modified Bessel function of the second kindlol
  int i;
  double k0=K0(x),k1=K1(x),kn=0;
  if(n==0)return k0;
  else if(n==1)return k1;
  for(i=1;i<n;i++){
    kn = k0+(2./x)*i*k1;
    k0=k1;
    k1=kn;
  }
  return kn;
}

double LBess(double z){//L = K2/z
  return KN(2,z)/z;
}
double MBess(double z){//M = (3K_3/4+K_1/4)/z
  return 1./z*(3.*KN(3,z)/4.+K1(z)/4);
}
double NBess(double z){
  return 0;//might not need this. see bessel.h
}
