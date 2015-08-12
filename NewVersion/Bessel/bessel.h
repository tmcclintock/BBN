#ifndef BESS_H
#define BESS_H
#include <math.h>

double I0(double z);
double I1(double z);
double IN(int n,double z);
double K0(double z);
double K1(double z);
double KN(int n,double z);
double LBess(double z);
double MBess(double z);
double NBess(double z);//I might not need this... in the paper it is used for the derivative of LBess, which I am just doing based on the KNs
#endif
