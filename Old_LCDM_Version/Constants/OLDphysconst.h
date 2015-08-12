///FOR BBN CALC

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define Q   1.607e26 //3*c^2/(8*pi*G); J s^2 m^-3
#define aR  7.5657e-16 //radiation constant; J m^-3 K^-4
#define H0  2.27e-18//Hubble Constant; s^-1
#define pcr Q*H0*H0 //critical density
#define Tnow  2.7 //temperature now; Kelvin

//DENSITIES SEEN TODAY
#define pm1 .27*pcr // matter 
#define pb1 .04*pcr // baryonic matter
#define pd1 .73*pcr // dark matter
#define pr1 aR*Tnow*Tnow*Tnow*Tnow //radiation


//INITIAL CONDITIONS
#define t0    100
#define T0    1.e10 // 10 GK
#define a0    2.70385e-10
#define pm0   1.136e19 // J/m^3
#define pb0   pm0*.148148148 // J/m^3
#define pr0   7.5716751534e24// J/m^3
#define pd0   .73*pcr//from calculated values above
#define p0    pm0+pr0+pd0
#define pres0 pr0/3. - pd0
#define z0    3.703175e9 //get these from the dynamics simulation
#define aD0   2.693e-11 // da/dt average

#define kb 8.6173324e4 //Boltzmann constant in ev/K
#define kb2 .086173324 //Boltzmann constant in MeV/GK
#define kb3 8.6173324e-11 //Boltzmann constant in MeV/K
#define QM  1.293 //neutron mass - proton mass * c^2 in MeV
#define tau 885.7 //neutron lifetime sec
#define q  2.531  //integration variable for np rates

#define avog   6.022e23 //avogadro's number
#define avogSQ avog*avog //avog squared. Used for 2 incoming particles; TOOBIG
#define mbc2 1.504314e-10 //average mass of a baryon in Joules
#define c    2.99792458e8 //speed of light in m/s


//old constants
#define T3 6.51e6 //T at 3 minutes in Ke9, in LCDM
#define ts 55000000
#define tmax ts
//#define z 5.929862032115561
