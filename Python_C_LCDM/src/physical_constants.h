//location: "../Constants/"

#ifndef PHYSICALCONST_H
#define PHYSICALCONST_H

#include <math.h>
#include "program_constants.h"

#define T0 100 //GK; initial temperature
#define XNU 3.28  //Effective number of neutrino species
#define T9Final 1e-1//test final temp//0.1e9 //0.1GK; final temperature

#define Ymin 1e-30 //minimum abundance
#define CT 0.01 //fractional Temperature timestep change
#define CY 0.1 //fractional Abundance timestep change
#define TIMECHANGE 1.5 //factor that the timestep can increase by

//MeV constants
#define QPNMeV 1.293 //MeV; neutron mass - proton mass
#define GMeV 1.19e-34 //cm^3/s^2/(MeV/c/c); Newton's gravitational constant
#define kbMeV 0.086173 //MeV/GK; Boltzmann constant
#define hbarMeV 6.582119e-22 //s*MeV; h_{bar}; Plank constant
#define meccMeV 0.5109989 //MeV; electron mass*c*c
#define arMeV 4.7222e27 //MeV/cm/cm/cm/GK/GK/GK/GK; stefan-boltzmann constant (aka the radiation constant)

#define QPN 2.072e-6 //g*cm*cm/s/s; neutron mass - proton mass
#define tau 880.1 //sec; neutron lifetime, not currently precise
#define G 6.67e-8 //cm*cm*cm/s/s/g; Newton's Gravitational constant in cgs
#define kb 1.38065e-7 //g*cm*cm/s/s/GK; Boltzmann constant
#define NA 6.022141e23 //number/mol; Avogadro number
#define hbar 1.054572e-27 //g*cm*cm/s; h_{bar}; Plank constant
#define c 2.998e10 //cm/s; speed of light
#define ar 8.418 //g/cm/cm/cm/GK/GK/GK/GK; radiation constant/c/c
#define zeta 1.388e-4 //1/GK; constant used for calculating baryon energy dens
#define M_u 931.454 //MeV/c/c; one atomic mass unit
#define M_ug 1.660468e-24 //grams; 1 amu in grams
#define me 9.10938e-28 //g; electron mass
#define mecc 8.187105e-7 //g*cm*cm/s/s; electron mass*c*c
#define PI M_PI

#define eta0 6.19e-10 //unitless; initial baryon to photon ratio
//extern double eta0; //Initial baryon to photon ratio, unitless

extern double rhob0; //initial baryon mass density in g/cm/cm/cm

/* Nuclide information */
extern double deltaM[totalnnuc];//MeV; mass excess for each nuclide
extern double Z[totalnnuc];//charge number of each nuclear speciesw
extern double Qvals[totalnreac];//energy differences used in reverse reactions
extern double A[totalnnuc];

/* Reaction Information */
/* There are 11 types of reactions,
   so we need to know the format of
   all 11 types.
*/
extern double nni[11];
extern double nnj[11];
extern double nnk[11];
extern double nnl[11];

/* A list of all reactions and their details */
extern double reaction_details[totalnreac][8];
#endif
