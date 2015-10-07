#ifndef DENSITIES_H
#define DENSITIES_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "bessel.h"
#include "expansion_rate.h"
#include "physical_constants.h"
#include "program_constants.h"

/*Densities*/
double getrhoGamma(double T9);//photon energy density
double getrhoNu(double h,double T9);//neutrino energy density
double getrhoe(double phie,double T9);//electron + positron energy density
double getrhob(double h,double T9,double *Y);//baryon energy density

/*Pressures*/
double getPGamma(double T9);
double getPNu(double T9);
double getPe(double phie,double T9);
double getPb(double h, double T9, double *Y);

/*Density Derivatives*/
double drhoedt(double rhoe,double phie,double H,double h,double T9,double S,double dSdt,double N,double M);//lepton en. dens. dt
double drhobdt(double h, double T9, double *dYdt);//baryon en. dens. dt
double drhoGdT(double T9);//photon en dens. dT
double drhoNdT(double T9);//neutrino en dens dT
double drhoedT(double rhoe, double phie, double T9);//lepton en. dens dT
double drhobdT(double h, double T9, double *Y);//baryon en dens dT
double drhoedphie(double phie,double T9);//lepton en dens d(chem potential)

#endif
