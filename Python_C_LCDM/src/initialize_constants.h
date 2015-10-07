#ifndef INITCONST_H
#define INITCONST_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "physical_constants.h"
#include "program_constants.h"
//#include "initialAbundances.h"
//#include "initialDynamics.h"


//void initializeAll(double *T9,double *h,double *phie,double *t,double *Y);
void initialize_constants(double *deltaM,double *Z,double *Qvals,
			  double *A,double *nni,double *nnj,double *nnk,
			  double *nnl,double reaction_details[totalnreac][8]);

#endif
