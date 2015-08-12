//location "../InitialConditions/"

#ifndef INITALL_H
#define INITALL_H

#include "../Constants/programConstants.h"

void initializeAll(double *T9,double *h,double *phie,double *t,double *Y);
void initializeConstants(double *deltaM,double *Z,double *Qvals,double *A,double *nni,double *nnj,double *nnk,double *nnl,double reactionDetails[totalnreac][8]);

#endif
