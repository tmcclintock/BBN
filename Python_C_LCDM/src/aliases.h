/* These are global aliases for the constants, so that they can be
   accessed in python.
*/

#ifndef ALIASES_H
#define ALIASES_H

//From program_constants.h
#include "program_constants.h"
int itNum_alias = itNum;
int nnuc_alias = nnuc;
int totalnnuc_alias = totalnnuc;
int nreac_alias = nreac;
int totalnreac_alias =totalnreac;


//From physical_constants.h
#include "physical_constants.h"
double T0_alias = T0;
double eta0_alias = eta0;
double XNU_alias = XNU;
double T9Final_alias = T9Final;
double Ymin_alias = Ymin;
double CT_alias = CT;
double CY_alias = CY;
double TIMECHANGE_alias = TIMECHANGE;


#endif
