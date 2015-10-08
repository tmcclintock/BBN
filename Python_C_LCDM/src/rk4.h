#ifndef RK4_H
#define RK4_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "collect_derivs.h"
#include "physical_constants.h"
#include "program_constants.h"
#include "timestep.h"

int rk4(double*,double*,double*,double*,double*,int,int);

#endif
