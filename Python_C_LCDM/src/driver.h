#ifndef DRIVER_H
#define DRIVER_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "initialize.h"
#include "rk4.h"

/* THIS IS A SPECIAL INCLUDE, IT IS ALWAYS AT TOP LEVEL, WILL MOVE UP SOON */
#include "aliases.h"

int driver(int,int,double);

#endif

