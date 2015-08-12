#include "abundanceChanges.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../Constants/physicalConstants.h"
#include "../Constants/programConstants.h"
#include "../ReactionRates/reactionRates.h"

void getdYdt(double h,double T9,double *Y,double *dYdt){


  /*Declare array containing reaction rates*/
  double rates[totalnreac];

  /*Calculate reaction rates*/
  reactionRates(T9,h,Y,rates);

  /*Calculate Abundance derivatives*/
  /* p reactions include
     n -> p + e + v
     p + e -> n + v
     p + n -> d + g
     d + g -> n + p
     d + d -> p + t
     t + p -> d + d
     d + p -> g + 3He
     3He + g -> p + d
     3He + n -> p + t
     t + p -> n + 3He
     3He + d -> p + 4He
     4He + p -> d + 3He
     7Be + n -> p + 7Li
     7Li + p -> n + 7Be
     7Li + p -> 4He + 4He
     4He + 4He -> p + 7Li
  */
  dYdt[0]=
    -Y[0]*rates[0]
    +Y[1]*rates[1]
    -Y[0]*Y[1]*rates[2]
    +Y[2]*rates[3]
    +Y[2]*Y[2]*rates[4]/2
    -Y[3]*Y[0]*rates[5]
    -Y[2]*Y[0]*rates[6]
    +Y[4]*rates[7]
    +Y[4]*Y[1]*rates[10]
    -Y[3]*Y[0]*rates[11]
    +Y[4]*Y[2]*rates[14]
    -Y[5]*Y[0]*rates[15]
    +Y[6]*Y[1]*rates[20]
    -Y[7]*Y[0]*rates[21]
    -Y[7]*Y[0]*rates[22]
    +Y[5]*Y[5]*rates[23]/2;
  printf("dYndt parts: +Y0*r0=%e\t-Y1*r1=%e\t-Y0*Y1*r2=%e\t+Y2*r3=%e\n",
	 /*-Y[0]**/rates[0],/*+Y[1]**/rates[1],/*-Y[0]*Y[1]**/rates[2],/*+Y[2]**/rates[3]);
  ////////////////////////////////////////////////////////////

  /* n reactions include
     n -> p + e + v
     p + e -> n + v
     p + n -> d + g
     d + g -> n + p
     d + d -> n + 3He
     3He + n -> d + d
     3He + n -> p + t
     t + p -> n + 3He
     7Be + n -> p + 7Li
     7Li + p -> n + 7Be
  */
  dYdt[1]=//neutrons
    +Y[0]*rates[0]
    -Y[1]*rates[1]
    -Y[0]*Y[1]*rates[2]
    +Y[2]*rates[3]
    +Y[2]*Y[2]*rates[8]/2
    -Y[4]*Y[1]*rates[9]
    -Y[4]*Y[1]*rates[10]
    +Y[3]*Y[0]*rates[11]
    -Y[6]*Y[1]*rates[20]
    +Y[7]*Y[0]*rates[21];
  ////////////////////////////////////////////////////////////

  /* d reactions include
     p + n -> d + g
     d + g -> n + p
     d + d -> p + t
     t + p -> d + d
     d + p -> g + 3He
     3He + g -> p + d
     d + d -> n + 3He
     3He + n -> d + d
     3He + d -> p + 4He
     4He + p -> d + 3He
  */
  dYdt[2]=
    +Y[0]*Y[1]*rates[2]
    -Y[2]*rates[3]
    -Y[2]*Y[2]*rates[4]
    +Y[3]*Y[0]*rates[5]
    -Y[2]*Y[0]*rates[6]
    +Y[4]*rates[7]
    -Y[2]*Y[2]*rates[8]
    +Y[4]*Y[1]*rates[9]
    -Y[4]*Y[2]*rates[14]
    +Y[5]*Y[0]*rates[15];
  ////////////////////////////////////////////////////////////

  /* t reactions include
     d + d -> p + t
     t + p -> d + d
     3He + n -> p + t
     t + p -> n + 3He
     t + 4He -> g + 7Li
     7Li + g -> 4He + t
  */
  dYdt[3]=
    +Y[2]*Y[2]*rates[4]/2
    -Y[3]*Y[0]*rates[5]
    +Y[4]*Y[1]*rates[10]
    -Y[3]*Y[0]*rates[11]
    -Y[3]*Y[5]*rates[18]
    +Y[7]*rates[19];
  ////////////////////////////////////////////////////////////

  /* 3He reactions include
     d + p -> g + 3He
     3He + g -> p + d
     d + d -> n + 3He
     3He + n -> d + d
     3He + n -> p + t
     t + p -> n + 3He
     3He + d -> p + 4He
     4He + p -> d + 3He
     3He + 4He -> g + 7Be
     7Be + g -> 4He + 3He
   */
  dYdt[4]=
    +Y[2]*Y[0]*rates[6]
    -Y[4]*rates[7]
    +Y[2]*Y[2]*rates[8]/2
    -Y[4]*Y[1]*rates[9]
    -Y[4]*Y[1]*rates[10]
    +Y[3]*Y[0]*rates[11]
    -Y[4]*Y[2]*rates[14]
    +Y[5]*Y[0]*rates[15]
    -Y[4]*Y[5]*rates[16]
    +Y[6]*rates[17];
  ////////////////////////////////////////////////////////////

  /* 4He reactions include
     t + d -> n + 4He
     4He + n -> d + t
     3He + d -> p + 4He
     4He + p -> d + 3He
     3He + 4He -> g + 7Be
     7Be + g -> 4He + 3He
     t + 4He -> g + 7Li
     7Li + g -> 4He + t
     7Li + p -> 4He + 4He
     4He + 4He -> p + 7Li
   */
  dYdt[5]=
    +Y[3]*Y[2]*rates[12]
    -Y[5]*Y[1]*rates[13]
    +Y[4]*Y[2]*rates[14]
    -Y[5]*Y[0]*rates[15]
    -Y[4]*Y[5]*rates[16]
    +Y[6]*rates[17]
    -Y[3]*Y[5]*rates[18]
    +Y[7]*rates[19]
    +Y[7]*Y[0]*rates[22]
    -Y[5]*Y[5]*rates[23];
  ////////////////////////////////////////////////////////////

  /* 7Be reactions include
     3He + 4He -> g + 7Be
     7Be + g -> 4He + 3He
     7Be + n -> p + 7Li
     7Li + p -> n + 7Be
  */
  dYdt[6]=
    +Y[4]*Y[5]*rates[16]
    -Y[6]*rates[17]
    -Y[6]*Y[1]*rates[20]
    +Y[7]*Y[0]*rates[21];
  ////////////////////////////////////////////////////////////

  /* 7Li reactions include
     t + 4He -> g + 7Li
     7Li + g -> 4He + t
     7Be + n -> p + 7Li
     7Li + p -> n + 7Be
     7Li + p -> 4He + 4He
     4He + 4He -> p + 7Li
   */
  dYdt[7]=
    +Y[3]*Y[5]*rates[18]
    -Y[7]*rates[19]
    +Y[6]*Y[1]*rates[20]
    -Y[7]*Y[0]*rates[21]
    -Y[7]*Y[0]*rates[22]
    +Y[5]*Y[5]*rates[23]/2;
  ////////////////////////////////////////////////////////////

/*   int i;printf("\n"); */
/*   for(i=0;i<nnuc;i++){ */
/*     printf("Y[%i]=%e\tdYdt[%i]=%e\t",i,Y[i],i,dYdt[i]); */
/*   } */
/*   printf("\n"); */

}
