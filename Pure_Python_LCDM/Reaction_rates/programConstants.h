//location: "../Constants/"

#ifndef PROGCONSTANTS_H
#define PROGCONSTANTS_H

#define itNum 5//number of times bessel functions are calculated in loops
/*number of nuclear species; note these are just the ones I am currently testing*/
#define nnuc      9//number of nuclei implemented in the simulation
#define totalnnuc 9/*total number of nuclear species I am considering
		     I can use this along with nnuc to cut out higher order
		     nuclides
		   */

//These are now static variables
#define nreac_const  34//number of forward reactions
#define totalnreac_const 34//total number of reactions
int nreac = nreac_const;
int totalnreac =totalnreac_const;


#endif
