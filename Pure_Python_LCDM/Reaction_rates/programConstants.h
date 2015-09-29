//location: "../Constants/"

#ifndef PROGCONSTANTS_H
#define PROGCONSTANTS_H

#define itNum 5//number of times bessel functions are calculated in loops
/*number of nuclear species; 
  note these are just the ones I am currently testing
*/
#define nnuc      9//number of nuclei implemented in the simulation
#define totalnnuc 9
/*total number of nuclear species I am considering
  I can use this along with nnuc to cut out higher order
  nuclides
*/

//These are now static variables
#define nreac  34//number of forward reactions
#define totalnreac 34//total number of reactions

//These are global aliases for the constants, so that they can be
//accessed in python.
int itNum_alias = itNum;
int nnuc_alias = nnuc;
int totalnnuc_alias = totalnnuc;
int nreac_alias = nreac;
int totalnreac_alias =totalnreac;


#endif
