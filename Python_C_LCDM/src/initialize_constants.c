#include "initialize_constants.h"

/* This prepares the arrays that hold information to the types of reactions
   involved. It also initializes things like the arrays that hold mass 
   differences.
 */

void initialize_constants(double *deltaM,double *Z,double *Qvals,
			  double *A,double *nni,double *nnj,double *nnk,
			  double *nnl,double reaction_details[totalnreac][8]){
  /*both arrays go from 0 to nnuc but must be filled in by hand*/
  /*Mass excesses are in MeV from Ame2012*/
  deltaM[0]=8.071317144;//n
  deltaM[1]=7.288970591;//p
  deltaM[2]=13.135721741;//d
  deltaM[3]=14.9498061;//t
  deltaM[4]=14.9312155;//3He
  deltaM[5]=2.424915609;//4He
  deltaM[6]=15.769;//7Be
  deltaM[7]=14.907105;//7Li
  deltaM[8]=14.0898789;//6Li

  /*Charge numbers; 0:nnuc*/
  Z[0]=0;//n
  Z[1]=1;//p
  Z[2]=1;//d
  Z[3]=1;//t
  Z[4]=2;//He3
  Z[5]=2;//He4
  Z[6]=4;//7Be
  Z[7]=3;//7Li
  Z[8]=3;//6Li

  /*Atomic numbers of nuclides*/
  A[0]=1;//n
  A[1]=1;//p
  A[2]=2;//d
  A[3]=3;//t
  A[4]=3;//3He
  A[5]=4;//4He
  A[6]=7;//7Be
  A[7]=7;//7Li
  A[8]=6;//6Li

  /*Reaction types*/
  /*     0        1        2        3        4        5        6        7        8        9        10*/
  nni[0]=1;nni[1]=1;nni[2]=1;nni[3]=1;nni[4]=1;nni[5]=2;nni[6]=3;nni[7]=2;nni[8]=1;nni[9]=1;nni[10]=2;
  nnj[0]=0;nnj[1]=1;nnj[2]=1;nnj[3]=0;nnj[4]=1;nnj[5]=0;nnj[6]=0;nnj[7]=1;nnj[8]=1;nnj[9]=1;nnj[10]=0;
  nnk[0]=0;nnk[1]=0;nnk[2]=1;nnk[3]=0;nnk[4]=0;nnk[5]=1;nnk[6]=0;nnk[7]=0;nnk[8]=1;nnk[9]=0;nnk[10]=2;
  nnl[0]=1;nnl[1]=1;nnl[2]=1;nnl[3]=2;nnl[4]=2;nnl[5]=1;nnl[6]=1;nnl[7]=1;nnl[8]=2;nnl[9]=3;nnl[10]=1;

    double reaction_details_temp[totalnreac][8] = 
    {
      /* type: 0-10, each type has a unique (#n1,#n2,#n3,#n4) quartet */
      /* n1: incoming nuclide number */
      /* n2: incoming light nuclide number */
      /* n3: outgoing light nuclide number */
      /* n4: outgoing nuclide number */
      /* rev: reverse reaction coefficient */
      /* q9: energy release in reaction */
      
      /*   reac# type n1 n2 n3 n4  rev   q9   */
      {0.,0.,0.,-1.,-1.,1.,0.,0.},              // n -> p
      {1.,0.,3.,-1.,-1.,4.,0.,0.},		// H3 -> e- + v + He3
      {2.,3.,9.,-1.,-1.,5.,0.,0.},		// Li8 -> e- + v + 2He4
      {3.,0.,15.,-1.,-1.,16.,0.,0.},		// B12 -> e- + B + C12
      {4.,0.,20.,-1.,-1.,21.,0.,0.},		// C14 -> e- + v + N14
      {5.,3.,10.,-1.,-1.,5.,0.,0.},		// B8 -> e+ + v + 2He4
      {6.,0.,14.,-1.,-1.,13.,0.,0.},		// C11 -> e+ + v + B11
      {7.,0.,17.,-1.,-1.,16.,0.,0.},		// N12 -> e+ + v + C12
      {8.,0.,19.,-1.,-1.,18.,0.,0.},		// N13 -> e+ + v + C13
      {9.,0.,22.,-1.,-1.,21.,0.,0.},		// O14 -> e+ + v + N14
      {10.,0.,24.,-1.,-1.,23.,0.,0.},		// O15 -> e+ + v + N15
      /* Above are weak reactions. Below are strong reactions. */
      /* The next 11 reactions are the MOST important of the strong reactions. */
      {11.,1.,1.,0.,-1.,2.,0.471,25.82},        // p + n -> g + d
      {12.,1.,2.,1.,-1.,4.,1.63,63.750},	// H2 + p -> g + He3
      {13.,5.,2.,-1.,0.,4.,1.73,37.935},	// H2 + d -> n + He3
      {14.,5.,2.,-1.,1.,3.,1.73,46.798},        // H2 + d -> p + H3
      {15.,2.,4.,0.,1.,3.,1.002,8.863},	        // He3 + n -> p + H3
      {16.,2.,3.,2.,0.,5.,5.54,204.117},        // H3 + d -> n + He4
      {17.,2.,4.,2.,1.,5.,5.55,212.980},        // He3 + d -> p + He4
      {18.,1.,5.,4.,-1.,6.,1.11,18.423},        // He3 + a -> g + Be7
      {19.,1.,5.,3.,-1.,7.,1.11,28.640},        // H3 + a -> g + Li7
      {20.,2.,6.,0.,1.,7.,0.998,19.081},        // Be7 + n -> p + Li7
      {21.,4.,7.,1.,-1.,5.,4.69,201.291},       // Li7 + p -> a + He4

      /* Less important strong reactions. */
      {22.,1.,2.,0.,-1.,3.,1.63,72.62},		// H2 + n -> g + H3
      {23.,1.,4.,0.,-1.,5.,2.61,238.81},	// He3 + n -> g + He4
      {24.,4.,6.,0.,-1.,5.,4.70,220.39},        // Be7 + n -> a + He4
      {25.,1.,3.,1.,-1.,5.,2.61,229.932},       // H3 + p -> g + He4
      {26.,10.,4.,-1.,1.,5.,3.39,149.230},      // He3 + He3 -> 2p + He4
      {27.,8.,7.,2.,0.,5.,9.95,175.476},        // Li7 + d -> n + a + He4
      {28.,8.,6.,2.,1.,5.,9.97,194.557},        // Be7 + d -> p + a + He4
      {29.,1.,5.,2.,-1.,8.,1.531,17.118},       // H2 + a -> p + Li6
      {30.,1.,8.,0.,-1.,7.,1.19,84.17},         // Li6 + n -> g + Li7
      {31.,2.,8.,0.,3.,5.,1.070,55.494},        // Li6 + n -> a + H3 
      {32.,1.,8.,1.,-1.,6.,1.187,65.052},       // Li6 + p -> g + Be7
      {33.,2.,8.,1.,4.,5.,1.067,46.641}        // Li6 + p -> a + He3
      /* All reactions above here are the "most important" reactions
         according to Kawano's 1993 code.*/
    };
    
      /* ABBN functions to be added, one day
	35.,1.,7.,0.,-1.,9.,1.31,23.59,		// Li7 + n -> g + Li8
	36.,1.,12.,0.,-1.,13.,3.04,132.95,	// B10 + n -> g + B11
	37.,1.,13.,0.,-1.,15.,2.34,39.10,	// B11 + n -> g + B12
	38.,2.,14.,0.,1.,13.,1.002,32.080,	// C11 + n -> p + B11
	39.,2.,12.,0.,5.,7.,0.758,32.382,	// B10 + n -> a + Li7
	40.,1.,8.,1.,-1.,10.,1.30,1.595,	// Be7 + p -> g + B8
	41.,1.,11.,1.,-1.,12.,0.973,76.427,	// Be9 + p -> g + B10
	42.,1.,12.,1.,-1.,14.,3.03,100.840,	// B10 + p -> g + C11
	43.,1.,13.,1.,-1.,16.,7.01,185.173,	// B11 + p -> g + C12
	44.,1.,14.,1.,-1.,17.,2.33,6.975,	// C11 + p -> g + N12
	45.,2.,15.,1.,0.,16.,3.00,146.08,	// B12 + p -> n + C12
	46.,2.,11.,1.,5.,6.,0.618,24.674,	// Be9 + p -> a + Li6
	47.,2.,12.,1.,5.,8.,0.754,13.301,	// B10 + p -> a + Be7
	48.,2.,15.,1.,5.,11.,0.292,79.89,	// B12 + p -> a + Be9
	49.,1.,6.,5.,-1.,12.,1.58,51.753,	// Li6 + a -> g + B10
	50.,1.,7.,5.,-1.,13.,4.02,100.538,	// Li7 + a -> g + B11
	51.,1.,8.,5.,-1.,14.,4.02,87.539,	// Be7 + a -> g + C11
	52.,2.,10.,5.,1.,14.,3.08,86.00,	// B8 + a -> p + C11
	53.,2.,9.,5.,0.,13.,3.07,76.96,		// Li8 + a -> n + B11
	54.,2.,11.,5.,0.,16.,10.3,66.160,	// Be9 + a -> n + C12
	55.,2.,11.,2.,0.,12.,2.07,50.63,	// Be9 + d -> n + B10
	56.,2.,12.,2.,1.,13.,6.44,107.13,	// B10 + d -> p + B11
	57.,2.,13.,2.,0.,16.,14.9,159.36,	// B11 + d -> n + C12
	58.,7.,5.,0.,-1.,11.,0.584,18.260,	// He4 + a + n -> g + Be9
	59.,6.,5.,-1.,-1.,16.,2.00,84.420,	// He4 + 2a -> g + C12
	60.,8.,9.,1.,0.,5.,3.58,177.73,		// Li8 + p -> n + a + He4
	61.,8.,10.,0.,1.,5.,3.58,218.82,	// B8 + n -> p + a + He4
	62.,8.,11.,1.,2.,5.,0.807,7.555,	// Be9 + p -> d + a + He4
	63.,9.,13.,1.,-1.,5.,3.50,100.753,	// B11 + p -> 2a + Be4
	64.,9.,14.,0.,-1.,5.,3.49,132.83,	// C11 + n -> 2a + He4
	65.,1.,16.,0.,-1.,18.,0.886,57.41,	// C12 + n -> g + C13
	66.,1.,18.,0.,-1.,20.,3.58,94.88,	// C13 + n -> g + C14
	67.,1.,21.,0.,-1.,23.,2.71,125.74,	// N14 + n -> g + N15
	68.,2.,19.,0.,1.,18.,1.002,34.846,	// N13 + n -> p + C13
	69.,2.,21.,0.,1.,20.,3.003,7.263,	// N14 + n -> p + C14
	70.,2.,24.,0.,1.,23.,1.002,41.037,	// O15 + n -> p + N15
	71.,2.,24.,0.,5.,16.,0.709,98.661,	// O15 + n -> a + C12
	72.,1.,16.,1.,-1.,19.,0.884,22.553,	// C12 + p -> g + N13
	73.,1.,18.,1.,-1.,21.,1.19,87.621,	// C13 + p -> g + N14
	74.,1.,20.,1.,-1.,23.,0.900,118.452,	// C14 + p -> g + N15
	75.,1.,19.,1.,-1.,22.,3.57,53.706,	// N13 + p -> g + O14
	76.,1.,21.,1.,-1.,24.,2.70,84.678,	// N14 + p -> g + O15
	77.,2.,23.,1.,-1.,25.,3.62,140.734,	// N15 + p -> g + O16
	78.,2.,23.,1.,5.,16.,0.706,57.623,	// N15 + p -> a + C12
	79.,1.,16.,5.,-1.,25.,5.13,83.111,	// C12 + a -> g + O16
	80.,2.,12.,5.,1.,18.,9.36,47.16,	// B10 + a -> p + C13
	81.,2.,13.,5.,1.,20.,11.0,9.098,	// B11 + a -> p + C14
	82.,2.,14.,5.,1.,21.,3.68,33.915,	// C11 + a -> p + N14
	83.,3.,17.,5.,1.,24.,4.26,111.87,	// N12 + a -> p + O15
	84.,3.,19.,5.,1.,25.,5.81,60.557,	// N13 + a -> p + O16
	85.,2.,12.,5.,0.,19.,9.34,12.287,	// B10 + a -> n + N13
	86.,2.,13.,5.,0.,21.,3.67,1.835,	// B11 + a -> n + N14
	87.,2.,15.,5.,0.,23.,4.25,88.47,	// B12 + a -> n + N15
	88.,2.,18.,5.,0.,25.,5.79,25.711	// C13 + a -> n + O16
	};*/
  int i,j;
  for(i=0;i<totalnreac;i++){
    for(j=0;j<8;j++){
      reaction_details[i][j] = reaction_details_temp[i][j];
    }
  }//reaction_details now contains all of the information for the reactions
}
