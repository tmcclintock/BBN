#include "driver.h"

/* This drives the RK4 calls, beginning at T9_initial and ending when
   it reaches T9_final.
 */

/*These are global arrays and must be redeclared here*/
double deltaM[totalnnuc];
double Z[totalnnuc];
double Qvals[totalnreac];
double A[totalnnuc];
double nni[11];
double nnj[11];
double nnk[11];
double nnl[11];
double reaction_details[totalnreac][8];

/*This is the initial baryon to photon ratio*/
double eta0;


int driver(int to_print,int print_increment,double eta0_guess,double*Y_out){
  eta0 = eta0_guess;

  FILE*abundances_file = fopen("./output_files/abundances.txt","w");
  FILE*dynamics_file = fopen("./output_files/dynamics.txt","w");

  fprintf(abundances_file,"n\tp\tD\tT\tHe3\tHe4\tBe7\tLi7\tLi6\n");
  fprintf(dynamics_file,"t\tT9\th\tphie\n");

  double T9,h,phie,t;
  double Y[totalnnuc];
  double dt,dt_prev=1e-4;

  int inc=50,ip=50;//variables used for error correction in the abundance step

  /*Initialize calculation constants*/
  initialize_abundances_and_dynamics(&T9,&h,&phie,&t,Y);
  initialize_constants(deltaM,Z,Qvals,A,nni,nnj,nnk,nnl,reaction_details);  

  int fail=0;

  int n = 3 + totalnnuc;
  double vars[n],vars_new[n];
  int i,j;
  vars[0]=T9,vars[1]=h,vars[2]=phie;
  vars_new[0]=T9,vars_new[1]=h,vars_new[2]=phie;
  for(i=3;i<n;i++){
    vars[i]    =Y[i-3];
    vars_new[i]=Y[i-3];
  }

  i=0;
  while(T9>T9Final){
    fail = rk4(vars,vars_new,&t,&dt,&dt_prev,inc,ip);
    if(fail){printf("Error from driver!\n");return fail;}
    vars[0]=vars_new[0],vars[1]=vars_new[1],vars[2]=vars_new[2];
    for(j=3;j<n;j++)
      vars[j]=vars_new[j];
    t+=dt;
    dt_prev = dt;
    T9=vars[0],h=vars[1],phie=vars[2];
    if(to_print){
      if(i%print_increment==0){
	fprintf(dynamics_file,"%e\t%e\t%e\t%e\n",t,T9,h,phie);
	for(j=3;j<n;j++)
	  fprintf(abundances_file,"%e\t",vars[j]);
	fprintf(abundances_file,"\n");
      }
    }
    i++;
  }
  
  for(i=0;i<totalnnuc;i++)
    Y_out[i] = vars[i+3];

  fclose(abundances_file);
  fclose(dynamics_file);
  return fail;
}
/*
  int main(){
  int to_print = 1; //1 is true
  int print_increment = 1;
  double eta0_guess = 6.19e-10; //will be the MCMC variable
  
  int fail = 0;
  fail = driver(to_print,print_increment,eta0_guess);
  printf("main works from here!\n");
  }
*/
