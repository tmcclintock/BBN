#include "abundance_derivs.h"

/* This calculates the derivatives of the abundances.
*/

int get_abundance_derivs(double*dYdt,double*Y,double*v,
			 int ki,int inc,int ip,double*Y0,double dt){
  //First collect the reaction rates
  double T9=v[0],h=v[1],phie=v[2];
  double H=get_H_wrapper(T9,h,phie,Y);
  double f[totalnreac],r[totalnreac];
  reaction_rates(T9,f,r);

  //Calculate the baryon energy density
  double rhob=h*T9*T9*T9;

  int fail;//the variable to be returned
  int i,j,k,l,n,i1,j1,ind;
  int ierror=-1;
  int nord,test,icnvm;
  double cni=0,cnj=0,cnk=0,cnl=0,rni,rnj,rnk,rnl;
  int type[nreac],ni[nreac],nj[nreac],nk[nreac],nl[nreac];
  double rev[nreac],q9[nreac];
  double a[nnuc][nnuc],a0[nnuc][nnuc],b[nnuc],x[nnuc],yx[nnuc],yY[nnuc];
  double cx,sum;
  double xdy, t;
  /* Note: some variables are only needed to do the error correcting, but
     this process isn't needed as much since we are doing RK4 instead of
     RK2.*/
  
  /* Pull out information for each reaction from reactionDetails*/
  for(i=0;i<nreac;i++){
    type[i] =(int)reaction_details[i][1];
    ni[i]   =(int)reaction_details[i][2];
    nj[i]   =(int)reaction_details[i][3];
    nk[i]   =(int)reaction_details[i][4];
    nl[i]   =(int)reaction_details[i][5];
    rev[i]  =reaction_details[i][6];
    q9[i]   =reaction_details[i][7];
  }

  /* Initialize the coefficient array */
  for(i=0;i<nnuc;i++) for(j=0;j<nnuc;j++) a[i][j]=0;
  
  /* Calculate the contribution of each reaction to each coefficient */
  for(n=0;n<nreac;n++){
    ind = type[n];//the type of the current reaction
    i   = ni[n];//nuclide index of i,j,k and l
    j   = nj[n];
    k   = nk[n];
    l   = nl[n];
    
    if(i<nnuc && l <nnuc){//used to discount higher nuclei
      /* Find the number of each nuclei involved */
      rni=nni[ind];
      rnj=nnj[ind];
      rnk=nnk[ind];
      rnl=nnl[ind];
      
      switch(ind)//goto each relevant reaction type
	{
	  
	case 0:{/* (1,0,0,1) type */
	  cni=f[n];
	  cnj=0.;
	  cnk=0.;
	  cnl=r[n];
	  break;}

	case 1: { /* (1,1,0,1) type */
	  r[n]=rev[n]*1e10*pow(T9,1.5)*exp(-q9[n]/T9)*f[n];
	  f[n]=rhob*f[n];
	  cni=Y[j]*f[n]/2.;
	  cnj=Y[i]*f[n]/2.;
	  cnk=0.;
	  cnl=r[n];
	  break;}

	case 2:	{ /* (1,1,1,1) type */
	  f[n]=rhob*f[n];
	  r[n]=rev[n]*exp(-q9[n]/T9)*f[n];
	  cni=Y[j]*f[n]/2.;
	  cnj=Y[i]*f[n]/2.;
	  cnk=Y[l]*r[n]/2.;
	  cnl=Y[k]*r[n]/2.;
	  break;}

	case 3:	{ /* (1,0,0,2) type */
	  cni=f[n];
	  cnj=0.;
	  cnk=0.;
	  cnl=Y[l]*r[n]/2.;
	  break;}

	case 4:	{ /* (1,1,0,2) type */
	  f[n]=rhob*f[n];
	  r[n]=rev[n]*exp(-q9[n]/T9)*f[n];
	  cni=Y[j]*f[n]/2.;
	  cnj=Y[i]*f[n]/2.;
	  cnk=0.;
	  cnl=Y[l]*r[n]/2.;
	  break;}

	case 5:	{ /* (2,0,1,1) type */
	  f[n]=rhob*f[n];
	  r[n]=rev[n]*exp(-q9[n]/T9)*f[n];
	  cni=Y[i]*f[n]/2.;
	  cnj=0.;
	  cnk=Y[l]*r[n]/2.;
	  cnl=Y[k]*r[n]/2.;
	  break;}

	case 6:	{ /* (3,0,0,1) type */
	  r[n]=rev[n]*1.e20*pow(T9,1.5)*pow(T9,1.5)*exp(-q9[n]/T9)*f[n];
	  f[n]=rhob*rhob*f[n];
	  cni=Y[i]*Y[i]*f[n]/6.;
	  cnj=0.;
	  cnk=0.;
	  cnl=r[n];
	  break;}
		
	case 7:	{ /* (2,1,0,1) type */
	  r[n]=rev[n]*1.e20*pow(T9,1.5)*pow(T9,1.5)*exp(-q9[n]/T9)*f[n];
	  f[n]=rhob*rhob*f[n];
	  cni=Y[j]*Y[i]*f[n]/3.;
	  cnj=Y[i]*Y[i]*f[n]/6.;
	  cnk=0.;
	  cnl=r[n];
	  break;}

	case 8:	{ /* (1,1,1,2) type */
	  f[n]=rhob*f[n];
	  r[n]=rev[n]*1.e-10*pow(T9,-1.5)*rhob*exp(-q9[n]/T9)*f[n];
	  cni=Y[j]*f[n]/2.;
	  cnj=Y[i]*f[n]/2.;
	  cnk=Y[l]*Y[l]*r[n]/6.;
	  cnl=Y[k]*Y[l]*r[n]/3.;
	  break;}

	case 9:	{ /* (1,1,0,3) type */
	  f[n]=rhob*f[n];
	  r[n]=rev[n]*1.e-10*pow(T9,-1.5)*rhob*exp(-q9[n]/T9)*f[n];
	  cni=Y[j]*f[n]/2.;
	  cnj=Y[i]*f[n]/2.;
	  cnk=0.;
	  cnl=Y[l]*Y[l]*r[n]/6.;
	  break;}

	case 10:{ /* (2,0,2,1) type */
	  f[n]=rhob*f[n];
	  r[n]=rev[n]*1.e-10*pow(T9,-1.5)*rhob*exp(-q9[n]/T9)*f[n];
	  cni=Y[i]*f[n]/2.;
	  cnj=0.;
	  cnk=Y[l]*Y[k]*r[n]/3.;
	  cnl=Y[k]*Y[k]*r[n]/6.;}
	}//end switch(ind)

      /* Reverse indices when adding them to the a[][] array
         this is so the highest abundances are calculated first*/
      i=nnuc-1-i;
      j=nnuc-1-j;
      k=nnuc-1-k;
      l=nnuc-1-l;

      /* Add the contributions where appropriate*/
      if(j<nnuc) a[j][i]+=rnj*cni;
      if(k<nnuc) a[k][i]-=rnk*cni;
      a[i][i]+=rni*cni;
      a[l][i]-=rnl*cni;

      if (j<nnuc){
	a[j][j]+=rnj*cnj;
	if(k<nnuc) a[k][j]-=rnk*cnj;
	a[i][j]+=rni*cnj;
	a[l][j]-=rnl*cnj;
      }
      if (k<nnuc){
	if(j<nnuc) a[j][k]-=rnj*cnk;
	a[k][k]+=rnk*cnk;
	a[i][k]-=rni*cnk;
	a[l][k]+=rnl*cnk;
      }
      if(j<nnuc) a[j][l]-=rnj*cnl;
      if(k<nnuc) a[k][l]+=rnk*cnl;
      a[i][l]-=rni*cnl;
      a[l][l]+=rnl*cnl;
    }//end if i<nnuc && l<nnuc
  }//end for n:0,nreac
  
  /* Zero out small elements of a[][]*/
  for(i=0;i<nnuc;i++){
    i1=nnuc-1-i;//a reverse index
    for(j=0;j<nnuc;j++){
      j1=nnuc-1-j;//a reverse index
      if(fabs(a[j][i])<H*3.e-5*Y0[j1]/Y0[i1]) a[j][i]=0;//neglect small terms
      else a[j][i]*=dt;//or else multiply by timestep
    }//end for j

    /*Add the identity to a[][] and prepare the right hand side, b*/
    a[i][i]+=1;
    b[i1]=Y0[i];//b is inverted w.r.t. Y0
  }//end for i

  /* Prepare the final right hand side and the inverted solution array, yx */
  for(i=0;i<nnuc;i++){
    x[i] =b[i];//the right hand side
    yx[i]=0;   //will hold the final inverted solution
  }

  /* Save the current a[][] matrix for corrections later.
     We only do this one the inc'th iteration of the rk4 routine, which stands
     for "include corrections". Note that this is only done on the
     first step of the RK4 routine and at no later steps, hence
     the ki==1.*/
  if(ki==1) icnvm=ip; else icnvm=0;
  if(icnvm==inc) for(i=0;i<nnuc;i++) for(j=0;j<nnuc;j++) a0[j][i]=a[j][i];
  
  /* Prepare the fail variable. If it changes to !=0, then something is messed up */
  fail=0;

  /* Gaussian elimination*/
  for(i=0;i<nnuc;i++){
    if(a[i][i]==0.){//we cannot pivot on a 0, or else we would divide by 0
      fail=i; return fail;
    }

    for(j=i+1;j<nnuc;j++){
      if(a[j][i]!=0){         //this saves some flops
	cx = a[j][i]/a[i][i]; //the scaling factor
	for(k=i+1;k<nnuc;k++) 
	  a[j][k]-=cx*a[i][k];//subtract the scaled coefficient across the row
	a[j][i]=cx;           //put the scaling factor in below the pivot
	x[j]-=cx*x[i];        //reduce the solution vector accordingly
      }
    }//end for j
  }//end for i

  /* Do the back substitution */
  /* Note: the error correction parts are commented out below */
  nord=0; //this is if error corrections are included
  do{     //this is if error corrections are included
    x[nnuc-1]/=a[nnuc-1][nnuc-1];//the solution for neutrons (the bottom abundance)
    yx[nnuc-1]+=x[nnuc-1];
    for(i=nnuc-2;i>=0;i--){
      sum=0;
      for(j=i+1;j<nnuc;j++) sum+=a[i][j]*x[j];//add up all previous solutions
      x[i]=(x[i]-sum)/a[i][i];
      yx[i]+=x[i];//the general solution for Y(t+h)
    }//end for i
    
    test=1;
    if(icnvm==inc){
      for(i=0;i<nnuc;i++){
	if(yx[i]!=0){
	  xdy=fabs(x[i]/yx[i]);//relative error
	  if(xdy>2e-4){
	    if(nord<1){
	      nord++;
	      for(j=0;j<nnuc;j++){
		t=0;
		for(k=0;k<nnuc;k++)t+=a0[j][k]*yx[k];
		x[j]=b[j]-t;
	      }//end for j
	      for(j=0;j<nnuc;j++) for(k=j+1;k<nnuc;k++) x[k]-=a[k][j]*x[j];
	      break;//out of the i loop
	    }else{
	      fail=-1;
	      ierror=i;
	      return fail;//the corrections weren't good enough
	    }//end if nord<1
	  }else test=0;//end if xdy>2e-4
	}else test=0;//end if yx[i]==0
      }//end for i
    } else test=0;//end if icnvm==inc
  }while(test);
  
  /* Put the non-inverted solution into yY, and find dYdt*/
  for(i=0;i<nnuc;i++){
    yY[i]=yx[nnuc-1-i];
    dYdt[i]=(yY[i]-Y0[i])/dt;
  }

  /* Determine success of the algorithm */
  if(fail!=0){
    if(fail==-1) printf("y(%d) failed to converge\n",ierror);
    if(fail>=1) printf("%d th diagonal term equals zero\n",fail);
  }

  return fail;//if 0 is returned then we are all good
}
