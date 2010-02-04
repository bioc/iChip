/* Codes adapted from IsingChip8.c; only keep the necessary functions in IsingChip8.c */
/* Functions in this file are used for the paper "Bayesian Modeling of ChIP-chip data through a high-order Ising model */
/*The latest version, 12/20/2008 */
/* Qianxing Mo (moq@mskcc.org),Department of Epidemiology and Biostatistics, Memorial Sloan-Kettering Cancer Center */
/*
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
*/
#include <Rmath.h>
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>


/*d heterogenous variance 
  modified by 10/08/2009, add constraint mu0<mu1 
  modified on 1/20/2010, input chr and enrichement separately as vectors */
void iChip2(int *burning,int *size,int *nrow,int *chr,double *mydt,int *halfwin,double *sdCutoff,double *beta,
	    double *postX,int *X,double *mu0,double *lambda0,double *mu1,double *lambda1,int *verbose){
  int *nprobe,*position;
  int sampleSize,LB,RB,i,j,r,nnoise,nsignal,state,nrowm1,wtsum;
  double *score2,z0,z1,p,meanScore,varScore,sdScore,mean3sd,sdnoise,sdsignal;
  double sumnoise,sumsignal,sumsigx2,sumnoix2,Scale0,Scale1;
  nrowm1 = *nrow - 1;
  sampleSize = (*burning) + (*size);
  nprobe = (int *)R_alloc(*nrow, sizeof(int));
  position = (int *)R_alloc(*nrow, sizeof(int));
  score2 = (double *)R_alloc(*nrow, sizeof(double));

  if(nprobe==NULL || score2==NULL || position==NULL){
    error("Error: Fail to allocate memory!\n");
    /*    exit(1); */
  }

  meanScore = 0;
  varScore = 0;
  /*mydt, first column is chromosome; second is the score */
  for(j=0; j<(*nrow); j++){
    meanScore = meanScore + mydt[j];
  }
    
  meanScore = meanScore/(*nrow);
  for(i=0; i<(*nrow); i++){
    varScore = varScore + R_pow_di((mydt[i]-meanScore),2);
  }
  varScore = varScore/(*nrow-1);
  sdScore = sqrt(varScore);
  mean3sd = meanScore + (*sdCutoff)*sdScore;
  /*  Rprintf("overallMean = %f, overallVar = %f  \n",meanScore,varScore); */
  
  sumsignal = 0;
  nsignal = 0;
  sumsigx2 = 0;
  sumnoix2 = 0;
  for(i=0; i<(*nrow); i++){
    score2[i] = mydt[i] * mydt[i];
    if(mydt[i] > mean3sd){
      X[i] = 1;
      sumsignal = sumsignal + mydt[i];
      sumsigx2 = sumsigx2 + score2[i];
      nsignal = nsignal + 1;
    }else {
      X[i] = -1;
      sumnoix2 = sumnoix2 + score2[i];
    }
  }
  sumnoise = (*nrow)*meanScore - sumsignal;
  nnoise = (*nrow) - nsignal;
  lambda0[0] =(nnoise-1)/(sumnoix2 - sumnoise/nnoise*sumnoise);
  lambda1[0] = (nsignal-1)/(sumsigx2 - sumsignal/nsignal*sumsignal);
  /*  Rprintf("mNoise = %f, varNoise = %f \n",sumnoise/nnoise,1.0/lambda0[0]); */
  /*  Rprintf("mSignal = %f, varSignal = %f \n",sumsignal/nsignal,1.0/lambda1[0]); */
  /*Find the start and end location for each probe */
  for(i=0; i<(*nrow); i++){
    LB = i-1;
    RB = i+1;
    r = 0;
    while((LB > -1) && (chr[LB] == chr[i]) && (r < (*halfwin))){
      LB=LB-1;
      r++;
    }
    LB = LB+1;
    r = 0;
    while((RB < (*nrow)) && (chr[RB] == chr[i]) && (r < (*halfwin))){
      RB=RB+1;
      r++;
    }
    RB = RB-1; 
    nprobe[i] = (RB-LB+1);
    position[i] = LB;     /* Start position*/
  }
  /* Gibbs sampling */
  GetRNGstate();
  mu0[0] = rnorm(sumnoise/nnoise,sqrt(1.0/lambda0[0]/nnoise));
  mu1[0] = rnorm(sumsignal/nsignal,sqrt(1.0/lambda1[0]/nsignal));

  for(r=0; r<sampleSize;r++){
    if(r >0){
      mu0[r] = rnorm(sumnoise/nnoise, sqrt(1.0/lambda0[r-1]/nnoise));
      mu1[r] = rnorm(sumsignal/nsignal, sqrt(1.0/lambda1[r-1]/nsignal));
    }
    Scale0 = 2.0/(sumnoix2-2.0*mu0[r]*sumnoise+nnoise*mu0[r]*mu0[r]);
    lambda0[r] = rgamma(nnoise/2.0,Scale0);    
    Scale1 = 2.0/(sumsigx2-2.0*mu1[r]*sumsignal+nsignal*mu1[r]*mu1[r]);
    lambda1[r] = rgamma(nsignal/2.0,Scale1);
 
    for(i=0; i<(*nrow); i++){
      state = X[i];
      wtsum = 0 - X[i];
      /* count the number of 1 around probe[i] including probe[i] itself,so initialize wtsum = -X[i] */
      for(j=0; j< nprobe[i]; j++){
	wtsum = wtsum + X[j+position[i]];
      }
      z0 = R_pow_di((mydt[i]-mu0[r]),2)*lambda0[r]/2.0;
      z1 = R_pow_di((mydt[i]-mu1[r]),2)*lambda1[r]/2.0;
      p = 1.0/(1.0 + sqrt(lambda0[r]/lambda1[r])*exp(z1 - z0 - (*beta)*wtsum*2.0));
      if(runif(0,1) < p){
	X[i] = 1;
      }else {
	X[i] = -1;
      }
      if(state != X[i]){
	/*if state change from -1 to 1 */
	if(state == -1){
	  sumnoise = sumnoise - mydt[i];
	  sumnoix2 = sumnoix2 - score2[i];
	  nnoise = nnoise - 1;
	  sumsignal = sumsignal + mydt[i];
	  sumsigx2 = sumsigx2 + score2[i];
	  nsignal = nsignal + 1;
	}else {
	  sumnoise = sumnoise + mydt[i];
	  sumnoix2 = sumnoix2 + score2[i];
	  nnoise = nnoise + 1;
	  sumsignal = sumsignal - mydt[i];
	  sumsigx2 = sumsigx2 - score2[i];
	  nsignal = nsignal - 1;
	}
      }
      /* identifiability constraint: let mu0 < mu1 for; if not, -1 is binding site*/
      if(r >= (*burning)){
	if(mu0[r] < mu1[r]){
	  if(X[i] == 1){
	    postX[i] = postX[i] + 1;
	  }
	}else {
	  if(X[i] == -1){
	    postX[i] = postX[i] + 1;
	  }
	}
      }
    }
    if(r%2000 == 0){
      R_CheckUserInterrupt();
      if(*verbose == 1){
	Rprintf("%d  ",r);
      }
    }
  }
  PutRNGstate();
  /*
  Rprintf("\n mean0 = %f, var0 = %f, ",sumnoise/nnoise,(sumnoix2 - sumnoise/nnoise*sumnoise)/(nnoise-1));
  Rprintf("mean1 = %f, var1 = %f \n",sumsignal/nsignal,(sumsigx2 - sumsignal/nsignal*sumsignal)/(nsignal-1));
  Rprintf("Note: if mean1 = -inf or inf, no enriched region is found!\n");
  */
  if(*verbose == 1){
    Rprintf("\n");
  }
  if((nnoise <= 1) || (nsignal <= 1)){
    Rprintf("Warning: all probes are in the same state at the last MCMC iteration.\n NO enriched region is found!\n");
  }

  for(i=0; i<(*nrow); i++){
    postX[i] = (double)postX[i]/(*size);
    /*   free(mydt[i]); */
  }
  
  /*
  free(mydt);
  free(nprobe);
  free(position);
  free(score2);
 */
}


/* one parameter standard ising model */
/* treat beta as random variable, use MH algorithm to estimate it */
/* the variance of binding and non-binding sites are assumed the same */
void iChip1(int *burning,int *Size,int *nrow,double *dt,double *sdCutoff1,double *betaStart,double *minBeta,
	    double *maxBeta,double *ransd,double *postX,int *X,double *postBeta,double *mu0,
	    double *mu1,double *lambda,int *verbose){

  int sampleSize,i,j,r,nnoise,nsignal,state,nrowm1,sumxx,rm1;
  double *score2,z0,z1,p,meanScore,varScore,sdScore,mean3sd,beta,Scale0;
  double halfsize,sumnoise,sumsignal,sumsigx2,sumnoix2,newBeta,ratio;
  halfsize = *nrow/2.0;
  sampleSize = (*burning) + (*Size);
  score2 = (double *)R_alloc(*nrow,sizeof(double));
  nrowm1 = *nrow - 1;
  if(score2 == NULL){
    error("Error: Fail to allocate memory!\n");
  }
  meanScore = 0;
  varScore = 0;

  for(i=0; i<(*nrow); i++){
    meanScore = meanScore + dt[i]; 
  }
  meanScore = meanScore/(*nrow);
  for(i=0; i<(*nrow); i++){
    varScore = varScore + R_pow_di((dt[i]-meanScore),2);
  }
  varScore = varScore/(*nrow-1);
  sdScore = sqrt(varScore);
  mean3sd = meanScore + (*sdCutoff1)*sdScore;
  /*  Rprintf("overallMean = %f, overallVar = %f \n",meanScore, varScore); */
  
  sumsignal = 0;
  nsignal = 0;
  sumsigx2 = 0;
  sumnoix2 = 0;
  for(i=0; i<(*nrow); i++){
    score2[i] = dt[i] * dt[i];
    if(dt[i] > mean3sd){
      X[i] = 1;
      sumsignal = sumsignal + dt[i];
      sumsigx2 = sumsigx2 + score2[i];
      nsignal = nsignal + 1;
    }else {
      X[i] = -1;
      sumnoix2 = sumnoix2 + score2[i];
    }
  }

  sumnoise = (*nrow)*meanScore - sumsignal;
  nnoise = (*nrow) - nsignal;
  mu0[0] =  sumnoise/nnoise;
  mu1[0] =  sumsignal/nsignal;
  /*
  Rprintf("mNoise = %f, varNoise = %f \n",sumnoise/nnoise,(sumnoix2 - sumnoise/nnoise*sumnoise)/(nnoise-1)); 
  Rprintf("mSignal = %f, varSignal = %f \n",sumsignal/nsignal,(sumsigx2 - sumsignal/nsignal*sumsignal)/(nsignal-1));
  */

  beta = *betaStart;

  /* Gibbs sampling */
  GetRNGstate();
  Scale0 = (sumnoix2-2.0*mu0[0]*sumnoise+nnoise*mu0[0]*mu0[0])+(sumsigx2-2.0*mu1[0]*sumsignal+nsignal*mu1[0]*mu1[0]);
  Scale0 = 2.0/Scale0;
  lambda[0] = rgamma(halfsize,Scale0);
  for(r=0; r<sampleSize;r++){
    if(r > 0){
      rm1 = r-1;
      Scale0 = (sumnoix2-2.0*mu0[rm1]*sumnoise+nnoise*mu0[rm1]*mu0[rm1])
	+(sumsigx2-2.0*mu1[rm1]*sumsignal+nsignal*mu1[rm1]*mu1[rm1]);
      Scale0 = 2.0/Scale0;
      lambda[r] = rgamma(halfsize,Scale0);
    }
    mu0[r] = rnorm(sumnoise/nnoise,sqrt(1.0/lambda[r]/nnoise));
    mu1[r] = rnorm(sumsignal/nsignal,sqrt(1.0/lambda[r]/nsignal));

    /* for the first probe */
    i = 0;
    state = X[i];
    z0 = R_pow_di((dt[i]-mu0[r]),2);
    z1 = R_pow_di((dt[i]-mu1[r]),2);
    p = 1.0/(1.0 + exp((z1-z0)*lambda[r]/2.0-2.0*X[i+1]*beta));
    if(runif(0,1)<p){
      X[i] = 1;
    }else {
      X[i] = -1;
    }
  
    if(state != X[i]){
      /*if state change from -1 to 1 */
      if(state == -1){
	sumnoise = sumnoise - dt[i];
	sumnoix2 = sumnoix2 - score2[i];
	nnoise = nnoise - 1;
	sumsignal = sumsignal + dt[i];
	sumsigx2 = sumsigx2 + score2[i];
	nsignal = nsignal + 1;
      }else {
	sumnoise = sumnoise + dt[i];
	sumnoix2 = sumnoix2 + score2[i];
	nnoise = nnoise + 1;
	sumsignal = sumsignal - dt[i];
	sumsigx2 = sumsigx2 - score2[i];
	nsignal = nsignal - 1;
      }
    }
    if(X[i] == 1 && r >= (*burning)){
      postX[i] = postX[i] + X[i];
    }

    for(i=1; i<nrowm1; i++){
      state = X[i];
      z0 = R_pow_di((dt[i]-mu0[r]),2);
      z1 = R_pow_di((dt[i]-mu1[r]),2);
      p = 1.0/(1.0 + exp((z1-z0)*lambda[r]/2.0-2.0*(X[i-1]+X[i+1])*beta));
      if(runif(0,1)<p){
	X[i] = 1;
      }else {
	X[i] = -1;
      }
      
      if(state != X[i]){
	/*if state change from 0 to 1 */
	if(state == -1){
	  sumnoise = sumnoise - dt[i];
	  sumnoix2 = sumnoix2 - score2[i];
	  nnoise = nnoise - 1;
	  sumsignal = sumsignal + dt[i];
	  sumsigx2 = sumsigx2 + score2[i];
	  nsignal = nsignal + 1;
	}else {
	  sumnoise = sumnoise + dt[i];
	  sumnoix2 = sumnoix2 + score2[i];
	  nnoise = nnoise + 1;
	  sumsignal = sumsignal - dt[i];
	  sumsigx2 = sumsigx2 - score2[i];
	  nsignal = nsignal - 1;
	}
      }
      if(X[i] == 1 && r >= (*burning)){
	postX[i] = postX[i] + 1;
      }
    }
    /* For the last probe */
    i = nrowm1;
    state = X[i];
    z0 = R_pow_di((dt[i]-mu0[r]),2);
    z1 = R_pow_di((dt[i]-mu1[r]),2);
    p = 1.0/(1.0 + exp((z1-z0)*lambda[r]/2.0-2.0*X[i-1]*beta)); 
    if(runif(0,1)<p){
      X[i] = 1;
    }else {
      X[i] = -1;
    }
    
    if(state != X[i]){
      /*if state change from -1 to 1 */
      if(state == -1){
	sumnoise = sumnoise - dt[i];
	sumnoix2 = sumnoix2 - score2[i];
	nnoise = nnoise - 1;
	sumsignal = sumsignal + dt[i];
	sumsigx2 = sumsigx2 + score2[i];
	nsignal = nsignal + 1;
      }else {
	sumnoise = sumnoise + dt[i];
	sumnoix2 = sumnoix2 + score2[i];
	nnoise = nnoise + 1;
	sumsignal = sumsignal - dt[i];
	sumsigx2 = sumsigx2 - score2[i];
	nsignal = nsignal - 1;
      }
    }
    if(X[i]==1 && r >= (*burning)){
      postX[i] = postX[i] + X[i];
    }   
    /* Metropolis random walk update beta */
    newBeta = beta + rnorm(0, *ransd);
    if(newBeta > (*minBeta) && newBeta < (*maxBeta)){
      sumxx = 0;
      for(i=0; i<nrowm1; i++){
	sumxx = sumxx + X[i]*X[i+1];
      }
      ratio = nrowm1*(log(cosh(beta))-log(cosh(newBeta)))+(newBeta-beta)*sumxx;
      if(ratio >= 0){
	beta = newBeta;
      }else {
	if(runif(0,1) < exp(ratio)){
	  beta = newBeta;	
	}
      }
    }
    postBeta[r] = beta;

    if(r%2000 == 0){
      R_CheckUserInterrupt();
      if(*verbose == 1){
	Rprintf("%d  ",r);
      }
    }
  }
  PutRNGstate();
  /*
  Rprintf("mNoise = %f, varNoise = %f \n",sumnoise/nnoise,(sumnoix2-sumnoise/nnoise*sumnoise)/(nnoise-1));
  Rprintf("mSignal = %f, varSignal = %f \n",sumsignal/nsignal,(sumsigx2-sumsignal/nsignal*sumsignal)/(nsignal-1));
  */

  if(*verbose == 1){
    Rprintf("\n");
  }
  if((nnoise <= 1) || (nsignal <= 1)){
    Rprintf("Warning: all probes are in the same state at the last MCMC iteration. \n NO enriched region is found!\n");
  }

  for(i=0; i<(*nrow); i++){
    postX[i] = (double)postX[i]/(*Size);
  }
  /*  free(score2); */
}



/* X is the region needed to be merged. Merged region will be put in Y*/
void MergeRegion(int *X,int *xrow, int *xcol, int *ycol,int *maxgap,int *Y,int *nregion){
  int **x, **y,i,j;
  x = (int **)R_alloc(*xrow,sizeof(int *));
  y = (int **)R_alloc(*xrow,sizeof(int *));

  if(x == NULL || y == NULL){
    error("Error: Fail to allocate memory! \n");
  }
  for(i=0; i < *xrow; i++){
    x[i] = (int *)R_alloc(*xcol, sizeof(int));
    y[i] = (int *)R_alloc(*ycol, sizeof(int));
    if(x[i] == NULL || y[i] == NULL){
      error("Error: Fail to allocate memory! \n");
    }
  }

  for(j=0; j < *xcol; j++){
    for(i=0; i < *xrow; i++){
      x[i][j] = X[i+(*xrow)*j];
      y[i][j] = 0;
    }
  }
  
  /* Just copy the first row */
  for(i=0; i < *ycol; i++){
    y[0][i] = Y[i];
  }

  j = 0;
  for(i=1; i < *xrow; i++){
    if(x[i][0] == y[j][0]){ /* if the same chromosome */
      if((x[i][1]-y[j][2]) < (*maxgap)){
        y[j][2] = x[i][1];
        y[j][4] = x[i][2];
      }else {
	j = j+1;
	y[j][0] = x[i][0];
	y[j][1] = x[i][1];
	y[j][2] = x[i][1];
	y[j][3] = x[i][2];
	y[j][4] = x[i][2];
      }
    }else{
      j = j+1;    
      y[j][0] = x[i][0];
      y[j][1] = x[i][1];
      y[j][2] = x[i][1];
      y[j][3] = x[i][2];
      y[j][4] = x[i][2];
    }    
  }

  *nregion = j; 
  for(j=0; j < *ycol; j++){
    for(i=0; i < *xrow; i++){
      Y[i+(*xrow)*j] = y[i][j];
    }
  }
  *nregion = *nregion + 1; /*because R index is from 1:n */
  /*
  for(i=0; i<(*xrow);i++){
    free(x[i]);
    free(y[i]);
  }
  */

  /*
  free(x);
  free(y);
  */
}

void fdr(int *klen, double *kp, int *blen, double *beta, double *efdr){
  int k,b;
  int *J;
  J = (int *)R_alloc(*klen,sizeof(int));
  if(J==NULL){
    error("Error: Fail to allocate memory!\n");
    /*    exit(1); */
  }

  for(k=0; k<(*klen); k++){
    J[k] = 0;
  }

  for(k=0; k<(*klen); k++){
    for(b=0; b<(*blen);b++){
      if(beta[b]<=kp[k]){
	J[k] = J[k] + 1;
	efdr[k] = efdr[k] + beta[b];
      }      
    }
    efdr[k] = efdr[k]/J[k];
  }
  /*  free(J); */
}

R_NativePrimitiveArgType iChip2Args[15] = {INTSXP, INTSXP, INTSXP, INTSXP, REALSXP, INTSXP, REALSXP, REALSXP, REALSXP, 
					   INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, INTSXP};
R_NativePrimitiveArgType iChip1Args[16] = {INTSXP,INTSXP,INTSXP,REALSXP,REALSXP,REALSXP,REALSXP,REALSXP,REALSXP,REALSXP,
					   INTSXP,REALSXP,REALSXP,REALSXP,REALSXP,INTSXP};
R_NativePrimitiveArgType MergeRegionArgs[7] =  {INTSXP,INTSXP,INTSXP,INTSXP,INTSXP,INTSXP,INTSXP};
R_NativePrimitiveArgType fdrArgs[5] =  {INTSXP,REALSXP,INTSXP,REALSXP,REALSXP};

static const R_CMethodDef CEntries[] = {
  {"iChip2", (DL_FUNC)&iChip2, 15, iChip2Args},
  {"iChip2", (DL_FUNC)&iChip2, 16, iChip1Args},
  {"MergeRegion", (DL_FUNC)&MergeRegion, 7, MergeRegionArgs},
  {"fdr", (DL_FUNC)&fdr, 5, fdrArgs},
  {NULL, NULL, 0}
}; 

/*
static const R_CMethodDef CEntries[] = {
  {"iChip2", (DL_FUNC)&iChip2, 15},
  {"MergeRegion", (DL_FUNC)&MergeRegion, 7},
  {"fdr", (DL_FUNC)&fdr, 5},
  {NULL, NULL, 0}
}; 
*/

void R_init_iChip(DllInfo *dll){
  R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
}
