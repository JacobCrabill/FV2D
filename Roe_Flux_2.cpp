/*==========================================================
 * Roe_Flux.cpp - Roe's approximate Riemann solver for 2D
 *
 * The calling syntax is:
 *
 *		Fn = Roe_Flux(WL,WR,norm)
 *
 * This is a MEX-file for MATLAB.
 * Copyright 2007-2008 The MathWorks, Inc.
 *
 *========================================================*/
/* $Revision: 1.1.10.2 $ */

/* See /usr/local/MATLAB/extern/examples for additional examples */

#include <math.h>
#include "mex.h"

/* The computational routine */
void Roe_Flux(double *u_l, double *u_r, double *norm, double* fn)
{
  double p_l,p_r;
  double h_l, h_r;
  double sq_rho,rrho,hm,usq,am,am_sq,unm;
  double lambda0,lambdaP,lambdaM;
  double rhoun_l, rhoun_r,eps;
  double a1,a2,a3,a4,a5,a6,aL1,bL1;
  double v_l[2], v_r[2], um[2], du[4];

  double gamma = 1.4;

  /* velocities */
  for (int i=0;i<2;i++)  {
    v_l[i] = u_l[i+1]/u_l[0];
    v_r[i] = u_r[i+1]/u_r[0];
  }

  p_l=(gamma-1.0)*(u_l[3]- 0.5*u_l[0]*(v_l[0]*v_l[0] + v_l[1]*v_l[1]));
  p_r=(gamma-1.0)*(u_r[3]- 0.5*u_r[0]*(v_r[0]*v_r[0] + v_r[1]*v_r[1]));

  h_l = (u_l[3]+p_l)/u_l[0];
  h_r = (u_r[3]+p_r)/u_r[0];

  sq_rho = sqrt(u_r[0]/u_l[0]);
  rrho = 1./(sq_rho+1.);

  for (int i=0;i<2;i++) {
    um[i] = rrho*(v_l[i]+sq_rho*v_r[i]);
  }

  hm = rrho*(h_l+sq_rho*h_r);

  usq=0.;
  for (int i=0;i<2;i++)
    usq += 0.5*um[i]*um[i];

  am_sq   = (gamma-1.)*(hm-usq);
  am  = sqrt(am_sq);
  unm = 0.;
  for (int i=0;i<2;i++) {
    unm += um[i]*norm[i];
  }

  // Compute Euler flux (first part)
  rhoun_l = 0.;
  rhoun_r = 0.;
  for (int i=0;i<2;i++) {
    rhoun_l += u_l[i+1]*norm[i];
    rhoun_r += u_r[i+1]*norm[i];
  }

  fn[0] = rhoun_l + rhoun_r;
  fn[1] = rhoun_l*v_l[0] + rhoun_r*v_r[0] + (p_l+p_r)*norm[0];
  fn[2] = rhoun_l*v_l[1] + rhoun_r*v_r[1] + (p_l+p_r)*norm[1];
  fn[3] = rhoun_l*h_l   +rhoun_r*h_r;

  for (int i=0;i<4;i++) {
    du[i] = u_r[i]-u_l[i];
  }

  lambda0 = fabs(unm);
  lambdaP = fabs(unm+am);
  lambdaM = fabs(unm-am);

  // Entropy fix
  eps = 0.5*(abs(rhoun_l/u_l[0]-rhoun_r/u_r[0])+ abs(sqrt(gamma*p_l/u_l[0])-sqrt(gamma*p_r/u_r[0])));
  if(lambda0 < 2.*eps)
    lambda0 = 0.25*lambda0*lambda0/eps + eps;
  if(lambdaP < 2.*eps)
    lambdaP = 0.25*lambdaP*lambdaP/eps + eps;
  if(lambdaM < 2.*eps)
    lambdaM = 0.25*lambdaM*lambdaM/eps + eps;

  a2 = 0.5*(lambdaP+lambdaM)-lambda0;
  a3 = 0.5*(lambdaP-lambdaM)/am;
  a1 = a2*(gamma-1.)/am_sq;
  a4 = a3*(gamma-1.);

  a5 = usq*du[0]-um[0]*du[1]-um[1]*du[2]+du[3];
  a6 = unm*du[0]-norm[0]*du[1]-norm[1]*du[2];

  aL1 = a1*a5 - a3*a6;
  bL1 = a4*a5 - a2*a6;

  // Compute Euler flux (second part)
  fn[0] = 0.5*(fn[0] - (lambda0*du[0]+aL1));
  fn[1] = 0.5*(fn[1] - (lambda0*du[1]+aL1*um[0]+bL1*norm[0]));
  fn[2] = 0.5*(fn[2] - (lambda0*du[2]+aL1*um[1]+bL1*norm[1]));
  fn[3] = 0.5*(fn[3] - (lambda0*du[3]+aL1*hm   +bL1*unm));
}

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
  double *WL;       /* n_edges x 4 : input : Left state vector */
  double *WR;       /* n_edges x 4 : input : Right state vector  */
  double *NORM;     /* n_edges x 2 : input : Interface unit normal vector */
  double *AREA;     /* n_edges x 1 : input : Interface area */ 
  double *Fn;       /* n_edges x 4 : output : Normal flux vector */
  int nedges;

  /* check for proper number of arguments */
  if(nrhs!=4) {
      mexErrMsgIdAndTxt("Roe_Flux:nrhs","Four inputs required.");
  }
  if(nlhs!=1) {
      mexErrMsgIdAndTxt("Roe_Flux:nlhs","One output required.");
  }

  /* check that first input argument is size 4x1 */
  if(mxGetN(prhs[0])!=4) {
      mexErrMsgIdAndTxt("Roe_Flux:WrongSize","Input WL must be a n_edgesx4 matrix.");
  }

  /* create a pointer to the data in the input matrices  */
  WL = mxGetPr(prhs[0]);
  WR = mxGetPr(prhs[1]);
  NORM = mxGetPr(prhs[2]);
  AREA = mxGetPr(prhs[3]);
  
  /* get dimensions of the input matrix */
  nedges = mxGetM(prhs[0]);

  /* create the output matrix */
  plhs[0] = mxCreateDoubleMatrix(nedges,4,mxREAL);

  /* get a pointer to the real data in the output matrix */
  Fn = mxGetPr(plhs[0]);

  /* Additional temporary storage variables */
  double wl[4], wr[4], fn[4], norm[2];
  
  /* call the computational routine */
  for (int i=0; i<nedges; i++) {
    for (int j=0; j<4; j++) {
      // Matlab column-major: [j*M + i]
      wl[j] = WL[j*nedges + i];
      wr[j] = WR[j*nedges + i];
      fn[j] = 0;
    }
    norm[0] = NORM[i];
    norm[1] = NORM[nedges+i];
    
    Roe_Flux(wl,wr,norm,fn);

    for (int j=0; j<4; j++) {
      Fn[j*nedges+i] = AREA[i]*fn[j];
    }
  }
}






















