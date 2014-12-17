/*==========================================================
 * convective_flux.cpp - Calculate convective portion of
 * normal interface flux from the Euler equations
 *
 * The calling syntax is:
 *
 *		Qn = convective_flux(WL,WR,norm)
 *
 * This is a MEX-file for MATLAB.
 * Copyright 2007-2008 The MathWorks, Inc.
 *
 *========================================================*/

#include <math.h>
#include "mex.h"

/* The computational routine */
void Convective_Flux(double *u_l, double *u_r, double *norm, double area, double* fn)
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

  /* pressure */
  p_l=(gamma-1.0)*(u_l[3]- 0.5*u_l[0]*(v_l[0]*v_l[0] + v_l[1]*v_l[1]));
  p_r=(gamma-1.0)*(u_r[3]- 0.5*u_r[0]*(v_r[0]*v_r[0] + v_r[1]*v_r[1]));

  /* enthalpy */
  h_l = (u_l[3]+p_l)/u_l[0];
  h_r = (u_r[3]+p_r)/u_r[0];
  
  /* face-normal momentum */
  rhoun_l = 0.;
  rhoun_r = 0.;
  for (int i=0;i<2;i++) {
    rhoun_l += u_l[i+1]*norm[i];
    rhoun_r += u_r[i+1]*norm[i];
  }
  
  /* Euler flux */
  fn[0] = area*0.5*(rhoun_l + rhoun_r);
  fn[1] = area*0.5*(rhoun_l*v_l[0] + rhoun_r*v_r[0] + (p_l+p_r)*norm[0]);
  fn[2] = area*0.5*(rhoun_l*v_l[1] + rhoun_r*v_r[1] + (p_l+p_r)*norm[1]);
  fn[3] = area*0.5*(rhoun_l*h_l   +rhoun_r*h_r);
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
      mexErrMsgIdAndTxt("convective_flux:nrhs","Four inputs required.");
  }
  if(nlhs!=1) {
      mexErrMsgIdAndTxt("convective_flux:nlhs","One output required.");
  }

  /* check that first input argument is size 4x1 */
  if(mxGetN(prhs[0])!=4) {
      mexErrMsgIdAndTxt("convective_flux:WrongSize","Input WL must be a n_edgesx4 matrix.");
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
  double wl[4], wr[4], fn[4], norm[2], area;
  
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
    area = AREA[i];
    
    Convective_Flux(wl,wr,norm,area,fn);

    for (int j=0; j<4; j++) {
      Fn[j*nedges+i] = fn[j];
    }
  }
}






















