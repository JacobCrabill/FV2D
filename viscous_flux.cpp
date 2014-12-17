/*==========================================================
 * viscous_flux.cpp - Calculate viscous terms of normal
 * interface flux from the Navier-Stokes equations
 *
 * The calling syntax is:
 *
 *		Fn = viscous_flux(tauL,tauR,phiL,phiR,norm,dA,cp,Pr,mu_inf)
 *
 * This is a MEX-file for MATLAB.
 * Copyright 2007-2008 The MathWorks, Inc.
 *
 *========================================================*/

#include <math.h>
#include "mex.h"

/* The computational routine */
void Viscous_Flux(double *grad_u, double *phi, const double &Cp, const double &Pr, 
    const double &R, double *F, double* G)
{
  double du_dx, du_dy, dv_dx, dv_dy;
  double dT_dx, dT_dy;
  double mu, T;
  double tauxx, tauxy, tauyy, diag;
  double Qx, Qy;
  double u, v;
  
  // from HiFiLES: double T_ref = 291.15; double mu_ref = 1.827E-5; 
  // From aerojet.ucdavis.edu/fluenthelp/html/ug/node337.htm:
  double T_s = 110.56;      // K - Sutherland's ref. temperature
  double T_ref = 273.11;    // K
  double mu_ref = 1.716E-5; // kg/m-s
  
  u = phi[1];
  v = phi[2];
  T = phi[3]/(R*phi[0]);

  if (isnan(T)) {
    mexPrintf("Error! NaN T in viscous_flux.cpp!\n");
  }

  //mexPrintf("T = %f\n",T);

  mu = mu_ref*pow(T/T_ref,1.5)*((T_ref+T_s)/(T+T_s));
  
  du_dx = grad_u[0];
  du_dy = grad_u[1];
  dv_dx = grad_u[2];
  dv_dy = grad_u[3];
  dT_dx = grad_u[4];
  dT_dy = grad_u[5];

  diag = (du_dx + dv_dy)/3.0;

  tauxx = 2.0*mu*(du_dx-diag);
  tauxy = mu*(du_dy + dv_dx);
  tauyy = 2.0*mu*(dv_dy-diag);
  Qx = mu*Cp/Pr*dT_dx;
  Qy = mu*Cp/Pr*dT_dy;

  F[0] = 0;
  F[1] = -tauxx;
  F[2] = -tauxy;
  F[3] = -u*tauxx - v*tauxy - Qx;

  G[0] = 0.0;
  G[1] = -tauxy;
  G[2] = -tauyy;
  G[3] = -u*tauxy - v*tauyy - Qy;
  
  /*mexPrintf("mu = %f\n",mu);
  mexPrintf("tauxx = %f\n",tauxx);
  mexPrintf("tauxy = %f\n",tauxy);
  mexPrintf("tauyy = %f\n",tauyy);*/
}

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
  double *tauL;          /* n_edges x 6 : input : Left state vector of gradients */
  double *tauR;          /* n_edges x 6 : input : Right state vector of gradients */
  double *phiL;          /* n_edges x 4 : input : Left state vector of primitives */
  double *phiR;          /* n_edges x 4 : input : Right state vector of primitives */
  double *NORM;          /* n_edges x 2 : input : Interface unit normal vector */
  double *AREA;          /* n_edges x 1 : input : Interface area */ 
  double Cp, Pr, R;      /* scalars : input : Physical constants */
  double *Fn;            /* n_edges x 4 : output : Normal flux vector */
  int nedges;

  /* check for proper number of arguments */
  if(nrhs!=9) {
      mexErrMsgTxt("viscous_flux.cpp : 9 inputs required.");
  }
  if(nlhs!=1) {
      mexErrMsgTxt("viscouse_flux.cpp : One output required.");
  }

  /* check that first input argument is size 4x1 */
  if(mxGetN(prhs[0])!=6) {
      mexErrMsgTxt("viscous_flux.cpp : First input (tauL) must be a n_edges x 6 matrix.");
  }

  /* create a pointer to the data in the input matrices  */
  tauL = mxGetPr(prhs[0]);
  tauR = mxGetPr(prhs[1]);
  phiL = mxGetPr(prhs[2]);
  phiR = mxGetPr(prhs[3]);
  NORM = mxGetPr(prhs[4]);
  AREA = mxGetPr(prhs[5]);
  Cp = mxGetScalar(prhs[6]);
  Pr = mxGetScalar(prhs[7]);
  R = mxGetScalar(prhs[8]);

  /* get dimensions of the input matrix */
  nedges = mxGetM(prhs[0]);

  /* create the output matrix */
  plhs[0] = mxCreateDoubleMatrix(nedges,4,mxREAL);

  /* get a pointer to the real data in the output matrix */
  Fn = mxGetPr(plhs[0]);

  /* Additional temporary storage variables */
  double taul[6], taur[6], phil[4], phir[4], norm[2], area;
  double FL[4], FR[4], GL[4], GR[4];
  //double PHI[4], TAU[6];
  
  //mexPrintf("\n");
  /* call the computational routine */
  // Matlab column-major: [i + M*j + N*k]
  for (int ie=0; ie<nedges; ie++) {
    for (int j=0; j<4; j++) {
      phil[j] = phiL[j*nedges + ie];
      phir[j] = phiR[j*nedges + ie];
      //PHI[j] = 0.5*(phil[j]+phir[j]);
    }
    for (int j=0; j<6; j++) {
      taul[j] = tauL[j*nedges + ie];
      taur[j] = tauR[j*nedges + ie];
      //TAU[j] = 0.5*(taul[j]+taur[j]);
    }
    
    Viscous_Flux(taul,phil,Cp,Pr,R,FL,GL);
    Viscous_Flux(taur,phir,Cp,Pr,R,FR,GR);
    //Viscous_Flux(TAU,PHI,Cp,Pr,R,FR,GR);
    
    norm[0] = NORM[ie+nedges*0];
    norm[1] = NORM[ie+nedges*1];
    area = AREA[ie];
    for (int j=0; j<4; j++) {
      Fn[j*nedges+ie] = 0.5*area*((FL[j]+FR[j])*norm[0]+(GL[j]+GR[j])*norm[1]);
      //Fn[ie+j*nedges] = area*(FR[j]*norm[0]+GR[j]*norm[1]);
    }
  }
}
