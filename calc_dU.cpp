/*==========================================================
 * Roe_Flux.cpp - Roe's approximate Riemann solver for 2D
 *
 * The calling syntax is:
 *
 *		deltaU = calc_dU_2(dUdX,DX,n_cells)
 *
 * This is a MEX-file for MATLAB.
 * Copyright 2007-2008 The MathWorks, Inc.
 *
 *========================================================*/
/* $Revision: 1.1.10.2 $ */

/* See /usr/local/MATLAB/extern/examples for additional examples */

#include <math.h>
#include "mex.h"

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
  double *dUdX; /* n_fields x n_faces x n_cells_total input matrix (primitive variables gradient in each element) */
  double *dXE;  /* n_dims x n_faces x n_cells input matrix (vector from cell center to each cell edge)  */
  double *dU;   /* n_fields x n_faces x n_cells output matrix (del(U) dot r_ij) */
  int n_cells;  /* input : number of [non-ghost] cells in the mesh */

  //int *size1, *size2;                   /* size of input matrix */
  int n_faces, n_cells_tot, n_fields, n_dims;

  /* check for proper number of arguments */
  if(nrhs!=3) {
    mexErrMsgTxt("calc_dU: 3 inputs required.");
  }
  if(nlhs!=1) {
    mexErrMsgTxt("calc_dU: One output required.");
  }

  /* create a pointer to the data in the input matrices  */
  dUdX = mxGetPr(prhs[0]);
  dXE = mxGetPr(prhs[1]);
  n_cells = (int)mxGetScalar(prhs[2]);

  /* get dimensions of the input matrix */
  const int *size1 = mxGetDimensions(prhs[0]);
  const int *size2 = mxGetDimensions(prhs[1]);
  
  n_fields = size1[0];
  n_dims = size1[1];
  n_cells_tot = size1[2];
  
  n_faces = size2[1];
  
  /* check that first input argument is size 4x1 */
  if(size1[1]!=size2[0]) {
    mexErrMsgTxt("calc_dU.cpp: Input matrices must be of compatible sizes!");
  }

  if(n_cells>n_cells_tot) {
    mexErrMsgTxt("calc_dU.cpp: size of dUdX is less than the given n_cells!");
  }
  
  /* create the output matrix */
  const mwSize dims[3] = {n_fields, n_faces, n_cells};
  plhs[0] = mxCreateNumericArray(3,dims,mxDOUBLE_CLASS,mxREAL);

  /* get a pointer to the real data in the output matrix */
  dU = mxGetPr(plhs[0]);

  /* Calculate delta-U = dU_dX * dX 
   * Note: A[i,j,k] = A[i+ M*j + M*N*k] */
  int ijk1, ijk2, ijk3;
  for (int ic=0; ic<n_cells; ic++) {
    for (int face=0; face<n_faces; face++) {
      for (int field=0; field<n_fields; field++) {
        for (int dim=0; dim<n_dims; dim++) {
          ijk1 = field+n_fields*face+n_fields*n_faces*ic; // dU
          ijk2 = field+n_fields*dim+n_fields*n_dims*ic;   // dUdX
          ijk3 = dim+n_dims*face+n_dims*n_faces*ic;       // dXE
          dU[ijk1] += dUdX[ijk2]*dXE[ijk3];
        }
      }
    }
  }
}
