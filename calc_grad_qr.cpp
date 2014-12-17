/*==========================================================
 * calc_grad_qr.cpp - Calculate gradient of a solution vector 
 * in the cells of an unstructured grid using a QR-based 
 * least-squares approach [dX*dUdX = dU, dX=QR]
 *
 * The calling syntax is:
 *
 *		dUdX = calc_grad_qr(U,c2ac,c2nac,Q,R)
 *
 * This is a MEX-file for MATLAB.
 * Written 5/30/2014 by Jacob Crabill
 *
 *========================================================*/

#include <math.h>
#include "mex.h"

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
  double *U;     /* input : n_cells_total x n_fields : solution vector in each cell (including ghost cells) */
  double *c2ac;  /* input : n_cells x max_nc : list of all cells encircling each interior cell */
  double *c2nac; /* input : n_cells x 1 : number of cells encircling each interior cell */
  double *Q;     /* input : max_nc x max_nc x n_cells : Q' of dX = QR */
  double *R;     /* input : n_dims x max_nc x n_cells : R of dX = QR */

  double *dUdX;  /* output : n_fields x n_dims x n_cells_total : grad of primitive variables in each cell */

  int n_fields, n_dims, n_cells, n_cells_tot;
  int max_nc;

  /* check for proper number of arguments */
  if(nrhs!=4) {
      mexErrMsgTxt("calc_grad_qr.cpp: 4 inputs required.");
  }
  if(nlhs!=1) {
      mexErrMsgTxt("calc_grad_qr.cpp: One output required.");
  }

  /* create a pointer to the data in the input matrices  */
  U = mxGetPr(prhs[0]);
  c2ac = mxGetPr(prhs[1]);
  c2nac = mxGetPr(prhs[2]);
  dXinv = mxGetPr(prhs[3]);

  /* get dimensions of the input matrices */
  n_cells_tot = mxGetM(prhs[0]);
  n_cells = mxGetM(prhs[1]);
  n_fields = mxGetN(prhs[0]); 
  n_dims = mxGetM(prhs[3]);
  max_nc = mxGetN(prhs[1]);

  /* check a few argument sizes against expectations */
  if(n_dims!=2) {
    mexPrintf("Expected n_dims to be 2; got %i instead.",n_dims);
    mexErrMsgTxt("calc_grad.cpp: Wrong number of physical dimensions found.");
  }
  
  /* create the output matrix */
  //if (mxGetNumberOfElements(plhs[0]) != n_fields*n_dims*n_cells_tot) {
    const mwSize dims[3] = {n_fields, n_dims, n_cells_tot};
    plhs[0] = mxCreateNumericArray(3,dims,mxDOUBLE_CLASS,mxREAL);
  //}

  /* get a pointer to the real data in the output matrix */
  dUdX = mxGetPr(plhs[0]);
  
  int nc, ic2, iu1, iu2, igradu, ix, idu;
  double *dU;
  for (int ic=0; ic<n_cells; ic++) {
    nc = c2nac[ic];
    //mexPrintf("\nic = %i, nc = %i, max_nc = %i\n",ic,nc,max_nc);
    // Get the dU matrix
    //mexPrintf("ic=%i: Getting dU\n",ic);
    dU = new double[n_fields*nc*2];//*sizeof(double)];
    for (int j=0; j<nc; j++) {
      ic2 = c2ac[ic+j*n_cells]-1;
      //mexPrintf("%i: ic2=%i\n",ic,ic2);
      for (int f=0; f<n_fields; f++) {
        iu1 = ic  + f*n_cells_tot;
        iu2 = ic2 + f*n_cells_tot;
        dU[j+f*nc] = U[iu2]-U[iu1];
        //mexPrintf("%i: dU[%i,%i] = %f",ic,j,f,dU[j+f*nc]);
      }
      //mexPrintf("\n");
    }
    
    // Get dU2 = Q'*dU
    dU2 = new double[n_fields*nc*2];
    for (i=0; i<n_fields*nc*2; i++) 
      dU2[i] = 0;

    for (int j=0; j<nc; j++) {
      for (int k=0; k<nc; k++) {
        for (int f=0; f<n_fields; f++) {
          dU2[j+f*nc] += Q[j+k*max_nc]*dU[k+f*nc];
        }
      }
    }

    // Get dUdX = R\dU2
    //mexPrintf("ic=%i: Getting dUdX\n",ic);
    for (int f=0; f<n_fields; f++) {
      for (int k=0; k<n_dims; k++) {
        for (int j=0; j<nc; j++) {
          igradu = f + k*n_fields + ic*n_dims*n_fields;
          ix = k + j*max_nc + ic*n_dims*max_nc;
          idu = j+f*nc;
          dUdX[igradu] += dXinv[ix]*dU[idu];
        }
      }
    }
    //mexPrintf("ic=%i: deleting dU\n",ic);
    delete[] dU;
    delete[] dU2;
  }
}
