/*==========================================================
 * apply_bcs_grad_phi.cpp - Apply boundary conditions to
 * solution gradients at the cell-centers of an unstructured
 * mesh, using a pre-computed QR-type approace
 *
 * The calling syntax is:
 *
 *		dphidX = apply_bcs_grad_phi(dphidX,c2c,c2e,bc,unorm,n_cells)
 *
 * This is a MEX-file for MATLAB.
 * Written 7/24/2014 by Jacob Crabill
 *
 *========================================================*/

#include <math.h>
#include "mex.h"

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
  //double *dphidX_in; /* input : n_fields x n_dims x n_cells_total : grad of primitive vars in each cell (incl. ghost cells) */
  double *c2c;       /* input : n_cells_tot x max_nc : list of cells surrounding each cell */
  double *c2e;       /* input : n_cells_tot x max_nf : Edges for all cells in mesh */
  double *bc;        /* input : n_edges x 1 : Boundary condition for all edges in mesh */
  double *unorm;     /* input : n_edges x n_dims : Unit normal for all 'real' edges */

  double *dphidX;    /* output : n_fields x n_dims x n_cells_total : grad of primitive vars in each cell (incl. ghost cells) */

  int n_fields, n_dims, n_edges, n_cells, n_cells_tot, max_nc;

  // Misc. indices
  int icg, ici, bcg, ie;
  int i11,i12,i21,i22,i31,i32,i41,i42;

  // Dummy vector, Normal vector, Tangential vector
  double blah[2], nhat[2], that[2];
  double isqrt2 = 1/sqrt(2);

  double dxdn, dxdt, dydn, dydt;
  double dudn, dudt, dvdn, dvdt, dpdn, dpdt, drhodn, drhodt;
  double dundn, dundt, dutdn, dutdt, dundx, dundy, dutdx, dutdy, dundu, dundv, dutdu, dutdv;

  /* check for proper number of arguments */
  if(nrhs!=6) {
      mexErrMsgTxt("apply_bcs_grad_phi.cpp: 6 inputs required.");
  }
  if(nlhs!=1) {
      mexErrMsgTxt("apply_bcs_grad_phi.cpp: One output required.");
  }

  /* create a pointer to the data in the input matrices  */
  //dphidX_in = mxGetPr(prhs[0]);
  c2c = mxGetPr(prhs[1]);
  c2e = mxGetPr(prhs[2]);
  bc = mxGetPr(prhs[3]);
  unorm = mxGetPr(prhs[4]);
  n_cells = (int)mxGetScalar(prhs[5]);

  /* get dimensions of the input matrices */
  max_nc = mxGetN(prhs[1]);
  n_edges = mxGetM(prhs[4]);

  const mwSize *dims = mxGetDimensions(prhs[0]);
  n_fields = dims[0];
  n_dims = dims[1];
  n_cells_tot = dims[2];

  /* Setup Output Matrix */
  plhs[0] = mxDuplicateArray(prhs[0]);

  /* Get pointer to output matrix */
  dphidX = mxGetPr(plhs[0]);

  for (icg=n_cells; icg<n_cells_tot; icg++) {
    ici = c2c[icg];
    ie = c2e[icg]; // boundary edge is first edge for ghost cell
    bcg = bc[ie];

    // NOTE: ghost cells ALWAYS on 'right' of all boundary edges
    if (bcg<5) {
      // Inlet/Outlet - Copy
      for (int f=0; f<n_fields; f++) {
        for (int k=0; k<n_dims; k++) {
          dphidX[f+n_fields*(k+icg*n_dims)] = dphidX[f+n_fields*(k+n_dims*icg)];
        }
      }
    }else if (bcg==5) {
      //mexPrintf("In bcg==5\n");
      // Slip Wall Mirror velocities about wall
      // du/dx and dv/dy are the same
      //mexPrintf("dphidX_icg = %f\n",dphidX[1+n_fields*(0+n_dims*icg)]);
      //mexPrintf("dphidX_ici = %f\n",dphidX[1+n_fields*(0+n_dims*ici)]);
      dphidX[1+n_fields*(0+n_dims*icg)] = dphidX[1+n_fields*(0+n_dims*ici)];
      dphidX[2+n_fields*(1+n_dims*icg)] = dphidX[2+n_fields*(1+n_dims*ici)];
      // du/dy and dv/dx are mirrored
      dphidX[1+n_fields*(1+n_dims*icg)] = -dphidX[1+n_fields*(1+n_dims*ici)];
      dphidX[2+n_fields*(0+n_dims*icg)] = -dphidX[2+n_fields*(0+n_dims*ici)];
      // drho/dx and dp/dx are same
      dphidX[0+n_fields*(0+n_dims*icg)] = dphidX[0+n_fields*(0+n_dims*ici)];
      dphidX[3+n_fields*(0+n_dims*icg)] = dphidX[3+n_fields*(0+n_dims*ici)];
      // drho/dy and dp/dy are mirrored
      dphidX[0+n_fields*(1+n_dims*icg)] = -dphidX[0+n_fields*(1+n_dims*ici)];
      dphidX[3+n_fields*(1+n_dims*icg)] = -dphidX[3+n_fields*(1+n_dims*ici)];
      //mexPrintf("Leaving bcg==5");
    }else if (bcg==6 || bcg==7) {
      // Adiabatic Wall or Isothermal Wall
      // Do coordinate transformation at edges from x,y to  normal,
      // tangential coordinates
      nhat[0] = unorm[ie]; nhat[1] = unorm[ie+n_edges]; // or n_edges_tot?
      // Create dummy normalized vector
      blah[0] = isqrt2;    blah[1] = isqrt2;
      // Get tangential vector
      that[0] = blah[0] - (nhat[0]*blah[0]+nhat[1]*blah[1])*nhat[0];
      that[1] = blah[1] - (nhat[0]*blah[0]+nhat[1]*blah[1])*nhat[1];
      that[0] = that[0] / sqrt(that[0]*that[0]+that[1]*that[1]);
      that[1] = that[1] / sqrt(that[0]*that[0]+that[1]*that[1]);

      dundu = nhat[0]; dundv = nhat[1];
      dutdu = that[0]; dutdv = that[1];
      dxdn = 1./nhat[0]; dydn = 1./nhat[1];
      dxdt = 1./that[0]; dydt = 1./that[1];

      i11 = 0+n_fields*(0+n_dims*ici); i12 = 0+n_fields*(1+n_dims*ici);
      i21 = 1+n_fields*(0+n_dims*ici); i22 = 1+n_fields*(1+n_dims*ici);
      i31 = 2+n_fields*(0+n_dims*ici); i32 = 2+n_fields*(1+n_dims*ici);
      i41 = 3+n_fields*(0+n_dims*ici); i42 = 3+n_fields*(1+n_dims*ici);
      // dudn = dudx*dxdn + dudy*dydn, etc.
      dudn = dphidX[i21]*dxdn + dphidX[i22]*dydn;
      dvdn = dphidX[i31]*dxdn + dphidX[i32]*dydn;
      dudt = dphidX[i21]*dxdt + dphidX[i22]*dydt;
      dvdt = dphidX[i31]*dxdt + dphidX[i32]*dydt;
      drhodn = dphidX[i11]*dxdn + dphidX[i12]*dydn;
      dpdn = dphidX[i41]*dxdn + dphidX[i42]*dydn;
      drhodt = dphidX[i11]*dxdt + dphidX[i12]*dydt;
      dpdt = dphidX[i41]*dxdt + dphidX[i42]*dydt;

      // dun/dn = dun/du * du/dn + dun/dv * dv/dn, etc.
      // dut/dn, dun/dn are same
      // dut/dt, dun/dt are mirrored
      dundn = dundu*dudn + dundv*dvdn;
      dundt =-dundu*dudt - dundv*dvdt;
      dutdn = dutdu*dudn + dutdv*dvdn;
      dutdt =-dutdu*dudt - dutdv*dvdt;

      // Transform un, ut to x,y
      dundx = dundn/dxdn + dundt/dxdt;
      dundy = dundn/dydn + dundt/dydt;
      dutdx = dutdn/dxdn + dutdt/dxdt;
      dutdy = dutdn/dydn + dutdt/dydt;

      i11 = 0+n_fields*(0+n_dims*icg); i12 = 0+n_fields*(1+n_dims*icg);
      i21 = 1+n_fields*(0+n_dims*icg); i22 = 1+n_fields*(1+n_dims*icg);
      i31 = 2+n_fields*(0+n_dims*icg); i32 = 2+n_fields*(1+n_dims*icg);
      i41 = 3+n_fields*(0+n_dims*icg); i42 = 3+n_fields*(1+n_dims*icg);
      // Transform all derivatives back to x,y
      dphidX[i21] = (1/dundu)*dundx + (1/dutdu)*dutdx; // dudx
      dphidX[i22] = (1/dundu)*dundy + (1/dutdu)*dutdy; // dudy
      dphidX[i31] = (1/dundv)*dundx + (1/dutdv)*dutdx; // dvdx
      dphidX[i32] = (1/dundv)*dundy + (1/dutdv)*dutdy; // dvdy
      // drho/dt, dp/dt are same
      // drho/dn, dp/dn are mirrored
      dphidX[i11] = -drhodn/dxdn + drhodt/dxdt; // drhodx
      dphidX[i12] = -drhodn/dydn + drhodt/dydt; // drhody
      dphidX[i41] = -dpdn/dxdn + dpdt/dxdt;     // dpdx
      dphidX[i42] = -dpdn/dydn + dpdt/dydt;     // dpdy
    }
  }
}
