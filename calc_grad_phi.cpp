/*==========================================================
 * calc_grad.cpp - Calculate gradient of a solution vector 
 * in the cells of an unstructured grid using a vertex- and
 * edge-based approach
 *
 * The calling syntax is:
 *
 *		dUdX = calc_grad(U,v2c,v2nc,xv,e2v,v2ne,c2v,c2nv,gv2iv,gv2bc,gv2be,unorm)
 *
 * This is a MEX-file for MATLAB.
 * Written 5/22/2014 by Jacob Crabill
 *
 *========================================================*/

#include <math.h>
#include "mex.h"

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
  double *U;     /* input : n_cells_total x n_fields : solution vector in each cell (including ghost cells) */
  double *v2c;   /* input : n_verts_tot x max_nc : list of cells surrounding each vertex */
  double *v2nc;  /* input : n_verts_tot x 1 : # of cells surrounding each vertex */
  double *xv;    /* input : n_verts_tot x n_dims : [x,y] position of each vertex */
  double *e2v;   /* input : n_edges_tot x 2 : list of edges in mesh */
  double *v2ne;  /* input : n_verts_tot x 1 : # of edges using each vertex */
  double *c2v;   /* input : n_cells_tot x max_nv : Vertices for all cells in mesh */
  double *c2nv;  /* input : n_cells_tot x 1 : # of vertices for each cell */
  double *gv2iv; /* input : n_ghost_verts x 1 : Ghost nodes to 'mirrored' internal node */
  double *gv2bc; /* input : n_ghost_verts x 1 : Boundary condition of  ghost node */
  double *gv2be; /* input : n_ghost_verts x 1 : Boundary edge of ghost node */
  double *unorm; /* input : n_edges x n_dims : Unit normal for all 'real' edges */
    
  double *dUdX;  /* output : n_fields x n_dims x n_cells_total : grad of primitive variables in each cell */

  int n_fields, n_dims, n_edges, n_verts, n_cells_tot, n_verts_tot, n_edges_tot;
  int n_ghost_verts;

  /* check for proper number of arguments */
  if(nrhs!=12) {
      mexErrMsgTxt("calc_grad.cpp: 12 inputs required.");
  }
  if(nlhs!=1) {
      mexErrMsgTxt("calc_grad.cpp: One output required.");
  }

  /* create a pointer to the data in the input matrices  */
  U = mxGetPr(prhs[0]);
  v2c = mxGetPr(prhs[1]);
  v2nc = mxGetPr(prhs[2]);
  xv = mxGetPr(prhs[3]);
  e2v = mxGetPr(prhs[4]);
  v2ne = mxGetPr(prhs[5]);
  c2v = mxGetPr(prhs[6]);
  c2nv = mxGetPr(prhs[7]);
  gv2iv = mxGetPr(prhs[8]);
  gv2bc = mxGetPr(prhs[9]);
  gv2be = mxGetPr(prhs[10]);
  unorm = mxGetPr(prhs[11]);

  /* get dimensions of the input matrices */
  n_cells_tot = mxGetM(prhs[0]);
  n_verts_tot = mxGetM(prhs[1]);
  n_fields = mxGetN(prhs[0]); 
  n_dims = mxGetN(prhs[3]);
  n_edges_tot = mxGetM(prhs[4]);
  n_ghost_verts = mxGetM(prhs[9]);
  n_edges = mxGetM(prhs[11]);
  
  n_verts = n_verts_tot - n_ghost_verts;

  /* check a few argument sizes against expectations */
  if(n_dims!=2) {
    mexPrintf("Expected n_dims to be 2; got %i instead.",n_dims);
    mexErrMsgTxt("calc_grad.cpp: Wrong number of physical dimensions found.");
  }
  
  /* create the output matrix */
  const mwSize dims[3] = {n_fields, n_dims, n_cells_tot};
  plhs[0] = mxCreateNumericArray(3,dims,mxDOUBLE_CLASS,mxREAL);

  /* get a pointer to the real data in the output matrix */
  dUdX = mxGetPr(plhs[0]);

  // Solution vector at vertices.
  double *Uv = new double[n_verts_tot*n_fields*sizeof(double)];
  // Solution gradient at vertices.  Indexing: (vertex,field,dim)
  double *gradUv = new double[n_verts_tot*n_fields*n_dims*sizeof(double)];

  // Initialize temporary solution & gradient vectors to 0
  for (int iv=0; iv<n_verts_tot; iv++) {
    for (int f=0; f<n_fields; f++) {
      Uv[f+iv*n_fields] = 0;
      for (int dim=0; dim<n_dims; dim++) {
        gradUv[f+dim*n_fields+iv*n_dims*n_fields] = 0;
      }
    }
  }

  /* (1) Get value of solution at (interior) vertices by averaging surrounding 
   *     cell-center values */
  int icn;
  for (int iv=0; iv<n_verts; iv++) {
    for (int j=0; j<v2nc[iv]; j++) {
      icn = v2c[iv+j*n_verts_tot]-1;
      for (int f=0; f<n_fields; f++) {
        Uv[f+iv*n_fields] += U[icn+f*n_cells_tot];
      }
    }
  }

  /* (2) Use gv2iv to get 'good' values at the ghost vertices */
  int ivi, ivg, icg, bc, ie;
  double Un;
  for (int iv=0; iv<n_ghost_verts; iv++) {
    ivg = iv + n_verts;
    ivi = gv2iv[iv]-1;
    bc = gv2bc[iv];
  
    if (bc<=4) // inlet/outlet - use avg. of ghost cells
    {
      for (int j=0; j<v2nc[ivg]; j++) {
        icg = v2c[ivg+j*n_verts_tot]-1;
        for (int f=0; f<n_fields; f++) {
          Uv[f+ivg*n_fields] += U[icg+f*n_cells_tot]/v2nc[ivg];
        }
      }
    }
    else if (bc==5) // slip wall - reflect
    {
      ie = gv2be[iv]-1;
      Un = 0;
      for (int k=0; k<n_dims; k++) {
        Un += Uv[k+1+ivi*n_fields]*fabs(unorm[ie+n_edges*k]);
      }
      for (int k=0; k<n_dims; k++) {
        Uv[k+1+ivg*n_fields] -= 2*Un*fabs(unorm[ie+n_edges*k]);
      }
      Uv[0+ivg*n_fields] = Uv[0+ivi*n_fields];
      Uv[n_fields+ivg*n_fields] = Uv[n_fields+ivi*n_fields];
    }
    else if (bc==6 || bc==7) // no-slip wall - mirror
    {
      for (int k=0; k<n_dims; k++) {
        Uv[k+1+ivg*n_fields] = -Uv[k+1+ivi*n_fields];
      }
      Uv[0+ivg*n_fields] = Uv[0+ivi*n_fields];
      Uv[n_fields+ivg*n_fields] = Uv[n_fields+ivi*n_fields];
    }
  }
  
  /* (3) Get the gradient at each vertex by averaging gradients along each edge */
  int idu1, idu2, iv1, iv2, iu1, iu2, ix1, ix2;
  double gradU, dn, dUdn, dX[2];

  for (int ie=0; ie<n_edges; ie++) {
    iv1 = e2v[ie+n_edges_tot*0]-1;
    iv2 = e2v[ie+n_edges_tot*1]-1;

    // Find edge vector (direction for gradient)
    dn = 0;
    for (int dim=0; dim<n_dims; dim++) {
      ix1 = iv1 + dim*n_verts_tot;
      ix2 = iv2 + dim*n_verts_tot;
      dX[dim] = xv[ix2]-xv[ix1];
      dn += dX[dim]*dX[dim];
    }
    dn = sqrt(dn);

    // Normalize
    for (int dim=0; dim<n_dims; dim++) {
      dX[dim] /= dn;
    }

    // Get gradient along direction of edge
    for (int f=0; f<n_fields; f++) {
      iu1 = f + iv1*n_fields;
      iu2 = f + iv2*n_fields;
      dUdn = (Uv[iu2]-Uv[iu1])/dn;

      for (int dim=0; dim<n_dims; dim++) {
        idu1 = f + dim*n_fields + iv1*n_dims*n_fields;
        idu2 = f + dim*n_fields + iv2*n_dims*n_fields;

        gradU = dUdn * dX[dim];
        gradUv[idu1] += gradU/v2ne[iv1];
        gradUv[idu2] += gradU/v2ne[iv2];
      }
    }
  }

  /* (4) Get the gradient in each cell by averaging the gradients at their vertices */
  int i1, i2, iv;
  for (int ic=0; ic<n_cells_tot; ic++) {
    for (int j=0; j<c2nv[ic]; j++) {
      iv = c2v[ic+j*n_cells_tot]-1;
      for (int f=0; f<n_fields; f++) {
        for (int dim=0; dim<n_dims; dim++) {
          i1 = f + dim*n_fields + ic*n_dims*n_fields;
          i2 = f + dim*n_fields + iv*n_dims*n_fields;
          dUdX[i1] += gradUv[i2]/c2nv[ic];
        }
      }
    }
  }

  delete[] gradUv;
  delete[] Uv;
}
