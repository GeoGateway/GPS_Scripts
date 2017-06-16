/*******************************************************************************
MODULE NAME
label_using_gauss_hmm

AUTHOR
Robert Granat

DESCRIPTION
Mex interface to rdahmm

COMPILE
make

NOTES
Not yet available

RG
*******************************************************************************/
/* Mex library */
#include "mex.h"

/* C library */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

/* CP library */
#include "cp.h"

/* UT library */
#include "ut.h"

/* NR library */
#include "nr.h"

/* DA library */
#include "da.h"

/* local header files */
#include "rdahmm_arguments.h"
#include "rdahmm_setup.h"
#include "rdahmm_control.h"
#include "rdahmm_gauss_output.h"
#include "rdahmm_generic.h"
#include "rdahmm.h"

/* this program's header */

/*******************************************************************************
 MEXFUNCTION
*******************************************************************************/
/* The gateway function */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  /* inputs */
  double            **cont_data;
  Gaussian_HMM      hmm;
  RDAEM_Parameters  params;
  int               viterbi;

  /* outputs */
  double       L;
  int          *Q;

  /* problem dimensions */
  int     T;
  int     D;
  int     N;

  /* other variables */
  int    i;
  double *p;

  /* check for proper number of arguments */
  if (nrhs < 2) 
  {
    mexErrMsgIdAndTxt("MyToolbox:label_using_gauss_hmm:nrhs", 
                      "At least two inputs required.");
  }
  if (nlhs < 1) 
  {
    mexErrMsgIdAndTxt("MyToolbox:label_using_gauss_hmm:nlhs", 
                      "At least one output required.");
  }
  if (nlhs > 2)
  {
    mexErrMsgIdAndTxt("MyToolbox:label_using_gauss_hmm:nlhs", 
                      "Too many outputs (max 2).");
  }

  /* get array dimensions */
  T = (int) mxGetN(prhs[0]);
  D = (int) mxGetM(prhs[0]);
  N = (int) mxGetScalar(mxGetField(prhs[1], 0, "N"));

  /* allocate working memory */
  cont_data = NR_dmatrix(1, T, 1, D);
  allocate_gauss_hmm(&hmm, N, D);
  Q = NR_ivector(1, T);

  /* copy data matrix into working memory */
  memcpy(&cont_data[1][1], mxGetPr(prhs[0]), T*D*sizeof(double));

  /* extract parameter values from Matlab parameter structure */
  memcpy(&(hmm.pi[1]), mxGetPr(mxGetField(prhs[1], 0, "pi")), 
         N*sizeof(double));
  memcpy(&(hmm.A[1][1]), mxGetPr(mxGetField(prhs[1], 0, "A")), 
         N*N*sizeof(double));
  memcpy(&(hmm.mu[1][1]), mxGetPr(mxGetField(prhs[1], 0, "mu")), 
         N*D*sizeof(double));
  for (i = 1; i <= N; i++)
  {
    memcpy(&(hmm.sigma[i][1][1]), 
           &(mxGetPr(mxGetField(prhs[1], 0, "sigma"))[(D*D)*(i-1)]), 
           N*D*sizeof(double));
  }

  if (nrhs > 2)
  {
    viterbi = (int) mxGetScalar(prhs[2]);
  }
  else
  {
    viterbi = UT_FALSE;
  }

  /* train the hmm */
  test_gauss_output_hmm(cont_data, hmm, Q, &L, T, viterbi, UT_FALSE);

  /* create the outputs */

  /* log likelihood */
  plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
  mxGetPr(plhs[0])[0] = L;
  
  /* state sequence */
  if (nlhs == 2)
  {
    plhs[1] = mxCreateDoubleMatrix(1, T, mxREAL);
    for (i = 1; i <= T; i++)
    {
      mxGetPr(plhs[1])[i-1] = (double) Q[i];
    }
  }

  /* free working memory */
  NR_free_dmatrix(cont_data, 1, T, 1, D);
  free_gauss_hmm(hmm);
  NR_free_ivector(Q, 1, T);
}
