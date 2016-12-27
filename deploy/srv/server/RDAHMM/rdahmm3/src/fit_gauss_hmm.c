/*******************************************************************************
MODULE NAME
fit_hmm

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
  int               N;
  RDAEM_Parameters  params;
  int               viterbi;

  /* outputs */
  double       L;
  Gaussian_HMM hmm;
  int          *Q;

  /* problem dimensions */
  int     T;
  int     D;

  /* other variables */
  int    i;
  double *p;

  /* check for proper number of arguments */
  if (nrhs != 3) 
  {
    mexErrMsgIdAndTxt("MyToolbox:fit_gauss_hmm:nrhs", 
                      "Three inputs required.");
  }
  if (nlhs < 1) 
  {
    mexErrMsgIdAndTxt("MyToolbox:fit_gauss_hmm:nlhs", 
                      "At least one output required.");
  }
  if (nlhs > 6)
  {
    mexErrMsgIdAndTxt("MyToolbox:fit_gauss_hmm:nlhs", 
                      "Too many outputs (max 6).");
  }

  /* get array dimensions */
  T = (int) mxGetN(prhs[0]);
  D = (int) mxGetM(prhs[0]);
  N = (int) mxGetScalar(prhs[1]);

  /* allocate working memory */
  cont_data = NR_dmatrix(1, T, 1, D);
  allocate_gauss_hmm(&hmm, N, D);
  Q = NR_ivector(1, T);
  params.covgraph = NR_ivector(1, N);

  /* copy data matrix into working memory */
  memcpy(&cont_data[1][1], mxGetPr(prhs[0]), T*D*sizeof(double));

  /* extract parameter values from Matlab parameter structure */
  params.thresh = mxGetScalar(mxGetField(prhs[2], 0, "thresh"));
  params.peps = mxGetScalar(mxGetField(prhs[2], 0, "peps"));
  params.regularize = (int) mxGetScalar(mxGetField(prhs[2], 0, "regularize"));
  params.reg_type = (int) mxGetScalar(mxGetField(prhs[2], 0, "reg_type"));
  params.init_type = (int) mxGetScalar(mxGetField(prhs[2], 0, "init_type"));
  params.omega_Q1 = mxGetScalar(mxGetField(prhs[2], 0, "omega_Q1"));
  params.omega_Q2 = mxGetScalar(mxGetField(prhs[2], 0, "omega_Q2"));
  params.omega_Q3 = mxGetScalar(mxGetField(prhs[2], 0, "omega_Q3"));
  params.omega_sigma = mxGetScalar(mxGetField(prhs[2], 0, "omega_sigma"));
  params.anneal = mxGetScalar(mxGetField(prhs[2], 0, "anneal"));
  params.anneal_step = mxGetScalar(mxGetField(prhs[2], 0, "anneal_step"));
  params.anneal_factor = mxGetScalar(mxGetField(prhs[2], 0, "anneal_factor"));
  params.beta_min = mxGetScalar(mxGetField(prhs[2], 0, "beta_min"));
  params.ntries = (int) mxGetScalar(mxGetField(prhs[2], 0, "ntries"));
  params.maxiters = (int) mxGetScalar(mxGetField(prhs[2], 0, "maxiters"));
  params.seed = (int) mxGetScalar(mxGetField(prhs[2], 0, "seed"));
  for (i = 1; i <= N; i++)
  {
    p = mxGetPr(mxGetField(prhs[2], 0, "covgraph"));
    params.covgraph[i] = (int) p[i-1];
  }
  viterbi = (int) mxGetScalar(mxGetField(prhs[2], 0, "viterbi"));

  /* train the hmm */
  train_gauss_output_hmm(cont_data, hmm, params, Q, &L, T, viterbi, UT_FALSE);

  /* create the outputs */

  /* log likelihood */
  plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
  mxGetPr(plhs[0])[0] = L;
  
  /* model parameters */
  if (nlhs >= 2)
  {
    plhs[1] = mxCreateDoubleMatrix(1, N, mxREAL);
    memcpy(mxGetPr(plhs[1]), &(hmm.pi[1]), N*sizeof(double));
  }
  if (nlhs >= 3)
  {
    plhs[2] = mxCreateDoubleMatrix(N, N, mxREAL);
    memcpy(mxGetPr(plhs[2]), &(hmm.A[1][1]), N*N*sizeof(double));
  }
  if (nlhs >= 4)
  {
    plhs[3] = mxCreateDoubleMatrix(D, N, mxREAL);
    memcpy(mxGetPr(plhs[3]), &(hmm.mu[1][1]), N*D*sizeof(double));
  }
  if (nlhs >= 5)
  {
    plhs[4] = mxCreateDoubleMatrix(D, D*N, mxREAL);
    for (i = 1; i <= N; i++)
    {
      memcpy(&(mxGetPr(plhs[4])[(D*D)*(i-1)]), 
             &(hmm.sigma[i][1][1]), D*D*sizeof(double));
    }
  }

  /* state sequence */
  if (nlhs == 6)
  {
    plhs[5] = mxCreateDoubleMatrix(1, T, mxREAL);
    for (i = 1; i <= T; i++)
    {
      mxGetPr(plhs[5])[i-1] = (double) Q[i];
    }
  }

  /* free working memory */
  NR_free_dmatrix(cont_data, 1, T, 1, D);
  free_rdaem_params(params, hmm);
  free_gauss_hmm(hmm);
  NR_free_ivector(Q, 1, T);
}
