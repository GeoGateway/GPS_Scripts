/*******************************************************************************
MODULE HEADER:
rdahmm_setup.h
*******************************************************************************/
#ifndef _RDAHMM_SETUP_H_
#define _RDAHMM_SETUP_H_ 1

#ifndef lint
static char rdahmm_setup_hdr_rcsid[] = "$Id$";
#endif
/* This string variable allows the RCS identification info to be printed. */

/* 
 * $Log$
 * */

/*==============================================================================
Data Structures
==============================================================================*/
typedef struct {
  int      N;        /* number of states */
  int      D;        /* dimension of data space */
  double   **A;      /* state-to-state transition probability matrix */
  double   *pi;      /* initial state probability vector */
  double   **mu;     /* state means */
  double   ***sigma; /* state covariance matrices */
} Gaussian_HMM;

typedef struct {
  double  thresh;          /* convergence threshold */
  double  peps;            /* perturbation epsilon */
  int     regularize;      /* flag to perform regularized learning */
  int     reg_type;        /* type of regularization term to use */
  int     init_type;       /* type of parameter initialization to use */
  double  omega_Q1;        /* regularization weighting term */
  double  omega_Q2;        /* regularization weighting term */
  double  omega_Q3;        /* regularization weighting term */
  double  omega_sigma;     /* regularization weighting term (Gaussian output) */
  int     anneal;          /* flag to perform deterministic annealing */
  double  anneal_step;     /* temperature step for deterministic annealing */
  double  anneal_factor;   /* temperature factor for deterministic annealing */
  double  beta_min;        /* initial computational temperature for DA */
  int     *covgraph;       /* covariance matrix constraints (Gaussian output) */
  int     ntries;          /* number of times to retry the optimization */
  int     maxiters;        /* maximum number of EM iterations */
  int     seed;            /* seed for initializing random number generator */
} RDAEM_Parameters;

/*==============================================================================
Constants, Macros
==============================================================================*/

/*==============================================================================
Variables
==============================================================================*/

/*==============================================================================
Function Declarations
==============================================================================*/
int handle_rdahmm_input(int argc, char **argv, arg *args,
                        char **output_type_names, char **init_type_names,
                        char **reg_type_names, char **data_file, char **L_file, 
                        char **Q_file, char **pi_file, char **A_file,
                        char **B_file, char **min_val_file, char **max_val_file,
                        char **range_file, int *T, int *D, int *N,
                        int *output_type, int *init_type, double *thresh,
                        double *peps, int *regularize, int *reg_type,
                        double *omega_Q1, double *omega_Q2, double *omega_Q3,
                        double *omega_sigma, int *anneal, double *anneal_step,
                        double *anneal_factor, double *beta_min, int *eval_only,
                        int *add_state, int *viterbi, int *weighted_covars,
                        char **covars_weights_file, int *use_covgraph,
                        char **covgraph_file, int *ntries, int *maxiters,
                        int *seed, int *verbose,
                        double ***cont_data, double **min_val, double **max_val,
                        double **range);

int allocate_gauss_hmm(Gaussian_HMM *hmm, int N, int D);

int setup_gauss_hmm(Gaussian_HMM *hmm,
                    int add_state, int eval_only, int weighted_covars,
                    int use_covgraph,
                    double *min_val, double *range,
                    double **cont_data, int N, int D, int T,
                    char *pi_file, char *A_file, char *B_file,
                    char *covars_weights_file, char *covgraph_file,
                    int **covgraph);

int setup_rdaem_params(RDAEM_Parameters *params, double thresh, double peps,
                       int regularize, int init_type, int reg_type, 
                       double omega_Q1, double omega_Q2, double omega_Q3, 
                       double omega_sigma, int anneal, double anneal_step, 
                       double anneal_factor, double beta_min, int *covgraph, 
                       int ntries, int maxiters, int seed);

int rdahmm_setup(int argc, char **argv, arg *args,
                 char **output_type_names, char **init_type_names,
                 char **reg_type_names, char **data_file, char **L_file,
                 char **Q_file, char **pi_file, char **A_file,
                 char **B_file, char **min_val_file, char **max_val_file,
                 char **range_file, int *T, int *D, int *N,
                 int *output_type, int *init_type, double *thresh,
                 double *peps, int *regularize, int *reg_type,
                 double *omega_Q1, double *omega_Q2, double *omega_Q3,
                 double *omega_sigma, int *anneal, double *anneal_step,
                 double *anneal_factor, double *beta_min, int *eval_only,
                 int *add_state, int *viterbi, int *weighted_covars,
                 char **covars_weights_file, int *use_covgraph,
                 char **covgraph_file, int **covgraph, int *ntries,
                 int *maxiters, int *seed, int *verbose,
                 double ***cont_data, double **min_val, double **max_val,
                 double **range,
                 Gaussian_HMM *hmm,
                 RDAEM_Parameters *params);

int free_gauss_hmm(Gaussian_HMM hmm);

int free_rdaem_params(RDAEM_Parameters params, Gaussian_HMM hmm);

#endif
