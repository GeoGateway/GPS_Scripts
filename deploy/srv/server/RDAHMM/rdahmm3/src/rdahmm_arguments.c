/*******************************************************************************
MODULE NAME
rdahmm_arguments

AUTHOR
Robert Granat

DESCRIPTION
See Robert Granat PhD thesis for details and notation.

COMPILE
make

NOTES
Not yet available

RG
*******************************************************************************/
#ifndef lint
static char rcsid[] = "$Id: rdahmm_arguments.c,v 1.1 2003/08/27 18:02:50 granat Exp $";
#endif
/* This string variable allows the RCS identification info to be printed. */

/* 
 * $Log: rdahmm_arguments.c,v $
 * Revision 1.1  2003/08/27 18:02:50  granat
 * Initial revision
 *
 * */

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
#include "rdahmm.h"

/* this program's header */
#include "rdahmm_arguments.h"

/*******************************************************************************
 INIT_RDAHMM_OUTPUT_TYPE_NAMES
 Initialize the name strings corresponding to the available output type options.
*******************************************************************************/
int init_rdahmm_output_type_names (char **output_type_names)
{
  output_type_names[RDAHMM_GAUSS_OUTPUT] = 
    (char*) strdup (RDAHMM_GAUSS_OUTPUT_NAME);

  return(UT_OK);
}

/*******************************************************************************
 INIT_RDAHMM_INIT_TYPE_NAMES
 Initialize the name strings corresponding to the available initialization type  options.
*******************************************************************************/
int init_rdahmm_init_type_names (char **init_type_names)
{
  init_type_names[RDAHMM_RANDOM_INIT] = 
    (char*) strdup (RDAHMM_RANDOM_INIT_NAME);

  init_type_names[RDAHMM_KMEANS_INIT] = 
    (char*) strdup (RDAHMM_KMEANS_INIT_NAME);

  return(UT_OK);
}

/*******************************************************************************
 INIT_RDAHMM_REG_TYPE_NAMES
 Initialize the name strings corresponding to the available regularization type 
 options.
*******************************************************************************/
int init_rdahmm_reg_type_names (char **init_reg_names)
{
  init_reg_names[RDAHMM_EUCLID_REG] = 
    (char*) strdup (RDAHMM_EUCLID_REG_NAME);

  return(UT_OK);
}

/*******************************************************************************
 INIT_ARGS
 Initialize structures that store command-line argument information.
*******************************************************************************/
int init_args(arg *args, char **output_type_names, char **init_type_names,
              char **reg_type_names)
{
  char output_type_choices[RDAHMM_MAX_STRING];
  char init_type_choices[RDAHMM_MAX_STRING];
  char reg_type_choices[RDAHMM_MAX_STRING];

  /* build strings containing names of all options, to be printed below */
  build_options_string(output_type_choices, output_type_names,
                       RDAHMM_NUM_OUTPUT_TYPES, 
                       "type of HMM output distribution");
  build_options_string(init_type_choices, init_type_names,
                       RDAHMM_NUM_INIT_TYPES, 
                       "type of HMM parameter initialization");
  build_options_string(reg_type_choices, reg_type_names,
                       RDAHMM_NUM_REG_TYPES, 
                       "type of regularized learning to apply");

  /* main input and output: files */
  Init_Arg(args, ARG_DATA_FILE, "data", FALSE, 1,
           "input observation sequence file");
  Init_Arg(args, ARG_L_FILE, "L", TRUE, 1,
           "output model log likelihood file");
  Init_Arg(args, ARG_Q_FILE, "Q", TRUE, 1,
           "output optimal state sequence file");
  Init_Arg(args, ARG_PI_FILE, "pi", TRUE, 1,
           "output model initial state probability file");
  Init_Arg(args, ARG_A_FILE, "A", TRUE, 1,
           "output model transition probability file");
  Init_Arg(args, ARG_B_FILE, "B", TRUE, 1,
           "output model output distribution file");
  Init_Arg(args, ARG_MINVAL_FILE, "minvalfile", TRUE, 1,
           "data minimum value file");
  Init_Arg(args, ARG_MAXVAL_FILE, "maxvalfile", TRUE, 1,
           "data maximum value file file");
  Init_Arg(args, ARG_RANGE_FILE, "rangefile", TRUE, 1,
           "data range file");
  Init_Arg(args, ARG_COVARS_WEIGHTS_FILE, "covarsweightsfile", TRUE, 1,
           "covariance component weightings file");
  Init_Arg(args, ARG_COVGRAPH_FILE, "covgraphfile", TRUE, 1,
           "covariance graph connectivity file");

  /* main input and output: variables */
  Init_Arg(args, ARG_T, "T", FALSE, 1,
           "number of observations");
  Init_Arg(args, ARG_D, "D", FALSE, 1,
           "dimension of observations");
  Init_Arg(args, ARG_N, "N", FALSE, 1,
           "number of model states");

  /* hmm output type */
  Init_Arg(args, ARG_OUTPUT_TYPE, "output_type", FALSE, 1,
           (char*) strdup(output_type_choices));

  /* hmm init type */
  Init_Arg(args, ARG_INIT_TYPE, "init_type", TRUE, 1,
           (char*) strdup(init_type_choices));

  /* hmm regularization type */
  Init_Arg(args, ARG_REG_TYPE, "reg_type", TRUE, 1,
           (char*) strdup(reg_type_choices));

  /* option flags */
  Init_Arg(args, ARG_VITERBI, "viterbi", TRUE, 0,
           "flag to perform Viterbi calculation of optimal state sequence");
  Init_Arg(args, ARG_REGULARIZE, "regularize", TRUE, 0,
           "flag to perform regularized learning");
  Init_Arg(args, ARG_ANNEAL, "anneal", TRUE, 0,
           "flag to perform deterministic annealing");
  Init_Arg(args, ARG_EVAL, "eval", TRUE, 0,
           "flag to perform state sequence evaluation only; model must exist");
  Init_Arg(args, ARG_ADDSTATE, "addstate", TRUE, 0,
           "flag to add blank state to existing model for evaluation");
  Init_Arg(args, ARG_WEIGHTCOVARS, "weightcovars", TRUE, 0,
           "flag to weight covariances during evaluation only mode");
  Init_Arg(args, ARG_COVGRAPH, "covgraph", TRUE, 0,
           "flag to use graph-based constraints on covariance");
  Init_Arg(args, ARG_VERBOSE, "v", TRUE, 0,
           "flag to print verbose output during execution");

  /* tuning parameters */
  Init_Arg(args, ARG_THRESH, "thresh", TRUE, 1,
           "convergence threshold");
  Init_Arg(args, ARG_PEPS, "peps", TRUE, 1,
           "perturbation epsilon");
  Init_Arg(args, ARG_OMEGA, "omega", TRUE, 4,
           "regularization weights (4)");
  Init_Arg(args, ARG_ANNEALSTEP, "annealstep", TRUE, 1,
           "annealing computational temperature step");
  Init_Arg(args, ARG_ANNEALFACTOR, "annealfactor", TRUE, 1,
           "annealing computational temperature multiplicative factor");
  Init_Arg(args, ARG_BETAMIN, "betamin", TRUE, 1,
           "starting computational temperature for annealing");
  Init_Arg(args, ARG_NTRIES, "ntries", TRUE, 1,
           "number of times to restart the maximization algorithm");
  Init_Arg(args, ARG_MAXITERS, "maxiters", TRUE, 1,
           "maximum number of EM iterations");
  Init_Arg(args, ARG_SEED, "seed", TRUE, 1,
           "seed for random number generation");

  return (UT_OK);
}

/*******************************************************************************
 HANDLE_ARGS
 Main function for parsing command line arguments and assigning variables
*******************************************************************************/
int handle_args(int argc, char **argv, arg *args, char **output_type_names, 
                char **init_type_names, char **reg_type_names, char **data_file,
                char **L_file, char **Q_file, char **pi_file, char **A_file,
                char **B_file, char **min_val_file, char **max_val_file, 
                char **range_file, int *T, int *D, int *N, int *output_type,
                int *init_type, double *thresh, double *peps, int *regularize,
                int *reg_type, double *omega_Q1, double *omega_Q2,
                double *omega_Q3, double *omega_sigma, int *anneal,
                double *anneal_step, double *anneal_factor, double *beta_min,
                int *eval_only, int *add_state, int *viterbi,
                int *weighted_covars, char **covars_weights_file, 
                int *use_covgraph, char **covgraph_file, 
                int *ntries, int *maxiters, int *seed, int *verbose)
{
  char *output_type_name;
  char *init_type_name;
  char *reg_type_name;
  char *base_name;
  char *strptr;

  /* initialize arrays of name strings corresponding to program options */
  init_rdahmm_output_type_names(output_type_names);
  init_rdahmm_init_type_names(init_type_names);
  init_rdahmm_reg_type_names(reg_type_names);

  /* initialize command line argument options */
  init_args(args, output_type_names, init_type_names, reg_type_names);

  /* get command line arguments */
  if (Get_Args(argc-1, argv, args, RDAHMM_NUM_ARGS) != Get_Arg_OK) 
  {
    return(UT_ERROR);
  }
  
  /* parse command line arguments and assign parameter values */
  if (args[ARG_DATA_FILE].exists)
  {
    *data_file = args[ARG_DATA_FILE].arg_str;
  }
  base_name = (char*) malloc( strlen(*data_file) * sizeof(char) );
  get_base_name(*data_file, base_name);

  if (args[ARG_L_FILE].exists)
  {
    *L_file = args[ARG_L_FILE].arg_str;
  }
  else
  {
    *L_file = make_extended_name(base_name, ".L");
  }
  
  if (args[ARG_Q_FILE].exists)
  {
    *Q_file = args[ARG_Q_FILE].arg_str;
  }
  else
  {
    *Q_file = make_extended_name(base_name, ".Q");
  }

  if (args[ARG_PI_FILE].exists)
  {
    *pi_file = args[ARG_PI_FILE].arg_str;
  }
  else
  {
    *pi_file = make_extended_name(base_name, ".pi");
  }

  if (args[ARG_A_FILE].exists)
  {
    *A_file = args[ARG_A_FILE].arg_str;
  }
  else
  {
    *A_file = make_extended_name(base_name, ".A");
  }

  if (args[ARG_B_FILE].exists)
  {
    *B_file = args[ARG_B_FILE].arg_str;
  }
  else
  {
    *B_file = make_extended_name(base_name, ".B");
  }

  if (args[ARG_MINVAL_FILE].exists)
  {
    *min_val_file = args[ARG_MINVAL_FILE].arg_str;
  }
  else
  {
    *min_val_file = make_extended_name(base_name, ".minval");
  }

  if (args[ARG_MAXVAL_FILE].exists)
  {
    *max_val_file = args[ARG_MAXVAL_FILE].arg_str;
  }
  else
  {
    *max_val_file = make_extended_name(base_name, ".maxval");
  }

  if (args[ARG_RANGE_FILE].exists)
  {
    *range_file = args[ARG_RANGE_FILE].arg_str;
  }
  else
  {
    *range_file = make_extended_name(base_name, ".range");
  }

  if (args[ARG_COVARS_WEIGHTS_FILE].exists)
  {
    *covars_weights_file = args[ARG_COVARS_WEIGHTS_FILE].arg_str;
  }
  else
  {
    *covars_weights_file = make_extended_name(base_name, ".covarsweights");
  }

  if (args[ARG_COVGRAPH_FILE].exists)
  {
    *covgraph_file = args[ARG_COVGRAPH_FILE].arg_str;
  }
  else
  {
    *covgraph_file = make_extended_name(base_name, ".covgraph");
  }

  if (args[ARG_T].exists)
  {
    *T = atoi(args[ARG_T].arg_str);
  }

  if (args[ARG_D].exists)
  {
    *D = atoi(args[ARG_D].arg_str);
  }

  if (args[ARG_N].exists)
  {
    *N = atoi(args[ARG_N].arg_str);
  }

  if (args[ARG_OUTPUT_TYPE].exists)
  {
    output_type_name = args[ARG_OUTPUT_TYPE].arg_str;
    *output_type = match_option_name(output_type_name, output_type_names,
                                     RDAHMM_NUM_OUTPUT_TYPES, "output type");
  }

  if (args[ARG_INIT_TYPE].exists)
  {
    init_type_name = args[ARG_INIT_TYPE].arg_str;
    *init_type = match_option_name(init_type_name, init_type_names,
                                   RDAHMM_NUM_INIT_TYPES, "init type");
  }
  else
  {
    *init_type = RDAHMM_RANDOM_INIT;
  }

  if (args[ARG_THRESH].exists)
  {
    *thresh = atof(args[ARG_THRESH].arg_str);
  }
  else
  {
    *thresh = ARG_THRESH_DEFAULT;
  }

  if (args[ARG_PEPS].exists)
  {
    *peps = strtod(args[ARG_PEPS].arg_str, (char **) NULL);
  }
  else
  {
    *peps = ARG_PEPS_DEFAULT;
  }

  if (args[ARG_REGULARIZE].exists)
  {
    *regularize = UT_TRUE;
  }
  else
  {
    *regularize = UT_FALSE;
  }

  if (args[ARG_REG_TYPE].exists)
  {
    reg_type_name = args[ARG_REG_TYPE].arg_str;
    *reg_type = match_option_name(reg_type_name, reg_type_names,
                                  RDAHMM_NUM_REG_TYPES, "init type");
  }
  else
  {
    if (*output_type == RDAHMM_GAUSS_OUTPUT)
    {
      *reg_type = RDAHMM_EUCLID_REG;
    }
  }

  if (args[ARG_OMEGA].exists)
  {
    *omega_Q1 = strtod(args[ARG_OMEGA].arg_str, &strptr);
    args[ARG_OMEGA].arg_str = strptr;
    *omega_Q2 = strtod(args[ARG_OMEGA].arg_str, &strptr);
    args[ARG_OMEGA].arg_str = strptr;
    *omega_Q3 = strtod(args[ARG_OMEGA].arg_str, &strptr);
    args[ARG_OMEGA].arg_str = strptr;
    *omega_sigma = strtod(args[ARG_OMEGA].arg_str, (char **) NULL);
  }
  else
  {
    *omega_Q1 = ARG_OMEGA_Q1_DEFAULT;
    *omega_Q2 = ARG_OMEGA_Q2_DEFAULT;
    *omega_Q3 = ARG_OMEGA_Q3_DEFAULT;
    *omega_sigma = ARG_OMEGA_SIGMA_DEFAULT;
  }

  if (args[ARG_ANNEAL].exists)
  {
    *anneal = UT_TRUE;
  }
  else
  {
    *anneal = UT_FALSE;
  }

  if (args[ARG_ANNEALSTEP].exists)
  {
    *anneal_step = strtod(args[ARG_ANNEALSTEP].arg_str, (char **) NULL);
  }
  else
  {
    *anneal_step = ARG_ANNEALSTEP_DEFAULT;
  }

  if (args[ARG_ANNEALFACTOR].exists)
  {
    *anneal_factor = strtod(args[ARG_ANNEALFACTOR].arg_str, (char **) NULL);
  }
  else
  {
    *anneal_factor = ARG_ANNEALFACTOR_DEFAULT;
  }

  if (args[ARG_BETAMIN].exists)
  {
    *beta_min = strtod(args[ARG_BETAMIN].arg_str, (char **) NULL);
  }
  else
  {
    *beta_min = ARG_BETAMIN_DEFAULT;
  }

  if (args[ARG_EVAL].exists)
  {
    *eval_only = UT_TRUE;
  }
  else
  {
    *eval_only = UT_FALSE;
  }

  if (args[ARG_VITERBI].exists)
  {
    *viterbi =  UT_TRUE;
  }
  else
  {
    *viterbi = UT_FALSE;
  }

  if (args[ARG_ADDSTATE].exists)
  {
    *add_state =  UT_TRUE;
  }
  else
  {
    *add_state = UT_FALSE;
  }

  if (args[ARG_WEIGHTCOVARS].exists)
  {
    *weighted_covars =  UT_TRUE;
  }
  else
  {
    *weighted_covars = UT_FALSE;
  }

  if (args[ARG_COVGRAPH].exists)
  {
    *use_covgraph =  UT_TRUE;
  }
  else
  {
    *use_covgraph = UT_FALSE;
  }

  if (args[ARG_NTRIES].exists)
  {
    *ntries = atoi(args[ARG_NTRIES].arg_str);
  }
  else
  {
    *ntries = 1;
  }

  if (args[ARG_MAXITERS].exists)
  {
    *maxiters = atoi(args[ARG_MAXITERS].arg_str);
  }
  else
  {
    *maxiters = ARG_MAXITERS_DEFAULT;
  }

  if (args[ARG_SEED].exists)
  {
    *seed = atoi(args[ARG_SEED].arg_str);
  }
  else
  {
    *seed = UT_FALSE;
  }

  if (args[ARG_VERBOSE].exists)
  {
    *verbose = UT_TRUE;
  }
  else
  {
    *verbose = UT_FALSE;
  }

  return(UT_OK);
}
