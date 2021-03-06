/*******************************************************************************
MODULE HEADER:
rdahmm_arguments.h
*******************************************************************************/
#ifndef _RDAHMM_ARGUMENTS_H_
#define _RDAHMM_ARGUMENTS_H_ 1

#ifndef lint
static char rdahmm_arguments_hdr_rcsid[] = "$Id: rdahmm_arguments.h,v 1.1 2003/08/27 18:05:28 granat Exp $";
#endif
/* This string variable allows the RCS identification info to be printed. */

/* 
 * $Log: rdahmm_arguments.h,v $
 * Revision 1.1  2003/08/27 18:05:28  granat
 * Initial revision
 *
 * */

/*==============================================================================
Data Structures
==============================================================================*/

/*==============================================================================
Constants, Macros
==============================================================================*/

/*==============================================================================
Variables
==============================================================================*/

/*==============================================================================
Function Declarations
==============================================================================*/
int init_rdahmm_output_type_names (char **output_type_names);

int init_rdahmm_init_type_names (char **init_type_names);

int init_rdahmm_reg_type_names (char **init_reg_names);

int init_args(arg *args, char **output_type_names, char **init_type_names,
              char **reg_type_names);

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
                int *ntries, int *maxiters, int *seed, int *verbose);

#endif
