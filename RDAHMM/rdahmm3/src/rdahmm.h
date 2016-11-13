/*******************************************************************************
PROGRAM HEADER:
rdahmm.h
*******************************************************************************/
#ifndef _RDAHMM_H_
#define _RDAHMM_H_ 1

#ifndef lint
static char rdahmm_hdr_rcsid[] = "$Id: rdahmm.h,v 1.1 2003/08/27 18:05:02 granat Exp granat $";
#endif
/* This string variable allows the RCS identification info to be printed. */

/* 
 * $Log: rdahmm.h,v $
 * Revision 1.1  2003/08/27 18:05:02  granat
 * Initial revision
 *
 * */

/*==============================================================================
Data Structures
==============================================================================*/

/*==============================================================================
Constants, Macros
==============================================================================*/

/* arguments */

#define ARG_DATA_FILE      1
#define ARG_L_FILE         2
#define ARG_Q_FILE         3
#define ARG_PI_FILE        4
#define ARG_A_FILE         5
#define ARG_B_FILE         6
#define ARG_MINVAL_FILE    7
#define ARG_MAXVAL_FILE    8
#define ARG_RANGE_FILE     9
#define ARG_COVARS_WEIGHTS_FILE 10
#define ARG_COVGRAPH_FILE 11

#define ARG_T             12
#define ARG_D             13
#define ARG_N             14
#define ARG_OUTPUT_TYPE   15
#define ARG_INIT_TYPE     16
#define ARG_THRESH        17
#define ARG_PEPS          18
#define ARG_REGULARIZE    19
#define ARG_REG_TYPE      20
#define ARG_OMEGA         21
#define ARG_ANNEAL        22
#define ARG_ANNEALSTEP    23
#define ARG_ANNEALFACTOR  24
#define ARG_BETAMIN       25

#define ARG_EVAL          26
#define ARG_ADDSTATE      27
#define ARG_WEIGHTCOVARS  28
#define ARG_COVGRAPH      29
#define ARG_NTRIES        30
#define ARG_MAXITERS      31

#define ARG_SEED          32

#define ARG_VERBOSE       33

#define ARG_VITERBI       34

#define RDAHMM_NUM_ARGS   34

#define RDAHMM_NUM_OUTPUT_TYPES 1
#define RDAHMM_NUM_INIT_TYPES 2
#define RDAHMM_NUM_REG_TYPES 1

#define RDAHMM_GAUSS_OUTPUT 1
#define RDAHMM_GAUSS_OUTPUT_NAME "gauss"

#define RDAHMM_RANDOM_INIT 1
#define RDAHMM_RANDOM_INIT_NAME "random"
#define RDAHMM_KMEANS_INIT 2
#define RDAHMM_KMEANS_INIT_NAME "kmeans"

#define RDAHMM_EUCLID_REG 1
#define RDAHMM_EUCLID_REG_NAME "euclid"

#define ARG_THRESH_DEFAULT      1.0e-3
#define ARG_PEPS_DEFAULT        1.0e-3
#define ARG_OMEGA_Q1_DEFAULT    1.0e-6
#define ARG_OMEGA_Q2_DEFAULT    1.0e-6
#define ARG_OMEGA_Q3_DEFAULT    1.0e-6
#define ARG_OMEGA_SIGMA_DEFAULT 1.0e-6

#define ARG_ANNEALSTEP_DEFAULT   0.05
#define ARG_ANNEALFACTOR_DEFAULT 0
#define ARG_BETAMIN_DEFAULT      0

#define ARG_MAXITERS_DEFAULT     1000

#define RDAHMM_MAX_STRING 100
/*==============================================================================
Variables
==============================================================================*/

/*==============================================================================
Function Declarations
==============================================================================*/

#endif
