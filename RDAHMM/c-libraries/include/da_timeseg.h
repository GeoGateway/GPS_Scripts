/*******************************************************************************
MODULE HEADER:
da_timeseg.h
*******************************************************************************/

#ifndef _DA_TIMESEG_H_
#define _DA_TIMESEG_H_
/* Protects from multiple inclusion. */

#ifndef lint
static char da_timeseg_h_rcsid[] = "$Id: da_timeseg.h,v 1.4 1999/07/20 17:34:48 granat Exp $";
#endif
/* This string variable allows the RCS identification info to be printed. */

/* 
 * $Log: da_timeseg.h,v $
 * Revision 1.4  1999/07/20 17:34:48  granat
 * changed to double precision
 *
 * Revision 1.3  1997/03/06 02:25:31  agray
 * updated name of init_hmm_model_type_names().
 *
 * Revision 1.2  1997/01/29 21:41:01  agray
 * moved initialization-type variables to hmm driver.
 * added model-type names.
 * new format.
 *
 * Revision 1.1  1996/10/31 02:22:30  agray
 * Initial revision
 *
 * */

/*==============================================================================
Data Structures
==============================================================================*/

/* HMM */

struct da_hmm_struct {
  int       model_type;   /* discrete, single gaussian, or mixture */
  int       num_states;   /* N, number of states */
  int       num_symbols;  /* M, number of discrete output symbols */
  int       num_comps;    /* K, number of continuous mixture components */
  int       num_dims;     /* number of dimensions/features in data */
  double**   prob_trans;   /* A, state transition probability matrix */
  double**   prob_symbol;  /* B, output symbol probabilities for each state */
  double*    prob_init;    /* pi, initial probability of being in each state */
  double***  means;        /* mu's, means of continuous mixture model */
  double**** covars;       /* U's, covariance matrices of cont. mixture model */
  double**   weights;      /* c's, mixing proportions of cont. mixture model */
};

typedef struct da_hmm_struct HMM;

/* maybe later:
  def*    def;          * structure containing attributes' type & name info *
#define name_of_attr(H,I)     ((H->def)[I]) */

/*==============================================================================
Constants, Macros
==============================================================================*/

/* model types */

#define HMM_DISCRETE_MODEL   1
#define HMM_GAUSSIAN_MODEL   2
#define HMM_MIXTURE_MODEL    3

#define HMM_NUM_MODEL_TYPES  3

#define HMM_DISCRETE_MODEL_NAME "discrete"
#define HMM_GAUSSIAN_MODEL_NAME "gaussian"
#define HMM_MIXTURE_MODEL_NAME  "mixture"

/*==============================================================================
Variables
==============================================================================*/

/*==============================================================================
Function Declarations
==============================================================================*/

/* general HMM utilities */

HMM *hmm(int model_type, int num_states, int num_symbols, int num_comps,
         int num_dims);
void free_hmm(HMM *the_hmm);


/* learning/scoring discrete HMM's */

double learn_disc_hmm_using_em(HMM *the_hmm, HMM *last_hmm,
                              int *data, int num_data, 
                              double **alpha, double **beta, double *scale,
                              double **sum_t_xi, double **sum_t_gamma_by_symbol,
                              double *sum_t_gamma, double **zeta_t, 
                              double *gamma_t, double threshold);
double compute_disc_hmm_likelihood(HMM *the_hmm, int *data, int num_data, 
                                  double **alpha, double **beta, double *scale);
int estimate_disc_hmm_params(HMM *old_hmm, HMM *new_hmm,
                             int *data, int num_data, 
                             double **alpha, double **beta, 
                             double **sum_t_xi, double **sum_t_gamma_by_symbol,
                             double *sum_t_gamma, double **zeta_t,double *gamma_t);

/* learning/scoring continuous HMM's */

double learn_cont_hmm_using_em(HMM *the_hmm, HMM *last_hmm, 
                              double **data, int num_data,
                              double **alpha, double **beta, double *scale,
                              double **sum_t_xi, double **sum_t_gamma_by_comp,
                              double *sum_t_gamma, double **zeta_t,double *gamma_t,
                              double *mixture_densities, double **comp_densities, 
                              double *diff_vec, double threshold, double min_diag);
double compute_cont_hmm_likelihood(HMM *the_hmm, double **data, int num_data, 
                                  double **alpha, double **beta, double *scale,
                                  double min_diag);
int estimate_cont_hmm_params(HMM *old_hmm, HMM *new_hmm, 
                             double **data, int num_data, 
                             double **alpha, double **beta, 
                             double **sum_t_xi, double **sum_t_gamma_by_comp,
                             double *sum_t_gamma, double **zeta_t, double *gamma_t,
                             double *mixture_densities, double **comp_densities, 
                             double *diff, double min_diag);

/* simulating data from HMM's */

int simulate_data_using_disc_hmm(HMM* the_hmm, int *data, int *states,
                                 int num_data, double *cum_prob_states, 
                                 double *cum_prob_symbols, int seed);
int simulate_data_using_cont_hmm(HMM* the_hmm, double **data, int *states,
                                 int num_data, double *cum_prob_states,
                                 double *cum_prob_comps, double **chol_factor,
                                 double *diag, int seed);

/* state assignment */

int assign_states_using_iml(int *states, double **alpha, double **beta, 
                            int num_states, int num_data);

/* predefined names */

int init_hmm_model_type_names(char **model_type_names);

#endif /* _DA_TIMESEG_H_ */
