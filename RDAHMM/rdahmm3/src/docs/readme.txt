RDAHMM - HMM based time series analysis software
Robert Granat
granat@jpl.nasa.gov

RDAHMM fits a hidden Markov model (HMM) to time series data.
As input, it takes the time series (can be multidimensional),
model size (number of hidden states), and various algorithm
parameters.  As output, it produces the fitted model, the
log likelihood of that model with respect to the input time 
series, and a labeling of the time series observations
according to the estimated (individually most likely) hidden 
state assignments.

RDAHMM implements the regularized deterministic annealing 
expectation-maximization (RDAEM) algorithm for model fitting.
However, it can be run using only standard expectation-
maximization (EM) or deterministic annealing EM.  It is
possible to run it using only regularized EM, but this is
not recommended.  Detailed documentation of the mathematics
implemented can be found in [Granat04].

The method can be run in one of two modes. The default mode
is one in which the method both learns a model for the data 
and classifies it simultaneously.  Alternatively, if a
model has already been learned and stored in files using 
the default mode, that model can then be used to evaluate
other time series with the same number of feature dimensions.
This is referred to as the "evaluation mode."

Inputs and outputs to RDAHMM are flat ascii matrices and
vectors.  The method has numerous control parameters for
expert users but the default values should be suitable for
most problems and only a few command line arguments are 
mandatory.

RDAHMM is run from the command line as follows:

USAGE: ../bin/rdahmm  -data 'input observation sequence file'
     [-L 'output model log likelihood file']
     [-Q 'output optimal state sequence file']
     [-pi 'output model initial state probability file']
     [-A 'output model transition probability file']
     [-B 'output model output distribution file']
     [-minvalfile 'data minimum value file']
     [-maxvalfile 'data maximum value file file']
     [-rangefile 'data range file']
     [-covarsweightsfile 'covariance component weightings file']
     -T 'number of observations'
     -D 'dimension of observations'
     -N 'number of model states'
     -output_type 'type of HMM output distribution {gauss}'
     [-init_type 'type of HMM parameter initialization {random}']
     [-thresh 'convergence threshold']
     [-peps 'perturbation epsilon']
     [-regularize 'flag to perform regularized learning']
     [-reg_type 'type of regularized learning to apply {euclid}']
     [-omega 'regularization weights (4)']
     [-anneal 'flag to perform deterministic annealing']
     [-annealstep 'annealing computational temperature step']
     [-annealfactor 'annealing computational temperature multiplicative factor']
     [-betamin 'starting computational temperature for annealing']
     [-eval 'flag to perform state sequence evaluation only; model must exist']
     [-addstate 'flag to add blank state to existing model for evaluation']
     [-weightcovars 'flag to weight covariances during evaluation only mode']
     [-ntries 'number of times to restart the maximization algorithm']
     [-seed 'seed for random number generation']
     [-v 'flag to print verbose output during execution']

In the following section, these arguments will be discussed in detail.

     -data 'input observation sequence file'

The input data is expected to be in the form of space-delimited ascii
TxD matrix of time series observations, where T is the number of time
steps/observations, and D is the dimension of the data.  Currently 
only continuous-valued data is supported.  The data file name can be
anything, but the convention "filename.dat" is preferred.

     [-L 'output model log likelihood file']

This can be anything, but if it is not provided, and the data file
is called "filename.dat", this will default to "filename.L".  The
file contains a single asciii floating point value.

     [-Q 'output optimal state sequence file']

This can be anything, but if it is not provided, and the data file
is called "filename.dat", this will default to "filename.Q".  The
file contains a column vector of ascii integer values of length T.
These are the observation labels.

     [-pi 'output model initial state probability file']

This can be anything, but if it is not provided, and the data file
is called "filename.dat", this will default to "filename.pi".  The
file contains a column vector of ascii floating point values of
length N, where N is the number of model states.  This the 
probability distribution of the initial state at t=1.

     [-A 'output model transition probability file']

This can be anything, but if it is not provided, and the data file
is called "filename.dat", this will default to "filename.A".  The
file contains an NxN matrix of space delimited floating point
values.  The values in row i represents the probability 
distribution of the next state given that the current state is i. 

     [-B 'output model output distribution file']

This can be anything, but if it is not provided, and the data file
is called "filename.dat", this will default to "filename.B"  The
format of this file type will change depending on the type of 
output distribution.  Currently, however, only Gaussian output
distributions are supported.

Gaussian: The file contains a ((D+1)*N) x D matrix of space
delimited ascii floating point values.  These can be thought of
as N (D+1)xD blocks, each describing the output distribution for
a state.  The first row of the block is the mean for that state,
which the lower DxD matrix is the covariance matrix for that
state.

     [-minvalfile 'data minimum value file']

This can be anything, but if it is not provided, and the data file
is called "filename.dat", this will default to "filename.minval"  
The output file contains an ascii vector of space separated 
floating point values indicating the minimum values of each input
data dimension.

     [-maxvalfile 'data maximum value file file']

This can be anything, but if it is not provided, and the data file
is called "filename.dat", this will default to "filename.maxval"  
The output file contains an ascii vector of space separated 
floating point values indicating the maximum values of each input
data dimension.

     [-rangefile 'data range file']

This can be anything, but if it is not provided, and the data file
is called "filename.dat", this will default to "filename.range"  
The output file contains an ascii vector of space separated 
floating point values indicating the range of values (maximum
minus minimum) of each input data dimension.

     [-covarsweightsfile 'covariance component weightings file']

The file containing a vector of floating point values of length
equal to the number of feature dimensions of the data.  These
are used to weight the covariance matrices in the evaluation
mode ("-eval").  This argument only has effect if the
"-weightcovars" flag is used.

     -T 'number of observations'

The number of rows in the time series file.

     -D 'dimension of observations'

The number of columns in the time series file.

     -N 'number of model states'

The number of HMM states, or classification catagories, used to 
analyze the data.  Setting N=1 will not produce useful results.

     -output_type 'type of HMM output distribution {gauss}'

Currently only Gaussian outputs are supported.

     [-init_type 'type of HMM parameter initialization {random}']

Currently only random initialization is supported.

     [-thresh 'convergence threshold']

EM convergence threshold, used for all EM variations.  Defaults to
1.0e-3.  If the objective function has very flat plateaus, it may 
be useful to set this lower.  However, lower thresholds can reduce 
time to convergence considerably.

     [-peps 'perturbation epsilon']

The value by which the output distributions are perturbed if 
identical output distributions are found in two different states.
For Gaussian outputs, the mean of one is perturbed.  Defaults to
1.0e-3.  Do not change this value unless you have an intimate 
understanding of the algorithm and code.

     [-regularize 'flag to perform regularized learning']

Perform regularized EM.

     [-reg_type 'type of regularized learning to apply {euclid}']

Only Euclidean regularization is currently supported.  See 
[Granat04].

     [-omega 'regularization weights (4)']

The four regularization weights to be used with regularization.
These are omega_Q1, omega_Q2, omega_Q3, and omega_sigma in 
[Granat04].  Even if no other regularization is to be performed,
it is highly recommended that the weights be set to 0 0 0 1.0e-6,
as this will prevent computational errors due to inversion of
poorly conditioned matrices.  To do only euclidean regularization 
with the maximum weight at each iteration, set the weights to 0 0 
[large value] 1.0e-6.

     [-anneal 'flag to perform deterministic annealing']

Perform deterministic annealing EM.

     [-annealstep 'annealing computational temperature step']

The change in computational temperature at each annealing outer
loop.  The current default is 0.05, but 0.001 is safer and usually
returns better results, although it is slower.

     [-annealfactor 'annealing computational temperature multiplicative factor']

An alternative to "-annealstep," if specified it increases the
computational temperature by a multiplicative factor each time
through the outer loop.  This argument exists for testing and
demonstration purposes only; at this time experimental results
indicate that it returns inferior results to using a fixed step
increase in the computational temperature.  Use of this argument
is not recommended for most users.

     [-betamin 'starting computational temperature for annealing']

The starting computational temperature when performing annealing.
The value must be greater than zero if "-annealfactor" is used.
The use of this argument is not recommended for most users.

     [-eval 'flag to perform state sequence evaluation only; model must exist']

Run the program in an evaluation mode.  In this mode, existing
files output by a previous run of the program are used as input.
The necessary inputs are the model parameters, contained in ".pi", 
".A", and ".B" files, as well as ".minval", ".maxval", and ".range"
files describing the original training data.  The program uses
this information to evaluate the time series indicated by the 
"-data" argument and produces as its output a state sequence for
that time series classifying each output according to one of the
discrete model states.  The "-Q" argument is used to provide the
output file name for this information.

     [-addstate 'flag to add blank state to existing model for evaluation']

This argument is only valid in evaluation only mode ("-eval").  
It specifies that an extra state be added to the model used in
evaluation.  The purpose of this model is to account for observed
values in the data being evaluated that lie well outside of those
described by the states in the evaluation model.  The added state
is currently only supported for Gaussian outputs; it has a mean
and covariance of that of the entire time series.  The initial
and state-to-state probabilities for the added state are 1/(N+1),
where N is the number of states in the model being used for
evaluation.

     [-weightcovars 'flag to weight covariances during evaluation only mode']

When using the program in evaluation mode ("-eval"), this 
argument indicates the weights on the varies feature dimensions
are to be used to alter the covariance matrices of the states 
in order to emphasize or de-emphasize certain dimensions of the 
feature space.  The effects of using these weights are still not 
well tested and usage of this argument is not recommended for most
users.  Actual weight values are provided in the file given by
"-covarsweightsfile".

     [-ntries 'number of times to restart the maximization algorithm']

This argument instructs the program to repeat the optimization
procedure the input integer number of times, taking as the final
result the run with the greatest log likelihood.  If unspecified,
defaults to one.

     [-seed 'seed for random number generation']

Set this to an integer value if you don't want to use a random
seed.  This allows for replication of results.

     [-v 'flag to print verbose output during execution']

If set, prints information about the algorithm progress to stdout.

The following is an example of running the program on data from
the space station beta gimbal array motor current sensor.  First,
in default mode:

rdahmm -data P6PB08FC0476Cplus_fft_2_4_8_16snippet.dat -T 724 -D 5 -N 5 -Q P6PB08FC0476Cplus_fft_2_4_8_16snippet.Q1 -output_type gauss -anneal -annealfactor 1.1 -betamin 0.1 -regularize -omega 0 0 0.1 1.0e-6 -seed 1

The data file P6PB08FC0476Cplus_fft_2_4_8_16snippet.dat is a 724x5
flat ascii matrix of floating point values.

The model created by the above execution in default mode is then
used to classify a much larger time series using the evaluation
mode:

rdahmm -data P6PB08FC0476Cplus_fft_2_4_8_16.dat -T 123323 -D 5 -N 5 -output_type gauss -A P6PB08FC0476Cplus_fft_2_4_8_16snippet.A -B P6PB08FC0476Cplus_fft_2_4_8_16snippet.B -pi P6PB08FC0476Cplus_fft_2_4_8_16snippet.pi -minvalfile P6PB08FC0476Cplus_fft_2_4_8_16snippet.minval -maxvalfile P6PB08FC0476Cplus_fft_2_4_8_16snippet.maxval -rangefile P6PB08FC0476Cplus_fft_2_4_8_16snippet.range -Q P6PB08FC0476Cplus_fft_2_4_8_16.Q1 -eval -addstate
