'''
Usage: python dadi_bottle_schange.py
For running param optimization for demographic model divergence with an instantaneous
size change followed by gradual exponential size change.

This is a modified version of the 'dadi_Run_Optimizations.py' script in which
we run optimizations for 2D comparisons for a large set of models that have been
made available as part of published works. These models are stored in the
Models_2D.py script, and will be called directly here. The user can delete or
comment out models to analyze a subset of the models available. 

This script must be in the same working directory as Optimize_Functions.py, which
contains all the functions necessary, as well as the  Models_2D.py script, which
has all the model definitions.


General workflow:
 The optimization routine runs a user-defined number of rounds, each with a user-defined
 or predefined number of replicates. The starting parameters are initially random, but after
 each round is complete the parameters of the best scoring replicate from that round are
 used to generate perturbed starting parameters for the replicates of the subsequent round.
 The arguments controlling steps of the optimization algorithm (maxiter) and perturbation
 of starting parameters (fold) can be supplied by the user for more control across rounds.
 The user can also supply their own set of initial parameters, or set custom bounds on the
 parameters (upper_bound and lower_bound) to meet specific model needs. This flexibility
 should allow these scripts to be generally useful for model-fitting with any data set.

 
Outputs:
 For each model run, there will be a log file showing the optimization steps per replicate
 and a summary file that has all the important information. Here is an example of the output
 from a summary file, which will be in tab-delimited format:
 
 Model	Replicate	log-likelihood	AIC	chi-squared	theta	optimized_params(nu1, nu2, m, T)
 sym_mig	Round_1_Replicate_1	-1684.99	3377.98	14628.4	383.04	0.2356,0.5311,0.8302,0.182
 sym_mig	Round_1_Replicate_2	-2255.47	4518.94	68948.93	478.71	0.3972,0.2322,2.6093,0.611
 sym_mig	Round_1_Replicate_3	-2837.96	5683.92	231032.51	718.25	0.1078,0.3932,4.2544,2.9936
 sym_mig	Round_1_Replicate_4	-4262.29	8532.58	8907386.55	288.05	0.3689,0.8892,3.0951,2.8496
 sym_mig	Round_1_Replicate_5	-4474.86	8957.72	13029301.84	188.94	2.9248,1.9986,0.2484,0.3688


Notes/Caveats:
 The likelihood and AIC returned represent the true likelihood only if the SNPs are
 unlinked across loci. For ddRADseq data where a single SNP is selected per locus, this
 is true, but if SNPs are linked across loci then the likelihood is actually a composite
 likelihood and using something like AIC is no longer appropriate for model comparisons.
 See the discussion group for more information on this subject. 

Citations:
 If you use these scripts or the main diversification models for your work, please
 cite the following publication:
    Portik, D.M., Leache, A.D., Rivera, D., Blackburn, D.C., Rodel, M.-O.,
    Barej, M.F., Hirschfeld, M., Burger, M., and M.K. Fujita. 2017.
    Evaluating mechanisms of diversification in a Guineo-Congolian forest
    frog using demographic model selection. Molecular Ecology 26: 5245-5263.
    doi: 10.1111/mec.14266
 
 If you use the additional diversification models or the island models set please cite 
 the following publication:
    Charles, K.C., Bell, R.C., Blackburn, D.C., Burger, M., Fujita, M.K.,
    Gvozdik, V., Jongsma, G.F.M., Leache, A.D., and D.M. Portik. Sky, sea,
    and forest islands: diversification in the African leaf-folding frog
    Afrixalus paradorsalis (Order: Anura, Family: Hyperoliidae).
    Journal of Biogeography 45: 1781-1794. 
    doi: 10.1111/jbi.13365
        
 If you are interesting in contributing your models to this workflow, please email me!

-------------------------
Written for Python 2.7 and 3.7
Python modules required:
-Numpy
-Scipy
-dadi
-------------------------

Daniel Portik
daniel.portik@gmail.com
https://github.com/dportik
Updated September 2019
'''

import sys
import os
import numpy
import dadi
import pylab
from datetime import datetime
# add locations of dadi_pipeline functions
sys.path.append("/blue/soltis/kasey.pham/bin/dadi_pipeline")
sys.path.append("/blue/soltis/kasey.pham/bin/dadi_pipeline/Two_Population_Pipeline")
# location of model functions
sys.path.append("/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/demography/dadi/2D")
import Optimize_Functions
import all_snps_models


wdir = "/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/demography/dadi/2D"
name_stem = "globMR_cordMR_ns32-16"

#===========================================================================
# Import data to create joint-site frequency spectrum
#===========================================================================

#**************
snps = "/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/04.filter_snps/all_fil_biallelic.vcf"

#Create python dictionary from snps file
# dd = dadi.Misc.make_data_dict(snps)

#**************
#pop_ids is a list which should match the populations headers of your SNPs file columns
pop_ids=["glob_MR", "cord_MR"]

#**************
#projection sizes, in ALLELES not individuals
proj = [32, 16]

#Convert this dictionary into folded AFS object
#[polarized = False] creates folded spectrum object
# fs = dadi.Spectrum.from_data_dict(dd, pop_ids=pop_ids, projections = proj, polarized = False)

# Read site frequency spectrum
fs = dadi.Spectrum.from_file("{WDIR}/data/fs/{NAME}.fs".format(WDIR = wdir, NAME = name_stem))
ns = fs.sample_sizes

#print some useful information about the afs or jsfs
print("\n\n============================================================================")
print("\nData for site frequency spectrum:\n")
print("Projection: {}".format(proj))
print("Sample sizes: {}".format(fs.sample_sizes))
print("Sum of SFS: {}".format(numpy.around(fs.S(), 2)))
print("\n============================================================================\n")

#================================================================================
# Calling external 2D models from the Models_2D.py script
#================================================================================
'''
 We will use a function from the Optimize_Functions.py script for our optimization routines:
 
 Optimize_Routine(fs, pts, outfile, model_name, func, rounds, param_number, fs_folded=True, 
                        reps=None, maxiters=None, folds=None, in_params=None, 
                        in_upper=None, in_lower=None, param_labels=None, optimizer="log_fmin")
 
   Mandatory Arguments =
    fs:  spectrum object name
    pts: grid size for extrapolation, list of three values
    outfile:  prefix for output naming
    model_name: a label to help label the output files; ex. "no_mig"
    func: access the model function from within 'moments_Run_Optimizations.py' or from a separate python model script, 
          ex. after importing Models_2D, calling Models_2D.no_mig
    rounds: number of optimization rounds to perform
    param_number: number of parameters in the model selected (can count in params line for the model)
    fs_folded: A Boolean value (True or False) indicating whether the empirical fs is folded (True) or not (False).

   Optional Arguments =
     reps: a list of integers controlling the number of replicates in each of the optimization rounds
     maxiters: a list of integers controlling the maxiter argument in each of the optimization rounds
     folds: a list of integers controlling the fold argument when perturbing input parameter values
     in_params: a list of parameter values 
     in_upper: a list of upper bound values
     in_lower: a list of lower bound values
     param_labels: list of labels for parameters that will be written to the output file to keep track of their order
     optimizer: a string, to select the optimizer. Choices include: "log" (BFGS method), "log_lbfgsb" (L-BFGS-B method), 
                "log_fmin" (Nelder-Mead method), and "log_powell" (Powell's method).
'''


#create a prefix based on the population names to label the output files
#ex. Pop1_Pop2
#prefix = "_".join(pop_ids)
prefix = "globMR_cordMR"

#**************
# Define the grid points based on the sample size.
# Based on suggested sizes in dadi manual for a small sample size
pts = [max(ns)+20, max(ns)+30, max(ns)+40]

#**************
#Set the number of rounds here
rounds = 4

#define the lists for optional arguments
#you can change these to alter the settings of the optimization routine
reps = [15,30,50,100]
maxiters = [3,5,10,20]
folds = [3,2,2,1]

#**************
#Indicate whether your frequency spectrum object is folded (True) or unfolded (False)
fs_folded = True


'''
Diversification Model Set

This first set of models come from the following publication:

    Portik, D.M., Leache, A.D., Rivera, D., Blackburn, D.C., Rodel, M.-O.,
    Barej, M.F., Hirschfeld, M., Burger, M., and M.K.Fujita. 2017.
    Evaluating mechanisms of diversification in a Guineo-Congolian forest
    frog using demographic model selection. Molecular Ecology 26: 5245-5263.
    doi: 10.1111/mec.14266

'''

# Run optimization for model
# 02. Divergence with instantaneous size change and gradual exponential size change
for i in range(1,6):
   Optimize_Functions.Optimize_Routine(fs, pts, "{PREFIX}_run{NUM}".format(PREFIX = prefix, NUM = i), "bottle_schange", all_snps_models.bottle_schange, rounds, 8, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds, param_labels = "nu1i, nu2i, nu1m, nu2m, nu1f, nu2f, T1, T2")
