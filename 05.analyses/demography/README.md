# Coalescent Demographic History 

## Stairway Plot
### _E. globulus_
Generated Site Frequency Spectrum using [easySFS](https://github.com/isaacovercast/easySFS). Manually generated population list file.

```bash
module load python/3.8
ESFS_DIR="/blue/soltis/kasey.pham/bin/easySFS"
VCF_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/04.filter_snps"
"$ESFS_DIR"/easySFS.py -i "$VCF_DIR"/all_fil_biallelic.vcf.gz -a -p esfs_poplist.txt --preview
# globulus segregating sites peaked (1330554) at subsampling to 32 haplotypes
# cordata segregating sites peaked (650658) at subsampling to 16 haplotypes
"$ESFS_DIR"/easySFS.py -i "$VCF_DIR"/all_fil_biallelic.vcf.gz -p esfs_poplist.txt -o esfs_outp -a -f --order glob_MR,cord_MR --proj 32,16 -v
```

Plotted 1D SFS and 2D SFS calculated by `EasySFS` in `R`.
```R
glob_sfs <- data.frame(num_shared = c(1:20), freq = c(150720.2660535347, 149704.1333331226, 138528.5552475863, 120754.5050957728, 106197.5229067253, 93712.27048029556, 83316.25583310747, 73960.8358039187, 65172.14064693877, 59158.41872999133, 54850.43365601913, 51473.79357392833, 49497.51698651399, 49987.6716989067, 54399.83560117622, 29120.30070290052, 0, 0, 0, 0))
cord_sfs <- data.frame(num_shared = c(1:10), freq = c(101696.2171998626, 122582.5098039221, 144140.6759545933, 142122.6738906096, 136740.7892672861, 131559.7115927073, 150248.7395940837, 94617.52012383903, 0, 0))

plot(glob_sfs, type = "l", main = "E. globulus folded SFS projected down to 32", xlab = "Number of Alleles Shared", ylab = "Frequency")
plot(cord_sfs, type = "l", main = "E. cordata folded SFS projected down to 16", xlab = "Number of Alleles Shared", ylab = "Frequency")
```

Created `Stairway Plot 2` run files.

```bash
module load java/20
PROG_DIR="/blue/soltis/kasey.pham/bin/stairway_plot_v2.1.1"

java -cp "$PROG_DIR"/stairway_plot_es Stairbuilder glob_fold.blueprint
java -cp "$PROG_DIR"/stairway_plot_es Stairbuilder cord_fold.blueprint
```

Ran `Stairway Plot 2` and graphed results.
```bash
# Ran in UFRC queue system; see stairwayplot2_glob.job for more details.
# Resources used: 310 Mb, 1 hr
module load java/1.8.0_31
PROG_DIR="/blue/soltis/kasey.pham/bin/stairway_plot_v2.1.1"

bash glob_fold.blueprint.sh
bash glob_fold.blueprint.plot.sh
```

### _E. cordata_
Used Site Frequency Spectrum generated during _E. globulus_ analysis to make blueprint file.

```bash
# Ran in UFRC queue system; see stairwayplot2_cord.job for more details.
# Resources used: 230 Mb, 15 min
module load java/1.8.0_31
PROG_DIR="/blue/soltis/kasey.pham/bin/stairway_plot_v2.1.1"

bash cord_fold.blueprint.sh
bash cord_fold.blueprint.plot.sh
```

## Demography model fitting
Used [`dadi`](https://dadi.readthedocs.io/en/latest/) to fit different demography models to observed variant distribution.

### 2D Models
Created folded SFS for observed data and bootstrapped variants in `dadi`. Followed the manual's [2D demography example](https://dadi.readthedocs.io/en/latest/examples/basic_workflow/basic_workflow_2d_demographics/). Masked central cells (16,8) and (15,8) as a precaution against [artifacts from mis-mapping](https://groups.google.com/g/dadi-user/c/thIHbLj5zHQ).
```python
# dadi demography model fitting
# step: site frequency spectrum construction and bootstrapping
# relies on code from dadi manual examples: https://dadi.readthedocs.io/en/latest/examples/fs_from_data/fs_from_data/

# Pickle is used to save variables as files for future use
import pickle
# MatPlotLib is a libary dadi uses for plotting frequency spectrum
import matplotlib.pyplot as plt
import dadi
# import os

# names for output files
wdir = "/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/demography/dadi/2D"
name_stem = "globMR_cordMR_ns32-16"

# Parse the VCF file to generate a data dictionary
datafile = "/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/04.filter_snps/all_fil_biallelic.vcf.gz"
dd = dadi.Misc.make_data_dict_vcf(datafile, "/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/demography/stairway_plot/esfs_poplist.txt")

# Extract the spectrum for MR E. globulus and MR E. cordata from that dictionary, with 
# E. globulus projected down to 32 of 40 and E. cordata projected down to 16 of 20.
# haplotype counts at which the number of segregating sites was maximized
pop_ids, ns = ["glob_MR", "cord_MR"], [32, 16]
fs = dadi.Spectrum.from_data_dict(dd, pop_ids, ns, polarized = False) # folded since ancestral alleles not known

# Mask central cells
fs.mask[16,8] = True
fs.mask[15,8] = True

# Make directory for saving data dictionaries
# if not os.path.exists("{WDIR}/data/data_dicts".format(WDIR = wdir)):
#  os.makedirs("{WDIR}/data/data_dicts".format(WDIR = wdir))

# Saved data dictionary to disk
pick = open("{WDIR}/data/data_dicts/{NAME}.bpkl".format(WDIR = wdir, NAME = name_stem), "wb")
pickle.dump(dd, pick, 2)
# Saved extracted spectrum
fs.to_file("{WDIR}/data/fs/{NAME}.fs".format(WDIR = wdir, NAME = name_stem))

# Generate 100 bootstrap datasets, by dividing the genome into 500kb chunks and
# resampling from those chunks.
Nboot, chunk_size = 100, 5e5
chunks = dadi.Misc.fragment_data_dict(dd, chunk_size)
boots = dadi.Misc.bootstraps_from_dd_chunks(chunks, Nboot, pop_ids, ns, polarized = False)
# Saved bootstraps
for i in range(len(boots)):
    boots[i].to_file("{WDIR}/data/fs/bootstraps/{NAME}.boot_{NUM}.fs".format(WDIR = wdir, NAME = name_stem, NUM = str(i)))

# Plotted site frequency spectrum
fig = plt.figure(982342)
fig.clear()
dadi.Plotting.plot_single_2d_sfs(fs)
fig.savefig("{WDIR}/data/fs/{NAME}.png".format(WDIR = wdir, NAME = name_stem))
plt.close(fig)
```

Ran dadi using the `dadi_pipeline` implementation. 

Run scheme:
1. 15 replicates, 3 max iterations, 3-fold perturbation
2. 30 replicates, 5 max iterations, 2-fold perturbation
3. 50 replicates, 10 max iterations, 2-fold perturbation
4. 100 replicates, 20 max iterations, 1-fold perturbation

Each model ran using ~80 Mb RAM and between 2 hrs and 2.5 days. See python dadi python files under each model directory for more detail.

Models run:
|Index | Model                                  | parameters                | python file           |
| ---- | -------------------------------------- | ------------------------- | --------------------- |
| 1    | Divergence, gradual size change        | nu1i, nu2i, nu1f, nu2f, T | dadi_schange.py       |
| 2    | Divergence, instant size change and gradual size change | nu1i, nu2i, nu1m, nu2m, nu1f, nu2f, T1, T2 | dadi_bottle_schange.py |
| 3    | Divergence, instant size change, gradual size change | nu1i, nu2i, nu1m, nu2m, nu1f, nu2f, T1, T2, T3 | dadi_bottle_schange_thr_epoch.py |
| 4    | Divergence, asymmetric secondary contact and gradual size change | nu1i, nu2i, nu1f, nu2f, m12, m21, T1, T2 | dadi_sec_contact_schange.py |
| 5    | Divergence, instant size change and asymmetric secondary contact and gradual size change | nu1i, nu2i, nu1m, nu2m, nu1f, nu2f, m12, m21, T1, T2 | dadi_sec_contact_bottle_schange.py |
| 6    | Divergence, instant size change and asymmetric secondary contact, gradual size change | nu1i, nu2i, nu1m, nu2m, nu1f, nu2f, m12, m21, T1, T2 | dadi_sec_contact_bottle_schange_thr_epoch.py |

Extracted best scores from each run in `R`. Example code given below:
```R
name_pref <- "globMR_cordMR_run"
i_list <- c(1,2,3,4,5)
name_suff <- ".schange.optimized.txt"

tab_list <- lapply(i_list, function(i) read.table(paste(name_pref, i, name_suff, sep = ""), header = TRUE, sep = "\t")[c(96:195),])
best_scores <- lapply(tab_list, function(tab) tab[which(tab$AIC == min(tab$AIC)),])
```

Summary of results of best run for each model: 
| Index | ID                                   | Log Likelihood | AIC      |
| ----- | ------------------------------------ | -------------- | -------- |
| 1     | schange                              | -155089.2      | 310188.3 |
| 2     | bottle_schange                       | -168961.8      | 337939.6 |
| 3     | bottle_schange_thr_epoch             | -227180.9      | 454379.8 |
| 4     | sec_contact_schange                  | -140382.2      | 280780.4 |
| 5     | sec_contact_bottle_schange           | -225817.0      | 451654.0 |
| 6     | sec_contact_bottle_schange_thr_epoch | -158301.2      | 316624.3 |
More details and parameter values of each best run can be found in `dadi_2D_results_summary.xlxs`.

Plot SFS fit and likelihood ratio test for selecting the best model:
```python
import dadi
import numpy as np
import nlopt
import matplotlib.pyplot as plt
import pylab
import sys

# location of model functions
sys.path.append("/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/demography/dadi/2D")
import all_snps_models

wdir = "/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/demography/dadi/2D"
name_stem = "globMR_cordMR_ns32-16"

# import SFS
fs = dadi.Spectrum.from_file("{WDIR}/data/fs/{NAME}.fs".format(WDIR = wdir, NAME = name_stem))
ns = fs.sample_sizes

# establish grid size
pts = [max(ns)+20, max(ns)+30, max(ns)+40]

# import bootstraps
boots = []
for i in range(100):
    boots.append(dadi.Spectrum.from_file("{WDIR}/data/fs/bootstraps/{NAME}.boot_{IND}.fs".format(WDIR = wdir, NAME = name_stem, IND = i)))

# record results from each model's parameterization
# schange
smodel_bfps = [3.5617, 27.1096, 0.0256, 0.011, 0.0339]
# sec_contact_bottle_schange
cmodel_bfps = [4.3525, 26.3731, 0.0356, 0.0248, 10.7472, 1.0565, 0.3154, 0.0899]

nested_ind = [4,5,7]

# extrapolate function and make best fit model
smodel_func_ex = dadi.Numerics.make_extrap_log_func(all_snps_models.schange)
cmodel_func_ex = dadi.Numerics.make_extrap_log_func(all_snps_models.sec_contact_schange)

smodel = smodel_func_ex(smodel_bfps, ns, pts)
cmodel = cmodel_func_ex(cmodel_bfps, ns, pts)

## PLOT FIT
# simple model
fig = plt.figure(1, figsize=(10,6))
fig.clear()
dadi.Plotting.plot_2d_comp_multinom(smodel, fs, vmin=1, resid_range=3, pop_ids =('glob_MR','cord'), show=False)
pylab.savefig('schange_fit.png', dpi=250)

# complex model
fig = plt.figure(2, figsize=(10,6))
fig.clear()
dadi.Plotting.plot_2d_comp_multinom(cmodel, fs, vmin=1, resid_range=3, pop_ids =('glob_MR','cord'), show=False)
pylab.savefig('sec_contact_schange_fit.png', dpi=250)

## LIKELIHOOD RATIO TEST
adj = dadi.Godambe.LRT_adjust(func_ex=cmodel_func_ex, grid_pts=pts, all_boot=boots, p0=cmodel_bfps, data=fs, nested_indices=nested_ind, multinom=True)
# 0.0006056
D = adj*2*(-140382.2 - (-155089.2)) # 17.813
p_val = dadi.Godambe.sum_chi2_ppf(D, weights=[0.125, 0.375, 0.375, 0.125])
# P-value = 0.0001200, so Model 4 is significantly a better fit than Model 1.

## PARAMETER UNCERTAINTY ESTIMATION
uncert = dadi.Godambe.GIM_uncert(func_ex=cmodel_func_ex, grid_pts=pts, all_boot=boots, p0=cmodel_bfps, data=fs, log=True, multinom=True)
```

Using the Godambe Likelihood Adjustment, the P-value for Model 4 over Model 1 was 0.0001200, so Model 4 was significantly better of a fit.

Parameter uncertainty estimates for Model 4: Secondary Contact with gradual size change:
| Parameter | Point Estimate | Uncertainty |
| --------- | -------------- | ----------- |
| nu1i      | 4.3525         | 1.33653495  |
| nu2i      | 26.3731        | 1.57506129  |
| nu1f      | 0.0356         | 3.65428061  |
| nu2f      | 0.0248         | 3.22294423  |
| m12       | 10.7472        | 3.22662747  |
| m21       | 1.0565         | 4.17443425  |
| T1        | 0.3154         | 0.69228957  |
| T2        | 0.0899         | 3.03299481  |
| theta     | 404517.5       | 0.11088168  |

### 1D Models
Repeated SFS generation step for 1D SFS.
```python
# dadi demography model fitting
# step: site frequency spectrum construction and bootstrapping
# relies on code from dadi manual examples: https://dadi.readthedocs.io/en/latest/examples/fs_from_data/fs_from_data/

# Pickle is used to save variables as files for future use
import pickle
# MatPlotLib is a libary dadi uses for plotting frequency spectrum
import matplotlib.pyplot as plt
import dadi

# Parse the VCF file to generate a data dictionary
datafile = "/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/04.filter_snps/all_fil_biallelic.vcf.gz"
dd = dadi.Misc.make_data_dict_vcf(datafile, "/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/demography/stairway_plot/esfs_poplist.txt")

## GLOBULUS
# names for output files
wdir = "/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/demography/dadi/1D/glob"
name_stem_g = "globMR_ns32"

# Extract 1D spectrum for MR E. globulus from that dictionary, with 
# E. globulus projected down to 32 of 40.
# haplotype counts at which the number of segregating sites was maximized
pop_id_g = ["glob_MR"]
ns_g = [32]
fs_g = dadi.Spectrum.from_data_dict(dd, pop_id_g, ns_g, polarized = False) # folded

# Mask central cells
fs_g.mask[16] = True
fs_g.mask[15] = True

# Saved data dictionary to disk
pick = open("{WDIR}/data/data_dicts/{NAME}.bpkl".format(WDIR = wdir, NAME = name_stem_g), "wb")
pickle.dump(dd, pick, 2)
# Saved extracted spectrum
fs_g.to_file("{WDIR}/data/fs/{NAME}.fs".format(WDIR = wdir, NAME = name_stem_g))

# Generate 100 bootstrap datasets, by dividing the genome into 500kb chunks and
# resampling from those chunks.
Nboot, chunk_size = 100, 5e5
chunks = dadi.Misc.fragment_data_dict(dd, chunk_size)
boots = dadi.Misc.bootstraps_from_dd_chunks(chunks, Nboot, pop_id_g, ns_g, polarized = False)
# Saved bootstraps
for i in range(len(boots)):
    boots[i].to_file("{WDIR}/data/fs/bootstraps/{NAME}.boot_{NUM}.fs".format(WDIR = wdir, NAME = name_stem_g, NUM = str(i)))

# Plotted site frequency spectrum
fig = plt.figure(1)
fig.clear()
dadi.Plotting.plot_1d_fs(fs_g)
fig.savefig("{WDIR}/data/fs/{NAME}.png".format(WDIR = wdir, NAME = name_stem_g))
plt.close(fig)

## CORDATA
wdir = "/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/demography/dadi/1D/cord"
name_stem_c = "cordMR_ns16"

# Extract 1D spectrum for MR E. cordata from that dictionary, with 
# E. globulus projected down to 16 of 20.
# haplotype counts at which the number of segregating sites was maximized
pop_id_c = ["cord_MR"]
ns_c = [16]
fs_c = dadi.Spectrum.from_data_dict(dd, pop_id_c, ns_c, polarized = False) # folded

# Mask central cells
fs_c.mask[8] = True
fs_c.mask[7] = True

# Saved data dictionary to disk
pick = open("{WDIR}/data/data_dicts/{NAME}.bpkl".format(WDIR = wdir, NAME = name_stem_c), "wb")
pickle.dump(dd, pick, 2)
# Saved extracted spectrum
fs_c.to_file("{WDIR}/data/fs/{NAME}.fs".format(WDIR = wdir, NAME = name_stem_c))

# Generate 100 bootstrap datasets, by dividing the genome into 500kb chunks and
# resampling from those chunks.
Nboot, chunk_size = 100, 5e5
chunks = dadi.Misc.fragment_data_dict(dd, chunk_size)
boots = dadi.Misc.bootstraps_from_dd_chunks(chunks, Nboot, pop_id_c, ns_c, polarized = False)
# Saved bootstraps
for i in range(len(boots)):
    boots[i].to_file("{WDIR}/data/fs/bootstraps/{NAME}.boot_{NUM}.fs".format(WDIR = wdir, NAME = name_stem_c, NUM = str(i)))

# Plotted site frequency spectrum
fig = plt.figure(2)
fig.clear()
dadi.Plotting.plot_1d_fs(fs_c)
fig.savefig("{WDIR}/data/fs/{NAME}.png".format(WDIR = wdir, NAME = name_stem_c))
plt.close(fig)
```

Tested the same models for _E. globulus_ and _E. cordata_ from `dadi` 1D model set using Daniel Portik's `dadi_pipeline` scripts for multi-round parameter optimization.

Run scheme:
1. 15 replicates, 3 max iterations, 3-fold perturbation
2. 30 replicates, 5 max iterations, 2-fold perturbation
3. 50 replicates, 10 max iterations, 2-fold perturbation
4. 100 replicates, 20 max iterations, 1-fold perturbation

Models ran using ~80 Mb RAM and a maximum of 10 min.

Models tested:
| Index | Description                                | Parameters         | python file                  |
| ----- | ------------------------------------------ | ------------------ | ---------------------------- |
| 2     | Exponential size change                    | nu, T              | dadi_growth.job              |
| 3     | Neutral evolution, then instantaneous size change | nu, T       | dadi_two_epoch.job           |
| 4     | Instantaneous size change, then exponential size change | nuB, nuF, T | dadi_bottlegrowth.job  |
| 5     | Exponential size change, neutral evolution, exponential size change | nuB, nuF, TB, TF | dadi_three_epoch.job |

Extracted best scores from each run in `R`. Example code given below:
```R
name_pref <- "globMR_run"
i_list <- c(1,2,3,4,5)
name_suff <- ".growth.optimized.txt"

tab_list <- lapply(i_list, function(i) read.table(paste(name_pref, i, name_suff, sep = ""), header = TRUE, sep = "\t")[c(96:195),])
best_scores <- lapply(tab_list, function(tab) tab[which(tab$AIC == min(tab$AIC)),])
```

The optimal model for both _E. globulus_ and _E. cordata_ was the three epoch model.

| Parameter | _E. globulus_ | _E. cordata_ |
| --------- | ------------- | ------------ |
| theta     | 94640.24      | 91424.01     |
| nuB       | 23.9132       | 26.5665      |
| nuF       | 0.0462        | 0.0101       |
| TB        | 5.1996        | 23.9653      |
| TF        | 0.0140        | 0.0141       |

Calculated real-life parameters for optimal models, using L = 14,226,796.71 bp, generation time = 50 years, and mutation rate = 4.93e-09.

| Parameter | _E. globulus_ | _E. cordata_  |
| --------- | ------------- | ------------- |
| Nref      | 6746.71       | 6517.43       |
| nuB       | 161335.33     | 173145.24     |
| nuF       | 311.70        | 65.83         |
| TB        | 3508017.31    | 15619211.13   |
| TF        | 9445.39       | 9189.57       |

More details on calculations and values from other models can be found in the files `dadi_glob1D_results_summary.xlxs` and `dadi_cord1D_results_summary.xlsx`.

Plotted best fit model to evaluate fit to actual dataset and performed Godambe likelihood ratio testing on more complicated versus more simple models in both taxa.

For _E. globulus_:
```python
## ----------- ##
## E. GLOBULUS ##
## ----------- ##

import dadi
import numpy as np
import nlopt
import matplotlib.pyplot as plt
import pylab
import sys

wdir = "/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/demography/dadi/1D/glob"
name_stem = "globMR_ns32"

# import SFS
fs = dadi.Spectrum.from_file("{WDIR}/data/fs/{NAME}.fs".format(WDIR = wdir, NAME = name_stem))
ns = fs.sample_sizes

# establish grid size
pts = [max(ns)+20, max(ns)+30, max(ns)+40]

# import bootstraps
boots = []
for i in range(100):
    boots.append(dadi.Spectrum.from_file("{WDIR}/data/fs/{NAME}.fs".format(WDIR = wdir, NAME = name_stem)))

# record results from each model's parameterization
# growth (AIC = 8467.06)
g_bfps = [0.0224, 0.01]
# bottlegrowth (AIC = 8440.84)
bg_bfps = [0.8337, 0.023, 0.01]
# three epoch (AIC = 2052.18)
te_bfps = [23.747, 0.0659, 5.2131, 0.0201]

g_bg_nested = [1] # nuF = nuB, interior
bg_te_nested = [3] # TF = 0, boundary

# extrapolate function and make best fit model
g_func_ex = dadi.Numerics.make_extrap_log_func(dadi.Demographics1D.growth)
bg_func_ex = dadi.Numerics.make_extrap_log_func(dadi.Demographics1D.bottlegrowth)
te_func_ex = dadi.Numerics.make_extrap_log_func(dadi.Demographics1D.three_epoch)

g_model = g_func_ex(g_bfps, ns, pts)
bg_model = bg_func_ex(bg_bfps, ns, pts)
te_model = te_func_ex(te_bfps, ns, pts)

## PLOT FIT
# growth
fig = plt.figure(1, figsize=(10,6))
fig.clear()
dadi.Plotting.plot_1d_comp_multinom(g_model, fs)
pylab.savefig('./02.growth/glob_growth_fit.png', dpi=250)

# bottlegrowth
fig = plt.figure(2, figsize=(10,6))
fig.clear()
dadi.Plotting.plot_1d_comp_multinom(bg_model, fs,)
pylab.savefig('./04.bottlegrowth/glob_bottlegrowth_fit.png', dpi=250)

# three epoch
fig = plt.figure(3, figsize=(10,6))
fig.clear()
dadi.Plotting.plot_1d_comp_multinom(te_model, fs)
pylab.savefig('./05.three_epoch/glob_three_epoch_fit.png', dpi=250)

## LIKELIHOOD RATIO TEST
# growth vs. bottlegrowth
adj = dadi.Godambe.LRT_adjust(func_ex=bg_func_ex, grid_pts=pts, all_boot=boots, p0=bg_bfps, data=fs, nested_indices=g_bg_nested, multinom=True)
D = adj*2*(4231.53 - 4217.42)
p_val = dadi.Godambe.sum_chi2_ppf(D, weights=(0,1)) # Pval = 0.0; significantly different

# bottlegrowth vs. three epoch
# This isn't right, they're not actually nested models.
# adj = dadi.Godambe.LRT_adjust(func_ex=te_func_ex, grid_pts=pts, all_boot=boots, p0=te_bfps, data=fs, nested_indices=bg_te_nested, multinom=True)
# D = adj*2*(4231.53 - 1022.09)
# p_val = dadi.Godambe.sum_chi2_ppf(D, weights=(0.5,0.5)) # Pval = 0.0; significantly different

## PARAMETER UNCERTAINTY
uncert = dadi.Godambe.GIM_uncert(func_ex=te_func_ex, grid_pts=pts, all_boot=boots, p0=te_bfps, data=fs, log = True, multinom = True)
```
Conclusion: the three_epoch model is the best for _E. globulus_. Uncertainty estimates from Godambe Information Matrix:
| Parameter | Point Estimate | Uncertainty |
| --------- | -------------- | ----------- |
| nuB       | 23.747         | 0.00385118  |
| nuF       | 0.0659         | 0.00534308  |
| TB        | 5.2131         | 0.00213894  |
| TF        | 0.0201         | 0.00506429  |
| theta     | 58556.6        | 0.00198105  |

Did the same for _E. cordata_:
```python
## ----------- ##
## E. CORDATA  ##
## ----------- ##

import dadi
import numpy as np
import nlopt
import matplotlib.pyplot as plt
import pylab
import sys

wdir = "/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/demography/dadi/1D/cord"
name_stem = "cordMR_ns16"

# import SFS
fs = dadi.Spectrum.from_file("{WDIR}/data/fs/{NAME}.fs".format(WDIR = wdir, NAME = name_stem))
ns = fs.sample_sizes

# establish grid size
pts = [max(ns)+20, max(ns)+30, max(ns)+40]

# import bootstraps
boots = []
for i in range(100):
    boots.append(dadi.Spectrum.from_file("{WDIR}/data/fs/{NAME}.fs".format(WDIR = wdir, NAME = name_stem)))

# record results from each model's parameterization
# two epoch (AIC = 8343.58)
twe_bfps = [0.01, 0.01]
# three epoch (AIC = 6177.92)
the_bfps = [27.6535, 0.014, 24.8504, 0.0213]

twe_the_nested = [1,3] #nuF = nuB, TF = 0, 1 interior and 1 boundary 

# extrapolate function and make best fit model
twe_func_ex = dadi.Numerics.make_extrap_log_func(dadi.Demographics1D.two_epoch)
the_func_ex = dadi.Numerics.make_extrap_log_func(dadi.Demographics1D.three_epoch)

twe_model = twe_func_ex(twe_bfps, ns, pts)
the_model = the_func_ex(the_bfps, ns, pts)

## PLOT FIT
# growth
fig = plt.figure(1, figsize=(10,6))
fig.clear()
dadi.Plotting.plot_1d_comp_multinom(twe_model, fs)
pylab.savefig('./03.two_epoch/cord_two_epoch_fit.png', dpi=250)

# three epoch
fig = plt.figure(2, figsize=(10,6))
fig.clear()
dadi.Plotting.plot_1d_comp_multinom(the_model, fs)
pylab.savefig('./05.three_epoch/cord_three_epoch_fit.png', dpi=250)

## LIKELIHOOD RATIO TEST
# two epoch vs. three epoch
adj = dadi.Godambe.LRT_adjust(func_ex=the_func_ex, grid_pts=pts, all_boot=boots, p0=the_bfps, data=fs, nested_indices=twe_the_nested, multinom=True) # why is this negative?? D can't be negative
D = adj*2*(4169.79 - 3084.96)
p_val = dadi.Godambe.sum_chi2_ppf(D, weights=(0, 0.5, 0.5)) # Pval = 1.0; not significantly different

# PARAMETER UNCERTAINTY
uncert = dadi.Godambe.GIM_uncert(func_ex=twe_func_ex, grid_pts=pts, all_boot=boots, p0=twe_bfps, data=fs, log = True, multinom = True)
# uncert was all NAs
```

The Godambe adjustment factor was negative, which yields a not-possible D value. I think the fit is just very poor on all of these models and maybe there isn't enough information in the bootstraps as well.