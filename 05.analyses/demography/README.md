# Coalescent Demographic History 

## Stairway Plot
### _E. globulus_
Generated Site Frequency Spectrum using [easySFS](https://github.com/isaacovercast/easySFS). Manually generated population list file.

```bash
module load python/3.8
ESFS_DIR="/blue/soltis/kasey.pham/bin/easySFS"
VCF_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/04.filter_snps"
"$ESFS_DIR"/easySFS.py -i "$VCF_DIR"/all_fil.vcf.gz -a -p esfs_poplist.txt --preview
# globulus segregating sites peaked (876054) at subsampling to 32 haplotypes
# cordata segregating sites peaked (690386) at subsampling to 16 haplotypes
"$ESFS_DIR"/easySFS.py -i "$VCF_DIR"/all_fil.vcf.gz -p esfs_poplist.txt -o esfs_outp -a -f --order glob_MR,cord_MR --proj 32,16 -v
```

Created `Stairway Plot 2` input file and ran program.

```bash
# Ran in UFRC queue system; see stairwayplot2_glob.job for more details.
# Resources used:

module load java/20
PROG_DIR="/blue/soltis/kasey.pham/bin/stairway_plot_v2.1.1"

java -cp "$PROG_DIR"/stairway_plot_es Stairbuilder glob_fold.blueprint
bash glob_fold.blueprint.sh
```

### _E. cordata_
Used Site Frequency Spectrum generated during _E. globulus_ analysis to make blueprint file.

```bash
# Ran in UFRC queue system; see stairwayplot2_cord.job for more details.
# Resources used:
module load java/20
PROG_DIR="/blue/soltis/kasey.pham/bin/stairway_plot_v2.1.1"

java -cp "$PROG_DIR"/stairway_plot_es Stairbuilder cord_fold.blueprint
bash cord_fold.blueprint.sh
```

## Demography model fitting
Used [`dadi`](https://dadi.readthedocs.io/en/latest/) to fit different demography models to observed variant distribution.

Created folded SFS for observed data and bootstrapped variants in `dadi`. Followed the manual's [2D demography example](https://dadi.readthedocs.io/en/latest/examples/basic_workflow/basic_workflow_2d_demographics/). Masked central cells (8,16), (8,15), (7,16) as a precaution against [artifacts from mis-mapping](https://groups.google.com/g/dadi-user/c/thIHbLj5zHQ).
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
datafile = "/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/04.filter_snps/all_fil_biallelic.vcf"
dd = dadi.Misc.make_data_dict_vcf(datafile, "/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/demography/stairway_plot/esfs_poplist.txt")

# Extract the spectrum for MR E. globulus and MR E. cordata from that dictionary, with 
# E. globulus projected down to 32 of 40 and E. cordata projected down to 16 of 20.
# haplotype counts at which the number of segregating sites was maximized
pop_ids, ns = ["glob_MR", "cord_MR"], [32, 16]
fs = dadi.Spectrum.from_data_dict(dd, pop_ids, ns, polarized = False) # folded since ancestral alleles not known

# Mask central cells
fs.mask[8,16] = True
fs.mask[8,15] = True
fs.mask[7,16] = True

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
4. 80 replicates, 20 max iterations, 1-fold perturbation

See python dadi python files under each model directory for more detail.

Models run:
| Model                                  | parameters                | python file           |
| -------------------------------------- | ------------------------- | --------------------- |
| Divergence, gradual size change        | nu1i, nu2i, nu1f, nu2f, T | dadi_schange.py       |
| Divergence, instant size change and gradual size change | nu1i, nu2i, nu1m, nu2m, nu1f, nu2f, T1, T2 | dadi_bottle_schange.py |
| Divergence, asymmetric secondary contact and gradual size change | nu1i, nu2i, nu1f, nu2f, m12, m21, T1, T2 |
| Divergence, instant size change and asymmetric secondary contact and gradual size change | nu1i, nu2i, nu1m, nu2m, nu1f, nu2f, m12, m21, T1, T2 |

Plot SFS fit and likelihood ratio test for selecting the best model:
```python
import dadi
import numpy as np
import nlopt
import matplotlib.pyplot as plt
import pylab
import sys

sys.path.append("/blue/soltis/kasey.pham/bin/dadi_pipeline")
sys.path.append("/blue/soltis/kasey.pham/bin/dadi_pipeline/Two_Population_Pipeline")
# location of model functions
sys.path.append("/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/demography/dadi/all_snps")
import all_snps_models

wdir = "/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/demography/dadi/all_snps"
name_stem = "globMR_cordMR_ns32-16"

# import SFS
fs = dadi.Spectrum.from_file("{WDIR}/data/fs/{NAME}.fs".format(WDIR = wdir, NAME = name_stem))
ns = fs.sample_sizes

# establish grid size
pts = [max(ns)+20, max(ns)+30, max(ns)+40]

# import bootstraps
boots = []
for i in range(100):
    boots.append(dadi.Spectrum.from_file("{WDIR}/data/fs/{NAME}.fs".format(WDIR = wdir, NAME = name_stem)))

# record results from sec_contact_schange parameterization
smodel_bfps = [10.8079, 0.2501, 0.01, 0.0164, 0.0185]
cmodel_bfps = [0.4248, 17.0388, 0.0844, 0.0352, 8.1421, 1.9196, 4.9025, 0.1318]
nested_ind = [4,5,6]

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
fig = plt.figure(1, figsize=(10,6))
fig.clear()
dadi.Plotting.plot_2d_comp_multinom(cmodel, fs, vmin=1, resid_range=3, pop_ids =('glob_MR','cord'), show=False)
pylab.savefig('sec_contact_schange_fit.png', dpi=250)

## LIKELIHOOD RATIO TEST
# adj = dadi.Godambe.LRT_adjust(func_ex=func_ex, grid_pts=pts, all_boot=boots, p0=cmodel_bfps, data=fs, nested_indices=nested_ind, multinom=True)
```