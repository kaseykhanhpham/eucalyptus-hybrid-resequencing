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

### LD-pruned SNPs
First performed linkage pruning on the set of _biallelic SNPs only_ in `PLINK` (as only biallelic SNPs can be used to estimate site frequency spectrum).
```bash
# Performed on UFRC queue system; see link_prune_biall.job for more details.
# Resources used: 
module load plink/1.90b3.39 

INDIR="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/04.filter_snps"

# prune linked biallelic SNPs
plink --vcf "$INDIR"/all_fil_biallelic.vcf.gz --double-id --allow-extra-chr --set-missing-var-ids @:# --indep-pairwise 50 10 0.2 --vcf-half-call m --out all_fil_biallelic
# extract VCF file from pruned SNPs
plink --vcf "$INDIR"/all_fil_biallelic.vcf.gz --double-id --allow-extra-chr --set-missing-var-ids @:# --extract all_fil_biallelic.prune.in --vcf-half-call m --recode vcf-iid --out "$INDIR"/all_fil_biallelic.prune
```

Created folded SFS for observed data and bootstrapped variants in `dadi`. Followed the manual's [2D demography example](https://dadi.readthedocs.io/en/latest/examples/basic_workflow/basic_workflow_2d_demographics/).
```bash
module load dadi
python dadi_prep.py
```

Ran parameter optimization for all 2D demographic models as described in Portik et al. 2017 and the associated pipeline, [`dadi_pipeline`](https://github.com/dportik/dadi_pipeline/tree/master). See individual python files and `dadi_pipeline`'s documentation for more details. In short, ran suggested 4 round optimization, using the previous round's best replicate (judged by likelihood) as a seed for heuristic search starting points with decreasing amounts of perturbation.

I used the default pipeline suggested by Portik et al. -- replicates of [10, 20, 30, 40], max iterations of [3, 5, 10, 15], and perturbations of [3, 2, 2, 1] for each respective round. I ran this optimization process 5 times for each of the models tested to check that the same optima were being recovered consistently. For the best model, I ran the optimization process 10 total times.

([Refer to `dadi_pipeline` for the full set of models coded by Portik et al.](https://github.com/dportik/dadi_pipeline/blob/master/Two_Population_Pipeline/Models_2D.pdf). I tested all the divergence models and none of the vicariance/island models.

Models used:
| Model                                          | Parameters         | python run file          |
| ---------------------------------------------- | ------------------ | ------------------------ |
| Divergence, no migration                       | nu1, nu2, T        | dadi_no_mig.py           |
| Divergence with continuous symmetric migration  | nu1, nu2, m, T    | dadi_sym_mig.py          |
| Divergence with continuous asymmetric migration | nu1, nu2, m12, m21, T | dadi_asym_mig.py     |
| Divergence with no migration, size change      | nu1a, nu2a, nu1b, nu2b, T1, T2 | dadi_no_mig_size.py |
| Divergence with continuous symmetric migration, size change | nu1a, nu2a, nu1b, nu2b, m, T1, T2 | dadi_sym_mig_size.py |
| Divergence with continuous asymmetric migration, size change | nu1a, nu2a, nu1b, nu2b, m12, m21, T1, T2 | dadi_asym_mig_size.py |
| Divergence in isolation, continuous symmetrical secondary contact | nu1, nu2, m, T1, T2 | dadi_sec_contact_sym_mig.py |
| Divergence in isolation, continuous asymmetrical secondary contact | nu1, nu2, m12, m21, T1, T2 | dadi_sec_contact_asym_mig.job |
| Divergence with ancient continuous symmetrical migration, isolation | nu1, nu2, m, T1, T2 | dadi_anc_sym_mig.py |
| Divergence with ancient continuous asymmetrical migration, isolation | nu1, nu2, m12, m21, T1, T2 | dadi_anc_asym_mig.py |
| Divergence in isolation, continuous symmetrical secondary contact, size change | nu1a, nu2a, nu1b, nu2b, m, T1, T2 | dadi_sec_contact_sym_mig_size.py |
| Divergence in isolation, continuous asymmetrical secondary contact, size change | nu1a, nu2a, nu1b, nu2b, m12, m21, T1, T2 |
| Divergence with ancient continuous symmetrical migration, isolation with size change | nu1a, nu2a, nu1b, nu2b, m, T1, T2 |
| Divergence with ancient continuous asymmetrical migration, isolation with size change | nu1a, nu2a, nu1b, nu2b, m12, m21, T1, T2 |
| Divergence in isolation, continuous symmetrical secondary contact, isolation | nu1, nu2, m, T1, T2, T3 | dadi_sec_contact_sym_mig_three_epoch.py |
| Divergence in isolation, continuous asymmetrical secondary contact, isolation | nu1, nu2, m12, m21, T1, T2, T3 | dadi_sec_contact_asym_mig_three_epoch.py |
| Divergence in isolation, continuous symmetrical secondary contact with size change, isolation | nu1a, nu2a, nu1b, nu2b, m, T1, T2, T3 |
| Divergence in isolation, continuous asymmetrical secondary contact with size change, isolation | nu1a, nu2a, nu1b, nu2b, m12, m21, T1, T2, T3 |

### All SNPs
Created folded SFS for observed data and bootstrapped variants in `dadi`. Followed the manual's [2D demography example](https://dadi.readthedocs.io/en/latest/examples/basic_workflow/basic_workflow_2d_demographics/).
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
wdir = "/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/demography/dadi/all_snps"
name_stem = "globMR_cordMR_ns32-16"

# Parse the VCF file to generate a data dictionary
datafile = "/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/04.filter_snps/all_fil_biallelic.vcf"
dd = dadi.Misc.make_data_dict_vcf(datafile, "/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/demography/stairway_plot/esfs_poplist.txt")

# Extract the spectrum for MR E. globulus and MR E. cordata from that dictionary, with 
# E. globulus projected down to 32 of 40 and E. cordata projected down to 16 of 20.
# haplotype counts at which the number of segregating sites was maximized
pop_ids, ns = ["glob_MR", "cord_MR"], [32, 16]
fs = dadi.Spectrum.from_data_dict(dd, pop_ids, ns, polarized = False) # folded since ancestral alleles not known

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