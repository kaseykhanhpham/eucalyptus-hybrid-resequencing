# _Eucalyptus_ Hybrid Resequencing
Pipeline for analysis of project resequencing introgressant _Eucalyptus globulus_ and _Eucalyptus cordata_ in Tasmania.

Directories are labeled in the order steps were performed. Computation performed on UF's HiperGator supercomputing cluster (UFRC), so job files are included where relevant.

**Collaborators**: Rebecca Jones, Ariane Gélinas Marion, Brad Potts, René Vaillancourt, Pam Soltis, Doug Soltis

**With helpful input from**: Zhe Cai, Matias Kirst, Brad Barbazuk


## Project Design
**Background:**

This project was conceived based off of previous work documenting _Eucalyptus globulus_ and _Eucalyptus cordata_ genetic diversity in Tasmania, Australia. [McKinnon et al. 2004](https://doi.org/10.1111/j.1365-294X.2004.02364.x) first documented individuals of _E. globulus_ which possessed a chloroplast haplotype otherwise unique to _E. cordata_ at the site Meehan Range, perhaps the result of hybridization and chloroplast capture. [McKinnon et al. 2010](https://doi.org/10.1111/j.1365-294X.2010.04579.x) then characterized these individuals of _E. globulus_ alongside _E. cordata_ using AFLP markers and found evidence of a small amount of nuclear introgression from _E. cordata_ to a few geographically close individuals of _E. globulus_.

**Sampling:**

We now aim to characterize the regions of nuclear introgression in these previously-identified _E. globulus_ individuals with whole-genome resequencing. This project samples 10 individuals of _E. cordata_, 20 individuals of introgressant _E. globulus_ (either containing introgressed AFLP + chloroplast or just chloroplast markers), and 10 individuals of "pure" _E. globulus_ from outside Meehan Range. My UTas collaborators identified individuals, collected samples, and extracted DNA.


## References
* The whole-genome long read assembly of _Eucalyptus globulus_ (in prep): EGLOB-X46
* The annotation of the EGLOB-X46 (in prep)
* The chloroplast assembly of _Eucalyptus globulus_: [AY780259.1](https://www.ncbi.nlm.nih.gov/nuccore/AY780259.1/)


## Python Environment
Started using `conda` environments at the beginning of data analysis. 

On UFRC:
```bash
ENV_DIR="/blue/soltis/kasey.pham/conda/envs"
conda create -yp "$ENV_DIR"/euc_hyb_reseq python=3.10 pip numpy scipy pandas plotnine
pip install egglib
```

On local computer:
```bash
conda create -n euc_hyb_reseq python pip numpy scipy pandas
pip install pong plotnine
```