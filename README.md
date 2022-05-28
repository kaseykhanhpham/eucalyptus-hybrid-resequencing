# Eucalyptus Hybrid Resequencing
Pipeline for analysis of project resequencing introgressant Eucalyptus globulus and Eucalyptus cordata in Tasmania.

Directories are labeled in the order steps were performed. Computation performed on UF's HiperGator supercomputing cluster (UFRC), so job files are included where relevant.

## References
* The whole-genome long read assembly of _Eucalyptus grandis_: [ASM165458v1](https://www.ncbi.nlm.nih.gov/assembly/GCF_016545825.1/)
* The [annotation of ASM165458v1](https://www.ncbi.nlm.nih.gov/genome/annotation_euk/Eucalyptus_grandis/102/)
* The chloroplast assembly of _Eucalyptus globulus_: [AY780259.1](https://www.ncbi.nlm.nih.gov/nuccore/AY780259.1/)

## Python environment
Started using conda environments at the beginning of data analysis. 

On UFRC:
```bash
ENV_DIR="/blue/soltis/kasey.pham/conda/envs"
conda create -yp "$ENV_DIR"/euc_hyb_reseq python=3.10 pip numpy scipy pandas plotnine
```

On local computer:
```bash
conda create -n euc_hyb_reseq python pip numpy scipy pandas
pip install pong plotnine
```