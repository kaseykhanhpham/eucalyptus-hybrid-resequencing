# Annotate _E. globulus_ X46 Genome
The repeat library and MAKER pipeline used here is from [Daren Card's annotation pipeline for the boa constrictor genome](https://gist.github.com/darencard/bb1001ac1532dd4225b030cf0cd61ce2).

## Mask repeats
Built [`RepeatModeler`](http://www.repeatmasker.org/RepeatModeler/) database for assembly.

```bash
module load repeatmodeler/2.0

ASSEMB_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/refs/Eglobulus_genome_X46"
WDIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/annot_genome/repeat_lib/denovo"

ln -s "$ASSEMB_DIR"/EGLOB-X46.v1.0.fa "$WDIR"/EGLOB-X46.v1.0.fa
BuildDatabase -name eglob_x46 -engine ncbi EGLOB-X46.v1.0.fa
```

Ran `RepeatModeler` to make _de novo_ repeat library.
```bash
# Run in UFRC queue system, see repmod.job for more details.
# Resources used: 13 Gb, 1 day

module load repeatmodeler/2.0
RepeatModeler -pa 11 -engine ncbi -database eglob_x46 2>&1
```

Identified pre-identified repeats for _Eucalyptus_ from [Repbase](https://www.girinst.org/repbase).

```bash
# Run in UFRC queueing system, see repbase.job for more details.
# Resources: 11.6 Gb, 1.5 hrs

module load repeatmasker/4.0.9
ASSEMB_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/refs/Eglobulus_genome_X46"

RepeatMasker -pa 12 -e ncbi -species eucalyptus -dir eucalyptus "$ASSEMB_DIR"/EGLOB-X46.v1.0.fa
```

Masked _de novo_ repeats from `RepeatModeler` from `RepeatMasker` RepBase-masked genome.

```bash
# Run in UFRC queueing system, see repmask.job for more details.
# Resources: 12.4 Gb, 2.5 hrs

module load repeatmasker/4.0.9
REPMOD_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/repeat_lib/denovo/RM_100984.WedMay51749352021"
REPBASE_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/repeat_lib/repeatmasker/eucalyptus"

RepeatMasker -pa 12 -e ncbi -lib "$REPMOD_DIR"/consensi.fa.classified -dir denovo_mask "$REPBASE_DIR"/EGLOB-X46.v1.0.fa.masked
```