# Chloroplast Assembly

## Run `FastPlast`

### First round of runs

**Create job to run `FastPlast`**:
Wrote python script to generate a job file which runs [`FastPlast`](https://github.com/mrmckain/Fast-Plast) on all samples, using trimmed reads from both sequencing runs for each sample:

```bash
module load python/3.8
SCRIPT_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/scripts"

python "$SCRIPT_DIR"/generate_fp_job.py
```

**Run `FastPlast`:**
```bash
# Run via job on UFRC, see fastplast1.job, fastplast2.job, fastplast3.job, fastplast4.job for details
# Showing command for 1 out of 40 samples as commands are essentially the same for all, save for sample IDs.
# Resources used: 5 Gb, 9 hrs
# Bowtie2 always throws an error and coredumps at some point in mapping the subset of reads, but I am getting good looking chloroplasts anyways and it is not because my trimmed reads are corrupted (they map with Bowtie2 to the reference plastome just fine)... So I am not going to try to fix it.
# WE01 HAD A WEIRD GENE RECOVERY -- REDO

READ_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/reads/trimmed_reads"
REF_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/refs/organelle"

perl $HPC_FASTPLAST_DIR/Fast-Plast/fast-plast.pl -1 "$READ_DIR"/S320_R1_paired.fq,"$READ_DIR"/S1_R1_paired.fq -2 "$READ_DIR"/S320_R2_paired.fq,"$READ_DIR"/S1_R2_paired.fq -n WA01 --subsample 30000000 --threads 7 --user_bowtie "REF_DIR"/AY780259.1 --clean deep --skip trim --coverage_analysis --min_coverage 10
```

### WE01 Debugging

Sample WE01 assembled strangely and only recovered 40% of chloroplast genes in its final assembly, even though after FastPlast step 4 (`afin` assembly), 90% of chloroplast genes were recovered. I aligned the `afin` assemblies for WE01 and WB01 (another _E. cordata_ sample with an unproblematic assembly) and the identified regions for each (lsc, irb, ssc).

```bash
module load mafft/7.490
# Run via job in UFRC, see mafft_WB01WE01.job for details
# Resources used: 

mafft --auto --thread 8 afin_iter2_WB01WE01.fa > afin_iter2_WB01WE01_aligned.fa
mafft --auto --thread 8 final_lsc_WB01WE01.fa > final_lsc_WB01WE01_aligned.fa
mafft --auto --thread 8 final_irb_WB01WE01.fa > final_irb_WB01WE01_aligned.fa
mafft --auto --thread 8 final_ssc_WB01WE01.fa > final_ssc_WB01WE01_aligned.fa
```

From alignment, I could see that the large single copy region (LSC) was missing a large chunk, so I traced back to step 5 in `FastPlast` and took a look at `WE01_regions_split0.fsa`. It appeared that the assembly identified the following pattern: SC - IR - SC - IR - SC with the middle single copy region consistent with the size of the small single copy in other assemblies. So I assume that the large single copy region was split in two and half removed in further iterations. I then put the file through the rest of the plastome finishing step as advised by [the `FastPlast` troubleshooting documentation](https://github.com/mrmckain/Fast-Plast/blob/master/Troubleshooting.md). 

```bash
module load fastplast/1.2.8
WDIR="/blue/soltis/kasey.pham/euc_hyb_reseq/cp_assembly/WE01"

cd "$WDIR"/5_Plastome_Finishing
perl $HPC_FASTPLAST_DIR/Fast-Plast/bin/orientate_plastome_v.2.0.pl WE01_regions_split0.fsa WE01_regions_split0.fsa.blastn WE01

cd "$WDIR"/Final_Assembly
cp "$WDIR"/5_Plastome_Finishing/WE01_CP_pieces.fsa .
cp "$WDIR"/5_Plastome_Finishing/WE01_FULLCP.fsa .

rm -r "$WDIR"/Coverage_Analysis
```

Ran coverage analysis on re-oriented plastome assembly.

```bash
# Run via job in UFRC, see coverage_WE01.job for details
# Resources used: 

module load fastplast/1.2.8
READS_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/reads/trimmed_reads"
WDIR="/blue/soltis/kasey.pham/euc_hyb_reseq/cp_assembly"

perl $HPC_FASTPLAST_DIR/Fast-Plast/fast-plast.pl -1 "$READS_DIR"/S324_R1_paired.fq,"$READS_DIR"/S5_R1_paired.fq -2 "$READS_DIR"/S324_R2_paired.fq,"$READS_DIR"/S5_R2_paired.fq -n WE01 --subsample 30000000 --threads 7 --only_coverage "$WDIR"/WE01/Final_Assembly/WE01_FULLCP.fsa --skip trim --min_coverage 10
```

I noticed during alignment that the SSC and LSC regions for WE01 needed to be reverse-complemented. (Referenced [this](https://www.biostars.org/p/14614/) Biostars answer)

```bash
# working directory: /blue/soltis/kasey.pham/euc_hyb_reseq/cp_assembly/WE01/Final_Assembly
module load biopieces/2.0
mkdir $BP_DATA $BP_TMP $BP_LOG

mv WE01_CP_pieces.fsa WE01_CP_pieces_raw.fsa
# reverse complement just the LSC and SSC regions
read_fasta -i WE01_CP_pieces_raw.fsa | grab -p lsc,ssc | reverse_seq | complement_seq | write_fasta -x -o WE01_CP_pieces.fsa
# extract the IRB region as-is
read_fasta -i WE01_CP_pieces_raw.fsa | grab -p irb | write_fasta -x -o irb_temp.fsa
# merge IRB sequence with the reverse-complemented sequences
cat irb_temp.fsa >> WE01_CP_pieces.fsa
rm irb_temp.fsa
```

### Housekeeping
```bash
SCRIPTS_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/scripts"
for NAME in WA01 WA02 WA03 WA04 WA05 WB01 WB02 WB03 WB04 WB05
do
    "$SCRIPTS_DIR"/clean_fastplast_dir.sh $NAME
done

SCRIPTS_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/scripts"
for NAME in WC01 WC02 WC03 WC04 WC05 WD01 WD02 WD03 WD04 WD05
do
    "$SCRIPTS_DIR"/clean_fastplast_dir.sh $NAME
done

SCRIPTS_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/scripts"
for NAME in WE01 WE02 WE03 WE04 WE05 WF01 WF02 WF03 WF04 WF05
do
    "$SCRIPTS_DIR"/clean_fastplast_dir.sh $NAME
done

SCRIPTS_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/scripts"
for NAME in WG01 WG02 WG03 WG04 WG05 WH01 WH02 WH03 WH04 WH05
do
    "$SCRIPTS_DIR"/clean_fastplast_dir.sh $NAME
done
```