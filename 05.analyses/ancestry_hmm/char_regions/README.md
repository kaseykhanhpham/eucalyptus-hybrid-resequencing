# Characterize intervals of interest from Ancestry_HMM

## Delimit Introgressed Regions

Extracted introgressed regions identified by `Ancestry_HMM` for each sample using a custom R script. Introgressed regions were defined as:
* Containing no more than 750 bp in a row with less than 0.90 posterior probability of homozygous _E. cordata_ ancestry
* Containing at least one variant with a > 0.95 posterior probability of homozygous _E. cordata_ ancestry

```bash
# Performed in UFRC queue system; see get_intr_regions.job for more details.
# Resources used: 1.59 Mb, 1 min

module load R/4.2
LIST_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/ancestry_hmm"
POST_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/ancestry_hmm/posteriors"
SCRIPT_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/scripts"

while read NAME
do
    echo doing "$NAME"
    Rscript "$SCRIPT_DIR"/get_intr_bed.r "$POST_DIR"/"$NAME".posterior 0.95 0.90 750 0.85 ind_region_beds/"$NAME"_ahmm_regions.bed
done < "$LIST_DIR"/Eglobulus_MR.txt
```

For the regions on Chromosomes 6, 7, and 8 with many individuals sharing regions, got the intersection and union of regions between the samples.

Included the larger region for the first introgressed region in Chromosome 7 in the intersect set because all but two samples had it.

## genome scan statistics
Checked pi, dxy, and f-stats for intervals.

| Chromosome | Start (bp)  | End (bp)    | pi          | dxy (cord)   | dxy (glob pure) | fDm          | df           |
| ---------- | ----------- | ----------- | ----------- | ------------ | --------------- | ------------ | ------------ |
| Chr06      | 13,159,909  | 13,180,364  | 0.0576      | 0.0493       | 0.0772          | 0.217        | 0.0926       |
| Chr07      | 230,577     | 237,358     | N/A         | N/A          | N/A             | 0.157        | 0.0877       |
| Chr07      | 278,051     | 284,832     | N/A         | N/A          | N/A             | 0.183        | 0.0662       |
| Chr07      | 427,255     | 447,600     | 0.0584      | 0.0557       | 0.0636          | -0.0936      | -0.0143      |
| Chr07      | 7,074,553   | 7,123,386   | 0.0301      | 0.0251       | 0.0314          | 0.0879       | 0.0194       |
| Chr08      | 22,454,645  | 22,486,243  | 0.0906      | 0.115        | 0.0826          | 0.109        | 0.0263       |

All segments except for Chromosome 7 427255-447600bp generally fit expectations under introgression, though I haven't tested against the actual distributions of each statistic. There are several regions near the edge of the Chromosome 7 arm that did not have enough information to calculate pi/dxy, which is somewhat concerning and needs to be looked into.

Plotted distribution of pi, dxy, fDm, df and average value of those statistics within introgressed windows.

```R

```