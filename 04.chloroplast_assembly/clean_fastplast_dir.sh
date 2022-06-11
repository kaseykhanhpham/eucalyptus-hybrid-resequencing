# bash script to further clean fastplast result directories
# usage: clean_fastplast_dir.sh [dir_name]

# move to given directory
echo "Cleaning $1"
cd $1

# clean 1_Trimmed_Reads
echo "Removing 1_Trimmed_Reads"
rm -r 1_Trimmed_Reads

# clean 2_Bowtie_Mapping
echo "Removing 2_Bowtie_Mapping"
rm -r 2_Bowtie_Mapping

# clean 3_Spades_Assembly
echo "Cleaning 3_Spades_Assembly"
cd 3_Spades_Assembly/spades_iter1
rm -r K* # remove results for each K run
rm -r misc
rm -r pipeline_state
rm before_rr.fasta

cd ../..

# clean 4_Afin_Assembly
echo "Cleaning 4_Afin_Assembly"
cd 4_Afin_Assembly
rm Angiosperm_Chloroplast_Genes* # remove cp gene database files
rm filtered_spades_contigs.fsa
rm *_afin_iter0.fa
rm *_afin_iter1.fa

cd ..

# clean 5_Plastome_Finishing
echo "Cleaning 5_Plastome_Finishing"
cd 5_Plastome_Finishing
# remove BLAST database files
rm *.nhr
rm *nin
rm *nsq

cd ..

# clean Coverage_Analysis
echo "Cleaning Coverage_Analysis"
cd Coverage_Analysis
rm core* # remove bowtie coredumps
rm *.fq # remove mapped reads record
rm *.sam # remove mapped SAM file

cd ..

echo "Done cleaning $1"
