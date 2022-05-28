# Metadata
This directory contains all files with sample information produced from sampling, sequencing, and analysis, as well as sample lists for ease of coding. Files are organized into subdirectories based on the step of the project they came from.


## Field Sampling

| Filename                                | Data Source               | Description                                  |
| --------------------------------------- | ------------------------- | -------------------------------------------- |
| 40samples_KaseyUF_updateDec2020.txt.kmz | Rebecca Jones, U Tasmania | Google Earth file with provenance of samples |
| 40samples_KaseyUF_updateDec2020.xlsx    | Rebecca Jones, U Tasmania | Spreadsheet of provenance and genotype data for samples |


## RAPiD Sequencing

| Filename                                | Data Source               | Description                                  |
| --------------------------------------- | ------------------------- | -------------------------------------------- |
| UFL_014801_Normal_Plate_Layout_22-02-21_GOOD.xlsx | Ariane Gélinas Marion, U Tasmania | Plate layout of samples sent to RAPiD Genomics for sequencing |
| UFL_014801_SampleSheet_june08.csv  | RAPiD Genomics | Table with sample codes and barcodes for run 2 of sequencing |
| UFL_014801_SampleSheet_may04.csv   | RAPiD Genomics | Table with sample codes and barcodes for run 1 of sequencing |


## Sequence Analysis

| Filename                                | Data Source               | Description                                  |
| --------------------------------------- | ------------------------- | -------------------------------------------- |
| chromosome_list.txt  | _Eucalyptus grandis_ genome | List of chromosome names in reference genome for code looping |
| sample_ids.txt          | sample_sequencing_metadata_all.xlsx | List of RAPiD ID numbers for each sample sequenced |
| sample_sequencing_metadata.csv | sample_sequencing_metadata_all.xlsx | Table of compiled metadata from all sources for each sequencing run file provided by RAPiD |
| sample_sequencing_metadata_all.xlsx | Filenames and headers of FASTQ files sent by RAPiD, sample sheets sent by RAPiD | Metadata for each FASTQ file, manually extracted and compiled. Files are split by sample, run, and direction. |
| sample_spp_table.csv                    | sample_sequencing_metadata_all.xlsx, 40samples_KaseyUF_updateDec2020.xlsx | Table translating RAPiD sample codes to U Tasmania accession numbers and species designation |
| seq_ids.txt | sample_sequencing_metadata_all.xlsx | List of all IDs unique to sample and sequencing run for code looping |

### Columns of note in sample_sequencing_metadata_all.xlsx
* **Filename**: The raw name of the FASTQ file sent by RAPiD Genomics
* **Sample**: RAPiD's ID for each _sample_, each of which was subject to two runs of sequencing.
* **RunSample**: RAPiD's ID for each _run_, unique for each sample and run, should identify each individual library prep.
* **Customer_Code**: Original accession numbers for each sample in U Tasmania data
