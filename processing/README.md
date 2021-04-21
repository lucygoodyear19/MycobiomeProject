# Mycobiome Project Generalised Processing Scripts Repository

&nbsp;

*Author: Luke Goodyear*  
*Created: 04/02/20*

This repository contains all scripts required for prcoessing mycobiome data, which can then be imported into analysis scripts. Each script has been generalised so can easily be used for different datasets.

&nbsp;

## Requirements

All code has been created for Mac so there may be a few differences in commands with respect to Linux.

Linux packages:

csvkit (v.1.0.4)  
cutadapt (v.2.8)

R packages:

dada2 (v.1.14.1)   
ShortRead (v.1.44.3)  
Biostrings (v.2.54.0)  
phyloseq (v.1.30.0)  
dplyr (v.1.0.1)  


File architecture must a similar structure to the below:

```
.
└── analysis
    ├── scripts
    ├── results  
    └── runs_countries
        ├── CostaRica_Ecuador_2017
        │   ├── DADA2_results
        │   │   ├── CostaRica
        │   │   └── Ecuador
        │   │       ├── DADA2_args.R 
        │   │       ├── filtering_args.R          
        │   │       └── HPC_outputs
        │   └── original_data
        └── Taiwan_Vietnam_2016
            ├── DADA2_results
            │   ├── Taiwan
            │   └── Vietnam
            │       ├── DADA2_args.R   
            │       ├── filtering_args.R   
            │       └── HPC_outputs
            └── original_data
 ```

 There is actually a lot of flexibility with the structure. The only requirement is that the distances between subdirectories must be the same and the countries within the run directories must be named by country only.


&nbsp;

## Workflow

### 0.1) [DADA2_args.R](#01-dada2_args)

R-script containing the arguments required to run DADA2_HPC_pipeline.R

### 0.2) [filtering_args.R](#02-filtering_args)

R-script containing the arguments required to run DADA2_data_filtering.R

### 1) [sample_seqs_prep.sh](#1-sample_seqs_prep)

Given a tarball containing data from an Illumina MiSeq run, this script processes the sequences into a form ready to be run through the DADA2 pipeline.

### 2) [DADA2_pipeline.R](#2-dada2_pipeline)

Processes seqeunce data to produce useful outputs for analysis, including taxa and abundance tables, quality profiles, and error rates.

### 3) [DADA2_data_filtering.R](#3-dada2_data_filtering)

Filters the results from the DADA2_HPC_pipeline R-script, removing mocks, negative and positive controls. Result is a phyloseq object.

### 4) [run_DADA2_HPC_pipleine.sh](#4-run_dada2_hpc_pipleine)

Runs DADA2_HPC_pipeline.R and the Data_filtering.R on the HPC. Each R-script is easily recognisable and can be hashed out as required. Running the whole bash script will result a phyloseq object containing otu_table, taxa_table, sample_data.

### 5) [phyloseq_merges.R](#5-phyloseq_merges)

Combines phyloseq objects via the command line so different combinations can be analysed.

### 6) [calculate_elevation.R](#6-calculate_elevation)

Given latitudes and longitudes in csv format, it calculates elevation using the geonames database.

&nbsp;



## Scripts

### 0.1) DADA2_args

#### *Requirements*

#### Arguments: 

None.

#### Packages:

None.

#### *How to use*

This script is stored in the relevent country folder and contains 11 variables:  

REV - reverse primer  
FWD - forward primer  
base_prefix - base name of sequences  
patFs - string pattern for forwards reads
patRs - string pattern for reverse reads
run  - run directory
country - country directory  
root_path - your path to the directory containing all runs  
path - path to directory containing sequences, using the root_path  
path2 - path to filtN directory, using the root_path  
path_out - path to directory where output will be sent to, using root_path   
cutadapt - path to cutadapt directory, using root_path  
unite.ref - path to UNITE database, using root_path 



&nbsp;



### 0.2) filtering_args

#### *Requirements*

#### Arguments: 

None.

#### Packages:

None.

#### *How to use*

This script is stored in the relevent country folder and contains 6 variables:  

run  - run directory
country - country directory  
root_path - your path to the directory containing all runs  
path_in - path to data input directory using root_path  
path_out - path to directory where output will be sent to, using root_path  
metadata_path - path to metadata file, using root_path  



&nbsp;



### 1) sample_seqs_prep

#### *Requirements*

Bash:

csvkit v. 1.0.4

To use on the HPC, **csvkit** must be added to a conda environment called ```bash_env```, which is loaded during the script. Please comment out these lines before running if you are running on a machine that doesn't use conda environments.

#### *Arguments*: 

Mandatory arguments:

1) root name of .tar.gz file containing data

2) csv file containing the following columns (as a minimum):

      i) Sample_ID (matching directory IDs)

      ii) Plate (number)

      - To generate the plate number in Mac Numbers, I use the following numbers function: 

      ```IF(IFERROR(SEARCH("SA*",I7_Index_ID IGF0009100_1)>0,0),1,IF(IFERROR(SEARCH("SB*",I7_Index_ID IGF0009100_1)>0,0),2, IF(IFERROR(SEARCH("SC*",I7_Index_ID IGF0009100_1)>0,0),3,4)))```

      - where:

      1) the IFERROR statement acts as ISNUMBER in excel.

      2) I7_Index_ID IGF0009100_1 is the cell in the I7_Index_ID column corresponding to the sample of interest, which contains a code the beginning of which corresponds to the plate number.

      iii) Origin (country matching country directory or Mock or PosC or NC)

      - The country name must contain no spaces and captial letters for new words, e.g. CostaRica.
      If the sample is a control, this column will contain "NC", "Mock" or "PosC" accordingly.

      - An example of the commmand to find this (based on Sample_ID):

      ```IF(IFERROR(SEARCH("*_V20*",Sample_Name IGF100881)>0,0),"Vietnam",IF(IFERROR(SEARCH("HND*",Sample_Name IGF100881)>0,0),"Honduras", IF(IFERROR(SEARCH("mock*",Sample_Name IGF100881)>0,0),"Mock",IF(IFERROR(SEARCH("pos*",Sample_Name IGF100881)>0,0),"PosC", IF(IFERROR(SEARCH("NC*",Sample_Name IGF100881)>0,0), "NTC", IF(IFERROR(SEARCH("NTC*",Sample_Name IGF100881)>0,0), "NTC","Pyrenees"))))))```

#### *How to use*

The script should be run from the directory containing the tar ball. I have saved this within the relevent **run** directory under a subdirectory **original_data** (as can be seen on file structure tree).

An example command to run the script:

```bash ../../../scripts/processing/sample_seq_prep.sh NorthEastAsia_IGFfolders ../plate_data```

#### *Outputs*

The script will output subdirectories in the **original_data** directory corresponding to each country present, each containing a subdirectory called **sample_seqs**, which contains all the fastq files for that country and all the fastq files for the mocks/posc/nc for each plate containing at least one sample from that country.

It is this **sample_seqs** directory that will be accessed in the DADA2 pipeline script below.



&nbsp;



### 2) DADA2_pipeline

#### *Requirements*

#### Arguments: 

DADA2_args.R 

#### Packages:

R:

dada2 (v.1.14.1)   
ShortRead (v.1.44.3)  
Biostrings (v.2.54.0)  

Linux:

cutadapt (v.2.8)

The HPC uses anaconda so I have installed all the above pacakges into a conda environment called ```R_processing_env```.

#### *How to use*

This is a script following the DADA2 ITS2 tutorial with some modifications to use only the forward reads and to run on HPC. Thanks to Phil Jervis for giving me his script as a basis and to the creators of the DADA2 pipeline (https://benjjneb.github.io/dada2/ITS_workflow.html). 

**Pipeline**

The script performs a number of sequential steps to determine community composition of samples and to assign taxonomy of fastq files:

1) Sequences with ambiguous bases are filtered out  

2) Presence of the primers specified in the arguments script is verified and their orientation is checked

3) Primers are removed using cutadapt tool due to the variable length of the ITS region (note cutadapt prints a lot to screen)

4) Quality profiles are plotted for the forward reads (only forward reads are fully processed due to consistent lower quality of reverse reads)

5) Reads are filtered by:
 - maximum number of "expected errors" allowed in a read
 (rather than averaged quality scores)
 - minimum length of 50bp to remove spurious, very low length sequences
 - truncate reads at the first instance of a quality score of less than 2 (rather than to a fixed length due to the variable length of the ITS region)                                    
 - removing PhiX, which is a common control used by illumina

6) Error rates are learned for each possible transition and are plotted (ignore ```not all sequences were the same length error```)

7) Identical reads are dereplicated

8) Divisive Amplicon Denoising Algorithm is applied to remove errors based on error rate and therefore infer sample community composition.

9) ASV (amplicon sequence variant) table is generated and chimeras are removed

10) Taxonomy is assigned using UNITE database

Also ensure that your sample name is separated from the rest of the fasta file name by an underscore, _, otherwise you will not get the sample names you want in the output. Else, rename your fasta files before running the script. For example, i had t rename some of my samples with the following ```rename 'CHN17_' 'CHN17-' *``` commandline command.

**Outputs**

1) quality profile.pdf 
      - shows the quiality profiles of the forward reads
2) read_counts_during_pipeline_steps.csv
      - view number of reads that made it through each stage of the pipeline:
      input, filtered, denoisedF, nonchim
3) error_rates.pdf
      - shows learned error rates
4) tax_table.txt
      - taxonomy table with ASVs as rows and taxonomic ranks as columns
5) abun_table_samplesinrows.txt
      - abundance table with samples as rows and ASVs as columns


Screen output (or .o file in the case of running on HPC) should be reviewed before proceeding, as should the first four outputs listed above. These all contain checks to ensure nothing has gone wrong during te pipeline.

No modifications should be required in order for this script to run with any R-script argument containing the correct variables and with the sequences in the correct format.



&nbsp;



### 3) DADA2_data_filtering

#### *Requirements*

#### Arguments: 

filtering_args.R

#### Packages:

R:
   
phyloseq (v.1.30.0)  
dplyr (v.1.0.1)  

#### *How to use*

Can be run from local commandline with the following command: 

```Rscript  ../../../../scripts/processing/DADA2_data_filtering.R filtering_args.R```

Filters by removing all controls and creates phyloseq object.

Troubleshooting:  

```Error in dimnames(x) <- dn : length of 'dimnames' [1] not equal to array extent```
This has occurred beacuse the taxa table and the abundance table get emptied of all samples, most likely because of a sample name issue, causing all samples to be filtered out.
Check line 107 to confirm sample_names in plate data match those in the abundance table (these may not match where we previsouly replaced _ with - (or othewrwise) to enable the truncation of sample names in the DADA2 pipeline)
```plates$Sample_Name <- gsub("_", "-", plates$Sample_Name)```

Another name change error I had 
``` Error in `[.data.frame`(plate, , (paste0("posC", plate_no))) : 
  undefined columns selected ```
This is part of the function that removes any negative controls and mock leakage into posC
Fixed by the following:
```plates$Sample_Name <- gsub("-A", "", plates$Sample_Name)``` where I had changed the original names of the mocks and posc to fit with the script (mock1-A to mock1) of the actual samples (so in the taxa df and the abundance df) but not in the plates csv.

A common warning message: ```Warning message: In estimate_richness(dada2, split = TRUE, measures = "Shannon"): The data you have provided does not have
any singletons. This is highly suspicious. Results of richness
estimates (for example) are probably unreliable, or wrong, if you have already trimmed low-abundance taxa from the data. We recommended that you find the un-trimmed data and retry.```
Sometimes I got the above warning but I hadn't removed the singletons or done any of that kind of trimming so I ignored it.

Output is two phyloseq objects. One keeps the ASV names as the full sequences, which is for merging and phylogenetic tree building, and the other renames to "ASV1, ASV2,..." for convenience and saves sequences under refseqs option of phyloseq option. The latter object is for analysis on that specific country/run.

&nbsp;



### 4) run_DADA2_HPC_pipeline

#### *How to use*

This script has been written to run on a PBS HPC so adjustments may need to be made to accommodate other machines.

The following command should be used to run the script, making sure to check your own releative paths:

```qsub -v "run=run, country=country" ../../../../scripts/processing/run_DADA2_HPC_pipeline.sh ```

Where ```run/country``` should be replaced with the relevent path. For example:

```qsub -v "run=FrenchGuiana_Seychelles_2016, country=FrenchGuiana" ../../../../../scripts/processing/run_DADA2_HPC_pipeline.sh```


Relative paths should be adapted for your own file structure or use absolute paths. I always run this from the relevant **HPC_outputs** directory so the output and error files are also located here.



&nbsp;



### 5) phyloseq_merges

#### *Requirements*

#### Arguments: 

Phyloseq objects to be merged (no maximum) should be used as command line arguments. ASV names must be full sequences.

#### Packages:

phyloseq (v.1.30.0)

#### *How to use*

This script merges phyloseq objects and subsets sample data by those variables specified in the script. Before running, please confirm the list of varibales within the script in the "Set up" section. Sequences of ASVs are saved into the refseqs slot and ASVs are renamed "ASV1, ASV2, ..." for convenience.

Outputs a single phyloseq object.



&nbsp;



### 6) calculate_elevation

#### *Requirements*

#### Arguments:

None

#### Packages:

rgbif
tidyr

#### *How to use*

Calculates elevations from the geonames database given latitudes and longitudes in csv format. 

Outputs two csvs, one is the original saved under "metadata_pre-elevs.csv" and the other is the new metadata csv with elevations included

Note, do not run this twice in a row or original data will be overwritten.
