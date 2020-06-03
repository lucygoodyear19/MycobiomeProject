# Mycobiome Project Generalised Scripts Repository

&nbsp;

*Author: Lucy Goodyear*  
*Created: 04/02/20*

This repository contains all scripts required for prcoessing and analysing mycobiome data.

&nbsp;

## Requirements

All code has been created for Mac so there may be a few differences in commands with respect to Linux.

File architecture must follow this structure:

```
.
└── Analysis
    ├── Generalised_Scripts
    ├── Results  
    └── Runs_Countries
        ├── CostaRica_Ecuador_2017
        │   ├── DADA2_Results
        │   │   ├── CostaRica
        │   │   └── Ecuador
        │   │       ├── DADA2_args.R   
        │   │       └── HPC_outputs
        │   └── Esto_Results
        └── Taiwan_Vietnam_2016
        │   ├── DADA2_Results
        │   │   ├── Taiwan
        │   │   └── Vietnam
        │   │       ├── DADA2_args.R   
        │   │       └── HPC_outputs
        │   └── Esto_Results
 ```

 There is actually a lot of flexibility with the structure. The only requirement is that the distances between subdirectories must be the same and the countries within the run directories must be named by country only.