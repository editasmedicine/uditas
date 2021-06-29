UDiTaS v1.1
===========

Overview
--------

UDiTaS(TM) stands for UniDirectional Targeted Sequencing, a novel sequencing method useful for measuring small indels as well as
structural rearrangements, like translocations, in a single reaction.

See details of the method in Giannoukos et al. BMC Genomics (2018) 19:212, https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-018-4561-9


Systems Requirements
--------------------

UDiTaS has been tested in python 3.9 and requires the python packages and tools listed in the file uditas_env.yml

The code requires setting up two environmental variables

`BOWTIE2_INDEXES` contains the location of the bowtie2 index files, typically hg38.1.bt2, hg38.2.bt2, etc.

`GENOMES_2BIT` contains the location of the 2bit files for the genomes used in the analysis, eg hg38.2bit

To test the code create a virtual python environment with

`conda env create -f uditas_env.yml`

then activate using

`source activate uditas_env`

To install uditas as an executable run

`python setup.py install`

To test the installation run

`pytest`

This will process the test data in

```
data/fig2a
data/fig2b
data/fig2c
```

These data are subsamples of the data displayed in Fig 2 of the paper.

Usage
-----
`uditas` is the main command to launch the UDiTaS analysis. The required argument is: `dir_sample`, the path of the directory with the data to be analyzed.

`dir_sample` should contain the fastq.gz files for R1, R2, I1 and I2. Used in the demultiplexing step.

`dir_sample` should also contain the analysis sheet `sample_info.csv` containing the description of all the samples, with their barcodes and guides used. See examples in the folder `data`

Once the setup has been run, the code can be run as

`uditas ./data/fig2a`

The full list of options are:

```
usage: uditas [-h] [-folder_genome_2bit FOLDER_GENOME_2BIT]
              [-skip_demultiplexing SKIP_DEMULTIPLEXING]
              [-skip_trimming SKIP_TRIMMING]
              [-skip_genome_local_alignment SKIP_GENOME_LOCAL_ALIGNMENT]
              [-skip_genome_global_alignment SKIP_GENOME_GLOBAL_ALIGNMENT]
              [-process_amplicon PROCESS_AMPLICON]
              [-skip_amplicon_global_alignment SKIP_AMPLICON_GLOBAL_ALIGNMENT]
              [-check_plasmid_insertions CHECK_PLASMID_INSERTIONS]
              [-skip_plasmid_alignment SKIP_PLASMID_ALIGNMENT] [-ncpu NCPU]
              [-window_size WINDOW_SIZE]
              [-default_amplicon_window_around_cut DEFAULT_AMPLICON_WINDOW_AROUND_CUT]
              [-min_MAPQ MIN_MAPQ] [-min_AS MIN_AS]
              [-process_AMP_seq_run PROCESS_AMP_SEQ_RUN]
              dir_sample

Process UDiTaS data

positional arguments:
  dir_sample            Directory with the sample to be processed

optional arguments:
  -h, --help            show this help message and exit
  -folder_genome_2bit FOLDER_GENOME_2BIT
                        Folder containing the 2bit file(s) with the reference
                        genome being used (default: GENOMES_2BIT)
  -skip_demultiplexing SKIP_DEMULTIPLEXING
                        Skip demultiplexing? Options: 0, 1 (skip) (default: 0)
  -skip_trimming SKIP_TRIMMING
                        Skip adapter trimming? Options: 0, 1 (skip) (default:
                        0)
  -skip_genome_local_alignment SKIP_GENOME_LOCAL_ALIGNMENT
                        Skip genome-wide local alignment? Options: 0 , 1
                        (skip) (default: 1)
  -skip_genome_global_alignment SKIP_GENOME_GLOBAL_ALIGNMENT
                        Skip genome-wide global alignment? Options: 0 , 1
                        (skip) (default: 0)
  -process_amplicon PROCESS_AMPLICON
                        Select row number (0-based) of amplicon to process,
                        set to all to process all amplicons (default: all)
  -skip_amplicon_global_alignment SKIP_AMPLICON_GLOBAL_ALIGNMENT
                        Skip amplicon global alignment? Options: 0, 1 (skip)
                        (default: 0)
  -check_plasmid_insertions CHECK_PLASMID_INSERTIONS
                        Check for plasmid insertions. Options: 0 (skip), 1
                        plamid_name and plasmid_sequence required in
                        sample_info.csv (default: 1)
  -skip_plasmid_alignment SKIP_PLASMID_ALIGNMENT
                        Skip plasmid alignment? Note, just alignment. Counts
                        still evaluated. Options: 0, 1 (skip) (default: 0)
  -ncpu NCPU            Number of CPUs to use (default: 4)
  -window_size WINDOW_SIZE
                        Window size around cut sites used to grab UDiTaS reads
                        (default: 15)
  -default_amplicon_window_around_cut DEFAULT_AMPLICON_WINDOW_AROUND_CUT
                        Window size around cut sites used to create amplicons
                        (default: 1000)
  -min_MAPQ MIN_MAPQ    Minimum mapping quality to include a read (default: 5)
  -min_AS MIN_AS        Minimum alignment score to include a read (default:
                        -180)
  -process_AMP_seq_run PROCESS_AMP_SEQ_RUN
                        Set to 1 to process an AMP-seq run using GUIDE-seq
                        adapters (default: 0)
```
