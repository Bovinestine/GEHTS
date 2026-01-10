# Prerequisite

The default setting is preparing a directory for a batch of GEHTS-chip.
It includes 2 code subdirectories: './files' and './files_for_onewell'.

To run this pipeline, two custom files are required: raw .tif images in './raw_image' and a .csv table for DLP particlemap.

You need to make a './raw_image' directory for saving partitioned .tif images per channel (including merged channel) per sequencing cycle of in situ sequencing raw images (recommend 'export' function in NIS-Elements Analysis software). 
The file name would be like 'base1xy001c1.tif' or 'base4xy400c6'.
For GEHTS-Chip format, base: 1~4, xy: 1~400, c: 1~6 (including A, C, G, nucleus, and merged)

Also, you need to make a .csv file (e.g. './chip1_for_processing_100.csv') that is a particlemap of drug laden particles (DLPs). You can make a particlemap using custom look-up table or DLP decoding tools (ref: Adv Sci (Weinh). 2019 Feb 6;6(3):1970014. doi: 10.1002/advs.201970014, Biomicrofluidics 16, 061101 (2022) https://doi.org/10.1063/5.0131733)

The version of CellProfiler must be 2.x. since the code modules were deprecated in the newer version.

# How to run

1. run './files_for_onewell/before_cell_profiler.m' by setting the current directory in the headup lines: 
cd('./')

2. run './files_for_onewell/Seq_Standard_v2.cpproj'


3. run './files_for_onewell/after_cell_profiler.m' by setting the directory to the current directory.
directory = './'

Then, you will find the raw gene count table './in_situ/well_gene_counts38.csv'.

