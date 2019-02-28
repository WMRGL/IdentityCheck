# Downloading VCFs from GEL
VCF files for specified patients are downloaded from openCGA.

`run_download_vcfs.py` is the main script for this process. It looks for a file named "input.txt" at the storage location specified in the config file. This file should be a tab delimited text file consiting of GEL participant IDs and the local Lab patient ID for that sample:

`<GEL participant ID>  <Lab patient ID>`

The scpit converts this to a dictionary consisting of GEL participant ID keys and local lab patient ID values.

This dictionary is used to download the VCF files from the CIP-API and openCGA using many of the functions developed in the JellyPy project: https://github.com/NHS-NGS/JellyPy/tree/master/pyCIPAPI (The scripts used have been copied to this project and referenced).

The VCF files are identified and downloaded.
The following storage directories are expected within the storage location specified in the config file:
* <storage_location>/VCFs
* <storage_loction>/Intersected-VCFs

The VCFs are converted to bgzip format and tabixed to allow bcftools to perform an intersect with the text file of SNP locations.
The following files are expected at within the storage location specified in the config file:
* <storage_location>/hg19.nochr.txt
* <storage_location>/hg38.nochr.txt

These text files contain the genomic coordinates of the SNP locations for builds GRCh37 and GRCh38.

The intersected VCF files for each sample can be found within the Intersected-VCFs folder. The files will be named with the local lab patient ID number.
