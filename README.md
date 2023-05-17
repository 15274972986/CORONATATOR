CORONATATOR
CORONAvirus annoTATOR, version 1.0.0
#####

#Brief tool guide of CORONATATOR
## 1. prepare *.bam file:
the bam file used for analysis should be sorted, and put in the directory in 'input' in form as 'input/PRJ_NAME/SAMPLE_NAME/align.out.sorted.bam'

## 2. edit Metadata.csv
fullfill the table containing sample info, in 'bin/Metadata.csv'. This is an example:
PRJTEST0001,SAMPTEST01,SAMPTEST01,,,,RNA-Seq,,,PAIRED,ILLUMINA,SARS2,in vivo,

The header of this table is 'PRJ,SRX,SRR,Spots,Bases,Average_spot_length,Assay,Selection,Library_source,Layout,Platform,Virus,in vivo/in vitro,Center_name'
info of first three column and column 12 must be provided(those are column 'PRJ,SRX,SRR' and 'Virus'), other info is not needed for pipeline.
note: 'Virus' column could be SARS, MERS, SARS2 etc, all species name are available for looking up in 'bin/species_LTRS.tsv'

## 3. run script
cd bin/
bash Wrapper.sh PRJ_NAME

## 4. field of output files (use files in 'output/PRJTEST0001/SAMPTEST01/' as example)
