# CORONATATOR
CORONAvirus annoTATOR, version 1.0.0

# Brief tool guide of CORONATATOR
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
BPcalled：intermediate output, from breakpoint identification procedure. \
(read_id, is_forward, is_reverse, is_supplementary_alignment, mapping_start, CIGAR, CIGAR_summary, breakpoint_position, breakpoint_position_of_mate, context_of_breakpoint) \
\
SARS2.c.sg: canonical sgRNA called, final outputs \
(project_id, sample_id, strand, breakpoint1, breakpoint2, context_of_breakpoint, count) \
\
SARS2.nc.sg: non-canonical sgRNA called, final outputs \
(project_id, sample_id, strand, breakpoint1, breakpoint2, context_of_breakpoint, count) \
\
SGRNAcalled_noncanonical: non-canonical sgRNA read info, intermediate output, unfiltered. \
(read_id, strand, breakpoint1, breakpoint2, context_of_breakpoint, CIGAR_summary) \
\


 
