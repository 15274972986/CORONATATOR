#!/bin/bash

#set -e
#set -u
#set -o pipefail

CORONATATOR=/data/LyuLin/SARS2/WorkDirectory/CORONATATOR
PREPROC=$CORONATATOR/bin/preproc.sh
CALLBP=$CORONATATOR/bin/s1.callBP.pl
CALLSG=$CORONATATOR/bin/s2.callSgRNA.pl
BAM_PATH=$CORONATATOR/input
META_PATH=$CORONATATOR/bin/Metadata.csv
REF_GENOME_PATH=$CORONATATOR/Refgenome
REF_GBK_PATH=$CORONATATOR/RefGBK
POSTPROC=$CORONATATOR/bin/s3.postProc.pl
OUTPUT=$CORONATATOR/output
PRJNA=$1

cd $OUTPUT

if [ ! -d "./$PRJNA" ]
then
	mkdir -p $PRJNA
fi

echo "$PRJNA started..."

cd $PRJNA

for SRA in `ls $BAM_PATH/$PRJNA`
do
	if [ ! -d "./$SRA" ]
	then
	mkdir $SRA
fi

echo "$SRA started"

cd $SRA

if [ ! -f "$BAM_PATH/$PRJNA/$SRA/align.out.sorted.bai" ]
then
	samtools index -@ 16 $BAM_PATH/$PRJNA/$SRA/align.out.sorted.bam $BAM_PATH/$PRJNA/$SRA/align.out.sorted.bai
fi

SPECIES="`grep $SRA $META_PATH | cut -d ',' -f 12`"

echo "Species is $SPECIES"

header="`basename -s '.fasta' $REF_GENOME_PATH/$SPECIES/*.fasta`"

echo "Header is $header"

if [ `tail -13 $CORONATATOR/bin/species_LTRS.tsv | grep -w $SPECIES | cut -d ' ' -f 3 | wc -l` -eq 1  ]
then
	leader="`tail -13 $CORONATATOR/bin/species_LTRS.tsv | grep -w $SPECIES | cut -d ' ' -f 3`"
else
	leader=59
fi

echo "Leader is $leader"
echo "Generating consensus sequence....."
bash $PREPROC $BAM_PATH/$PRJNA/$SRA/align.out.sorted.bam $REF_GENOME_PATH/$SPECIES/*.fasta $REF_GBK_PATH/$SPECIES.gb $header
echo "Generating BP file......"
samtools view -@ 16 $BAM_PATH/$PRJNA/$SRA/align.out.sorted.bam | perl $CALLBP -r consensus.fasta -l $leader > BPcalled
echo "Generating SGRNA......"

sort BPcalled | perl $CALLSG > SGRNAcalled_noncanonical

PRJNUM="`echo "$PRJNA" | cut -d "-" -f 1`"

echo "PRJ ID is $PRJNUM"

SRX="`cat $META_PATH | grep "$SRA" | grep -o '[ESDCT]RX[0-9]*'`"

echo "SRX ID is $SRX"
echo "Generating sg file"

perl $POSTPROC SGRNAcalled_noncanonical -c -l $leader -p $PRJNUM -s $SRX > $SPECIES.c.sg
perl $POSTPROC SGRNAcalled_noncanonical -p $PRJNUM -s $SRX > $SPECIES.nc.sg

cd ..
done

exit

#5425 -l 61
#131656 -l 58
#BY140568 -l 81
#HK140714 -l 65
#EPI_ISL_410542 -l 59
#EPI_ISL_410544 -l 59
#EPI_ISL_410543 -l 47
#RaTG13 -l 59
#MERS -l 66
#SARS2 -l 74
#SARS -l 71
#160660 -l 59
#HKU1 -l 46
#HKU3 -l 69
