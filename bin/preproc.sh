#!/usr/bash
#SAMTOOLS=/usr/bin/samtools
BCFTOOLS=/home/vin/anaconda3/bin/bcftools
TABIX=/home/vin/anaconda3/bin/tabix
ANNOTATOR=/data/LyuLin/Scripts/Download/vcf-annotator/vcf-annotator.py
#VHEADER=MT121215.1
VHEADER=$4
BAM=$1
REF=$2
#REF_GB=/nis/home/xxqing/ncov/ref/NC_045512.gb
REF_GB=$3
STRAND=/data/LyuLin/WorkDirectory/2020-4-23-preprocessing/strandness.pl
#$SAMTOOLS view -h $BAM  $VHEADER | $SAMTOOLS mpileup -t SP -t DP - > tmp.bcf
#$SAMTOOLS view -h $BAM  $VHEADER | $SAMTOOLS mpileup  -t SP -t DP -gf $REF - | $BCFTOOLS call -vmO z -o test.vcf.gz -

echo 'start'
samtools view -h $BAM $VHEADER | samtools mpileup -t SP -t DP -gf $REF - | $BCFTOOLS call -vcO z  - | $BCFTOOLS filter -O v -e 'QUAL<50' - | $BCFTOOLS view -g hom -Oz -o out.filter.vcf.gz
echo 'out.filter.vcf.gz generated'
$TABIX out.filter.vcf.gz
echo 'tabixed'
$BCFTOOLS consensus -f $REF out.filter.vcf.gz > consensus.fasta
echo 'consensus.fasta generated'
#cat $REF | /nis/home/xxqing/anaconda3/bin/bcftools consensus out.filter.vcf.gz > consensus.fasta
gzip -d -f out.filter.vcf.gz
python $ANNOTATOR --output annotated.vcf out.filter.vcf $REF_GB
echo 'gb file modified'
samtools view $BAM $VHEADER|perl $STRAND > strandness
#$SAMTOOLS view -h $BAM  $VHEADER | $SAMTOOLS mpileup  -t SP -t DP -gf $REF - | $BCFTOOLS call -vmO z  -o a.vcf.gz
#$BCFTOOLS filter -O v -o b_filter.vcf -e 'QUAL<50 ' --SnpGap 5 --set-GTs . a.vcf.gz
#$BCFTOOLS view b_filter.vcf -Oz -o out.filter.vcf.gz



