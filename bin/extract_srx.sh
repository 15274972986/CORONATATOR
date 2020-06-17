#!/bin/bash

set -e
set -u
#set -o pipefail

META_PATH=/data/LyuLin/WorkDirectory/2020-4-23-preprocessing/metainfo

PRJ=$1
SSR=$2

if [ -f "$META_PATH/${PRJ}_meta.txt" ]
then
	if [ `grep $SSR $META_PATH/${PRJ}_meta.txt | wc -l` -eq 1 ]
	then
		cat $META_PATH/${PRJ}_meta.txt | sed '1d' | grep "$SSR" | grep -o "[SET]RX[0-9]*"
	else
		echo "SRX0000SRX"
	fi
else
	echo "SRX0000PRJ"
fi
