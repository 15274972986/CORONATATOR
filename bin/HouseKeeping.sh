#!/bin/bash

#set -e
#set -u
#set -o pipefail

ls /data/LyuLin/WorkDirectory/2020-5-SARS2-BAM-NC | grep '^PRJ' | parallel bash Wrapper.sh {}
