#!/bin/bash

set -e
set -u

date >> test.out.tsv
/data/LyuLin/Scripts/NanoDecoder/nanodecoder_p3 -f /data/chenz/nanopore/20220125-BNP2193-P7-PAJ04283-sup.pass.fastq.gz -b /data/LyuLin/Scripts/NanoDecoder/3kbarcode.lst -m 2 -r 100 -t 5 >> test.out.tsv
date >> test.out.tsv
