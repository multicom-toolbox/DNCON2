#!/bin/bash
touch jhm-e-0.running
echo "running jackhmmer job jhm-e-0.."
/home/badri/hmmer-3.1b2-linux-intel-x86_64/binaries/jackhmmer --cpu 2 -N 5 -E 1 -A jhm-e-0.ali T0900.jh.fasta /home/badri/databases/uniref/uniref90pfilt > jhm-e-0-jackhmmer.log
if [ ! -f "jhm-e-0.ali" ]; then
   mv jhm-e-0.running jhm-e-0.failed
   echo "jackhmmer job jhm-e-0 failed!"
   exit
fi
/home/badri/DNCON2/scripts/reformat.pl -l 1500 -d 1500 sto a3m jhm-e-0.ali jhm-e-0.a3m
egrep -v "^>" jhm-e-0.a3m | sed 's/[a-z]//g' > jhm-e-0.aln
if [ -f "jhm-e-0.aln" ]; then
   rm jhm-e-0.running
   rm jhm-e-0-jackhmmer.log
   rm jhm-e-0.ali
   echo "jackhmmer job jhm-e-0 done."
   exit
fi
echo "Something went wrong! jhm-e-0.aln file not present!"
mv jhm-e-0.running jhm-e-0.failed
