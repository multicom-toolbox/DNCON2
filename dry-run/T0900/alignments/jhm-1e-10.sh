#!/bin/bash
touch jhm-1e-10.running
echo "running jackhmmer job jhm-1e-10.."
/home/badri/hmmer-3.1b2-linux-intel-x86_64/binaries/jackhmmer --cpu 2 -N 5 -E 1e-10 -A jhm-1e-10.ali T0900.jh.fasta /home/badri/databases/uniref/uniref90pfilt > jhm-1e-10-jackhmmer.log
if [ ! -f "jhm-1e-10.ali" ]; then
   mv jhm-1e-10.running jhm-1e-10.failed
   echo "jackhmmer job jhm-1e-10 failed!"
   exit
fi
/home/badri/DNCON2/scripts/reformat.pl -l 1500 -d 1500 sto a3m jhm-1e-10.ali jhm-1e-10.a3m
egrep -v "^>" jhm-1e-10.a3m | sed 's/[a-z]//g' > jhm-1e-10.aln
if [ -f "jhm-1e-10.aln" ]; then
   rm jhm-1e-10.running
   rm jhm-1e-10-jackhmmer.log
   rm jhm-1e-10.ali
   echo "jackhmmer job jhm-1e-10 done."
   exit
fi
echo "Something went wrong! jhm-1e-10.aln file not present!"
mv jhm-1e-10.running jhm-1e-10.failed
