#!/bin/bash
touch jhm-1e-20.running
echo "running jackhmmer job jhm-1e-20.."
/home/badri/hmmer-3.1b2-linux-intel-x86_64/binaries/jackhmmer --cpu 2 -N 5 -E 1e-20 -A jhm-1e-20.ali T0900.jh.fasta /home/badri/databases/uniref/uniref90pfilt > jhm-1e-20-jackhmmer.log
if [ ! -f "jhm-1e-20.ali" ]; then
   mv jhm-1e-20.running jhm-1e-20.failed
   echo "jackhmmer job jhm-1e-20 failed!"
   exit
fi
/home/badri/DNCON2/scripts/reformat.pl -l 1500 -d 1500 sto a3m jhm-1e-20.ali jhm-1e-20.a3m
egrep -v "^>" jhm-1e-20.a3m | sed 's/[a-z]//g' > jhm-1e-20.aln
if [ -f "jhm-1e-20.aln" ]; then
   rm jhm-1e-20.running
   rm jhm-1e-20-jackhmmer.log
   rm jhm-1e-20.ali
   echo "jackhmmer job jhm-1e-20 done."
   exit
fi
echo "Something went wrong! jhm-1e-20.aln file not present!"
mv jhm-1e-20.running jhm-1e-20.failed
