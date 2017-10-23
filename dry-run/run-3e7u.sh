#!/bin/bash

dtime=$(date +%Y-%b-%d)
../dncon2-v1.0.sh ./input/3e7u.fasta ./output/3e7u-$dtime/ &> ./output/3e7u-$dtime.log &

echo "Running in background.."
echo "Check log file.."
