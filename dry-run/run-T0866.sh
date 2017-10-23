#!/bin/bash

dtime=$(date +%Y-%b-%d)
../dncon2-v1.0.sh ./input/T0866.fasta ./output/T0866-$dtime/ &> ./output/T0866-$dtime.log &

echo "Running in background.."
echo "Check log file.."

