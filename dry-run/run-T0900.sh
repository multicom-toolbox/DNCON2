#!/bin/bash

dtime=$(date +%Y-%b-%d)
../dncon2-v1.0.sh ./input/T0900.fasta ./output/T0900-$dtime/ &> ./output/T0900-$dtime.log &

echo "Running in background.."
echo "Check log file.."

