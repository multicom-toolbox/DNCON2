#!/bin/bash -e
# Badri Adhikari, 9-9-2017
# The main script for making DNCON2 contact prediction

if [ $# != 2 ]; then
	echo "$0 <fasta> <output-directory>"
	exit
fi

ROOT=$(dirname $0)

$ROOT/scripts/dncon2-main.pl $1 $2
