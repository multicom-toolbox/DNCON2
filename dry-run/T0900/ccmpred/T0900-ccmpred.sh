#!/bin/bash
touch ccmpred.running
echo "running ccmpred .."
/home/badri/CCMpred/bin/ccmpred -t 8 T0900.aln T0900.ccmpred > ccmpred.log
if [ -s "T0900.ccmpred" ]; then
   mv ccmpred.running ccmpred.done
   echo "ccmpred job done."
   exit
fi
echo "ccmpred failed!"
mv ccmpred.running ccmpred.failed
