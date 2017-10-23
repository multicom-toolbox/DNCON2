#!/bin/bash
touch ccmpred.running
echo "running ccmpred .."
/home/badri/CCMpred/bin/ccmpred -t 8 3e7u.aln 3e7u.ccmpred > ccmpred.log
if [ -s "3e7u.ccmpred" ]; then
   mv ccmpred.running ccmpred.done
   echo "ccmpred job done."
   exit
fi
echo "ccmpred failed!"
mv ccmpred.running ccmpred.failed
