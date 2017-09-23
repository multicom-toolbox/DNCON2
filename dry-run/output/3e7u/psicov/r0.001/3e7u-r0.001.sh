#!/bin/bash
touch 3e7u-r0.001.running
echo "running 3e7u-r0.001 .."
date
/home/badri/psicov/psicov -o -r 0.001 3e7u.aln > 3e7u-r0.001.psicov
if [ -s "3e7u-r0.001.psicov" ]; then
   mv 3e7u-r0.001.running 3e7u-r0.001.done
   echo "3e7u-r0.001 job done."
   date
   exit
fi
mv 3e7u-r0.001.running 3e7u-r0.001.failed
echo "psicov job 3e7u-r0.001 failed!"
date
