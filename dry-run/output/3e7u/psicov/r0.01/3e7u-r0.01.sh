#!/bin/bash
touch 3e7u-r0.01.running
echo "running 3e7u-r0.01 .."
date
/home/badri/psicov/psicov -o -r 0.01 3e7u.aln > 3e7u-r0.01.psicov
if [ -s "3e7u-r0.01.psicov" ]; then
   mv 3e7u-r0.01.running 3e7u-r0.01.done
   echo "3e7u-r0.01 job done."
   date
   exit
fi
mv 3e7u-r0.01.running 3e7u-r0.01.failed
echo "psicov job 3e7u-r0.01 failed!"
date
