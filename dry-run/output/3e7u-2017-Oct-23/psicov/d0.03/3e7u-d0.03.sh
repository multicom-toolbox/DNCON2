#!/bin/bash
touch 3e7u-d0.03.running
echo "running 3e7u-d0.03 .."
date
/home/badri/psicov/psicov -o -d 0.03 3e7u.aln > 3e7u-d0.03.psicov
if [ -s "3e7u-d0.03.psicov" ]; then
   mv 3e7u-d0.03.running 3e7u-d0.03.done
   echo "3e7u-d0.03 job done."
   date
   exit
fi
mv 3e7u-d0.03.running 3e7u-d0.03.failed
echo "psicov job 3e7u-d0.03 failed!"
date
