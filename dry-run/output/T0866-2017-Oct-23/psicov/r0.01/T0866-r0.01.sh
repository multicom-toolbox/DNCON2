#!/bin/bash
touch T0866-r0.01.running
echo "running T0866-r0.01 .."
date
/home/badri/psicov/psicov -o -r 0.01 T0866.aln > T0866-r0.01.psicov
if [ -s "T0866-r0.01.psicov" ]; then
   mv T0866-r0.01.running T0866-r0.01.done
   echo "T0866-r0.01 job done."
   date
   exit
fi
mv T0866-r0.01.running T0866-r0.01.failed
echo "psicov job T0866-r0.01 failed!"
date
