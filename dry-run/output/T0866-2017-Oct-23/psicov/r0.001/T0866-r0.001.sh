#!/bin/bash
touch T0866-r0.001.running
echo "running T0866-r0.001 .."
date
/home/badri/psicov/psicov -o -r 0.001 T0866.aln > T0866-r0.001.psicov
if [ -s "T0866-r0.001.psicov" ]; then
   mv T0866-r0.001.running T0866-r0.001.done
   echo "T0866-r0.001 job done."
   date
   exit
fi
mv T0866-r0.001.running T0866-r0.001.failed
echo "psicov job T0866-r0.001 failed!"
date
