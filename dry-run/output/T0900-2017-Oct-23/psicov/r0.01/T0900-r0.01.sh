#!/bin/bash
touch T0900-r0.01.running
echo "running T0900-r0.01 .."
date
/home/badri/psicov/psicov -o -r 0.01 T0900.aln > T0900-r0.01.psicov
if [ -s "T0900-r0.01.psicov" ]; then
   mv T0900-r0.01.running T0900-r0.01.done
   echo "T0900-r0.01 job done."
   date
   exit
fi
mv T0900-r0.01.running T0900-r0.01.failed
echo "psicov job T0900-r0.01 failed!"
date
