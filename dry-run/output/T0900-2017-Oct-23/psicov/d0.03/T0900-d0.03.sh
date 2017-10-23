#!/bin/bash
touch T0900-d0.03.running
echo "running T0900-d0.03 .."
date
/home/badri/psicov/psicov -o -d 0.03 T0900.aln > T0900-d0.03.psicov
if [ -s "T0900-d0.03.psicov" ]; then
   mv T0900-d0.03.running T0900-d0.03.done
   echo "T0900-d0.03 job done."
   date
   exit
fi
mv T0900-d0.03.running T0900-d0.03.failed
echo "psicov job T0900-d0.03 failed!"
date
