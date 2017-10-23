#!/bin/bash
touch T0866-d0.03.running
echo "running T0866-d0.03 .."
date
/home/badri/psicov/psicov -o -d 0.03 T0866.aln > T0866-d0.03.psicov
if [ -s "T0866-d0.03.psicov" ]; then
   mv T0866-d0.03.running T0866-d0.03.done
   echo "T0866-d0.03 job done."
   date
   exit
fi
mv T0866-d0.03.running T0866-d0.03.failed
echo "psicov job T0866-d0.03 failed!"
date
