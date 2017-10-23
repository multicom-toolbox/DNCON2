#!/bin/bash
touch freecontact.running
echo "running freecontact .."
/usr/bin/freecontact < T0866.aln > T0866.freecontact.rr
if [ -s "T0866.freecontact.rr" ]; then
   mv freecontact.running freecontact.done
   echo "freecontact job done."
   exit
fi
echo "freecontact failed!"
mv freecontact.running freecontact.failed
