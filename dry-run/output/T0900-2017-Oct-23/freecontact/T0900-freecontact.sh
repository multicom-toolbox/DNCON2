#!/bin/bash
touch freecontact.running
echo "running freecontact .."
/usr/bin/freecontact < T0900.aln > T0900.freecontact.rr
if [ -s "T0900.freecontact.rr" ]; then
   mv freecontact.running freecontact.done
   echo "freecontact job done."
   exit
fi
echo "freecontact failed!"
mv freecontact.running freecontact.failed
