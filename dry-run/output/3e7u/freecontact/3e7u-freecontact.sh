#!/bin/bash
touch freecontact.running
echo "running freecontact .."
/usr/bin/freecontact < 3e7u.aln > 3e7u.freecontact.rr
if [ -s "3e7u.freecontact.rr" ]; then
   mv freecontact.running freecontact.done
   echo "freecontact job done."
   exit
fi
echo "freecontact failed!"
mv freecontact.running freecontact.failed
