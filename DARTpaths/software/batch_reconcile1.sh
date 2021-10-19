#!/bin/bash
set -e

# Batch script to run reconcile1.sh on multiple proteins (for overnight).
echo 'Running reconcile1 on '$@

for prot in "$@"
do
    echo '----------Working on '${prot}'----------'
    cd ./${prot}
    bash ./reconcile1.sh ${prot}
    cd ..
done
