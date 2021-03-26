#!/bin/bash
echo -e "#Sequence ID\tkingdom\tphylum\tclass\torder\tfamily\tgenus\tspecies"
cat $1 | cut -f1,3 | tr " ," "\t\t"
