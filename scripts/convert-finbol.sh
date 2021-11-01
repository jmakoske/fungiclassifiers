#!/bin/bash
echo -e "#Sequence ID\tkingdom\tphylum\tclass\torder\tfamily\tgenus\tspecies"
grep "^>" $1 | tr -d ">" | tr \| "\t"
grep "^>" $2 | tr -d ">" | tr \| "\t"
