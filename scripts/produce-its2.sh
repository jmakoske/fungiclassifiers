#!/bin/bash

if [ "$#" -ne 4 ]; then
    echo "USAGE: $0 prefix kmers level-as-number level-as-word"
    exit
fi

echo "Input files: $1-$3-$2.valid.txt, $1-$3-$2.probsnew.valid.txt"
echo "Output file: $1-$4-$2.results.txt"

nl -nln $1-$3-$2.valid.txt > a
cat $1-$3-$2.probsnew.valid.txt | cut -d, -f4- | tr , "\t" > b
echo -e "Index\tID\t$4\tprob_$4" > $1-$4-$2.results.txt 
paste a b >> $1-$4-$2.results.txt
rm a b
