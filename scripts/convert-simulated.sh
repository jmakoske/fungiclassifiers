#!/bin/bash
# Run as:
# scripts/convert-simulated.sh ../TaxonomicClassifier/Simulated_data/Ale_simulations/Scenario1/all.fasta > data/simulated-scenario1.classification
echo -e "#Sequence ID\tkingdom\tphylum\tclass\torder"
grep "^>" $1 | tr -d ">" | tr \| "\t"

