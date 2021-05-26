#!/usr/bin/env python
# coding: utf-8

import sys
import pandas as pd
from collections import defaultdict
from random import seed, sample, shuffle

if len(sys.argv)<6:
    print('USAGE: {} fastafile taxonomyfile level '
          'trainfile testfile'.format(sys.argv[0]))
    sys.exit()

df = pd.read_csv(sys.argv[2], sep="\t", index_col="#Sequence ID")
level = sys.argv[3]
df['isnull'] = df[level].isnull()

categories = defaultdict(list)
sequences = list()

with open(sys.argv[1]) as f:
    for i_line, line in enumerate(f.readlines()):
        if line.startswith(">"):
            seq = line[1:]
            seq = seq.rsplit()[0]
            cat = df.loc[seq][level]
            if not df.loc[seq]['isnull']:
                categories[cat].append(seq)
                sequences.append(seq)
n_seq = len(sequences)
n_test = int(n_seq*0.2)
print(n_seq, n_test)

test_sequences = dict()
seed(42)
while True:
    sampled_cat = sample(categories.keys(), 1)[0]
    sampled_seq = sample(categories[sampled_cat], 1)[0]
    if sampled_seq not in test_sequences:
        test_sequences[sampled_seq] = True
        if len(test_sequences) == n_test:
            break

with open(sys.argv[5], 'w') as f:
    for s in test_sequences.keys():
        f.write(s + '\n')

with open(sys.argv[4], 'w') as f:
    for s in sequences:
        if s not in test_sequences:
            f.write(s + '\n')
