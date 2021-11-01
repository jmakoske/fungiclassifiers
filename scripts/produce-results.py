#!/usr/bin/env python3
# coding: utf-8

import sys

if len(sys.argv)<2:
    print('USAGE: {} template'.format(sys.argv[0]))
    sys.exit()


levels = ['1', '2', '3', '4', '5', '6', '7']
outputtitle = 'Index\tID\tClass\tOrder\tFamily\tSubfamily\tTribe\tGenus\tSpecies\tprob_Class\tprob_Order\tprob_Family\tprob_Subfamily\tprob_Tribe\tprob_Genus\tprob_Species\n'

datadict = dict()
for lvl in levels:
    prefix = sys.argv[1].format(lvl)
    seqfile = prefix+".valid.txt"
    print('Processing', seqfile)
    seqs = list()
    with open(seqfile) as f:
        for i_line, line in enumerate(f.readlines()):
            line = line.rsplit()[0]
            seqs.append(line)
    resfile = prefix+".probsnew.valid.txt"
    print('Processing', resfile)
    with open(resfile) as f:
        for i_line, line in enumerate(f.readlines()):
            line = line.rsplit()[0]
            parts = line.split(",")
            if seqs[i_line] not in datadict:
                datadict[seqs[i_line]] = {'res-'+lvl: parts[3],
                                          'prob-'+lvl: parts[4]}
            else:
                datadict[seqs[i_line]]['res-'+lvl]=parts[3]
                datadict[seqs[i_line]]['prob-'+lvl]=parts[4]

#print(datadict)
print(len(datadict), 'sequences')

outprefix = sys.argv[1].format("all")
outfile = outprefix+".results.txt"
with open(outfile, 'w') as f:
    i = 1
    for k, v in datadict.items():
        f.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n"
                .format(i, k, v["res-1"], v["res-2"], v["res-3"], v["res-4"],
                        v["res-5"], v["res-6"], v["res-7"],
                        v["prob-1"], v["prob-2"], v["prob-3"], v["prob-4"],
                        v["prob-5"], v["prob-6"], v["prob-7"]))
        i = i+1
print('Wrote', outfile)
