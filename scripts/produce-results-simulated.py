#!/usr/bin/env python3
# coding: utf-8

import sys

if len(sys.argv)<3:
    print('USAGE: {} template outputfile'.format(sys.argv[0]))
    sys.exit()


levels = ['1', '2', '3', '4']
outputtitle = 'Index\tID\tSp1\tSp2\tSp3\tSp4\tprob_Sp1\tprob_Sp2\tprob_Sp3\tprob_Sp4\n'

datadict = dict()
for lvl in levels:
    prefix = sys.argv[1].format(lvl,lvl)
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

outfile = sys.argv[2]
with open(outfile, 'w') as f:
    f.write(outputtitle)
    i = 1
    for k, v in datadict.items():
        f.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n"
                .format(i, k, v["res-1"], v["res-2"], v["res-3"], v["res-4"],
                        v["prob-1"], v["prob-2"], v["prob-3"], v["prob-4"]))
        i = i+1
print('Wrote', outfile)
