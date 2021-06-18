#!/usr/bin/env python3
# coding: utf-8

import sys

if len(sys.argv)<3:
    print('USAGE: {} template its2|finbol'.format(sys.argv[0]))
    sys.exit()

if sys.argv[2] == "its2":
    levels = ['3', '4', '5', '6']
    outputtitle = 'Index\tID\tClass\tOrder\tFamily\tGenus\tSpecies\tprob_Class\tprob_Order\tprob_Family\tprob_Genus\tprob_Species\n'
elif sys.argv[2] == "finbol":
    levels = ['1', '2', '3', '4', '5', '6', '7']
    outputtitle = 'Index\tID\tClass\tOrder\tFamily\tSubfamily\tTribe\tGenus\tSpecies\tprob_Class\tprob_Order\tprob_Family\tprob_Subfamily\tprob_Tribe\tprob_Genus\tprob_Species\n'
else:
    print('Unsupported mode:', sys.argv[2])
    sys.exit()

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
    f.write(outputtitle)
    i = 1
    for k, v in datadict.items():
        if sys.argv[2] == "its2":
            for l in [3,4,5,6]:
                if "res-{}".format(l) not in v:
                    v["res-{}".format(l)] = "MISSING"
                if "prob-{}".format(l) not in v:
                    v["prob-{}".format(l)] = "0.0"
        f.write("{}\t{}".format(i, k))
        if sys.argv[2] == "finbol":
            f.write("\t{}\t{}".format(v["res-1"], v["res-2"]))
        f.write("\t{}\t{}\t{}\t{}".format(v["res-3"], v["res-4"],
                                          v["res-5"], v["res-6"]))
        if sys.argv[2] == "finbol":
            f.write("\t{}\t{}\t{}".format(v["res-7"], v["prob-1"], v["prob-2"]))
        f.write("\t{}\t{}\t{}\t{}".format(v["prob-3"], v["prob-4"],
                                          v["prob-5"], v["prob-6"]))
        if sys.argv[2] == "finbol":
            f.write("\t{}".format(v["prob-7"]))
        f.write("\n")
        i = i+1
print('Wrote', outfile)

