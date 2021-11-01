#!/usr/bin/env python
# coding: utf-8

import sys

with open(sys.argv[1]) as f:
    for i_line, line in enumerate(f.readlines()):
        if i_line == 0:
            print(line, end="")
            continue
        parts = line.rstrip().split("\t")
        for j in range(len(parts), 8):
            parts.append("Unknown_"+parts[-1])
        print("\t".join(parts))
