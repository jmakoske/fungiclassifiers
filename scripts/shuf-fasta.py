#!/usr/bin/env python
# coding: utf-8

import sys
import random

with open(sys.argv[1]) as f:
    data = f.read().splitlines()
    indices = list(range(0, len(data), 2))
    random.shuffle(indices)
    for i in indices:
        print(data[i])
        print(data[i+1])
