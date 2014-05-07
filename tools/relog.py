#!/usr/bin/env python

import sys

headers = None

try:
    log = open(sys.argv[1])
except:
    log = sys.stdin

i = 0
for l in log:
    if l[0] != '#':
        if not headers:
            headers = l.split()
            factor_index = headers.index(sys.argv[2])
            r1,r2,r3 = headers.index(sys.argv[3]), headers.index(sys.argv[4]), headers.index(sys.argv[5])
            print '\t'.join([state] + sys.argv[3:6])
        else:
            l = l.split()
            factor = float(l[factor_index])
            print '\t'.join([str(i), str(float(l[r1]) * factor), str(float(l[r2]) * factor), str(float(l[r3]) * factor)])
            i += 1

if log is not sys.stdin:
    log.close()
