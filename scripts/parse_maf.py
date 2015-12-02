#! /usr/bin/env python
# encoding: utf-8


def parse_maf(maf_path):
	with open(maf_path) as f:
		a = f.readline()
		while(a):
			ref = f.readline().split()[6].replace('-', '')
			read = f.readline().split()[6].replace('-', '')
			lf = f.readline()
			a = f.readline()

			yield((ref, read))

import sys
print(sys.argv[1])
pairs = parse_maf(sys.argv[1])

for pair in pairs: print(pair[0], pair[1])
