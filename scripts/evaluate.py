#! /usr/bin/env python
# encoding: utf-8

from util import *

# utility
def num(s):
	try: return int(s)
	except ValueError: return float(s)

if __name__ == '__main__':
	
	import sys
	# pbsim_path = '/Users/suzukihajime/docs/src/dl/pbsim-1.0.3/src/'
	# ref_path = '/Users/suzukihajime/docs/lab/oni/work/NC_000913.fna'

	lines = []
	with open(sys.argv[1], "r") as r: lines = r.readlines()

	with open(sys.argv[2], "w") as w:
		for line in lines:
			p = [num(s) for s in line.split()]
			res = evaluate(pbsim_path, ref_path, *p)
			w.write('{}\n'.format('\t'.join([str(r) for r in p+res])))
