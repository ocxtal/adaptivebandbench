#! /usr/bin/env python
# encoding: utf-8

from util import *

# utility
def num(s):
	try: return int(s)
	except ValueError: return float(s)

def parse_and_run(input_file, output_file, evaluate_function):

	lines = []
	with open(input_file, "r") as r: lines = r.readlines()

	with open(output_file, "w") as w:
		for line in lines:
			p = [num(s) for s in line.split()]
			# print(p)
			res = evaluate_function(pbsim_path, ref_path, bin_path, *p)
			w.write('{}\n'.format('\t'.join([str(r) for r in p+res])))


if __name__ == '__main__':

	import sys

	# pbsim_path = '/Users/suzukihajime/docs/src/dl/pbsim-1.0.3/src/'
	# ref_path = '/Users/suzukihajime/docs/lab/oni/work/NC_000913.fna'
	# bin_path = '..'

	parse_and_run(sys.argv[1], sys.argv[2], evaluate)

