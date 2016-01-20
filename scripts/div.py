#! /usr/bin/env python
# encoding: utf-8

from params import *
from util import *

def div(input, out1, out2, params_list = params_list):

	print(input, out1, out2)

	f = open(input, 'r')
	o1 = open(out1, 'w')
	o2 = open(out2, 'w')

	line = f.readline()
	while line:
		p = line.split('\t')
		h = p[:len(params_list)]
		r = p[len(params_list) + 1:]
		
		o1.write('\t'.join(h + [r[0].split('\n')[0]]) + '\n')
		o2.write('\t'.join(h + [r[1].split('\n')[0]]) + '\n')

		line = f.readline()


if __name__ == '__main__':
	import sys

	div(sys.argv[1], sys.argv[2], sys.argv[3])
	
