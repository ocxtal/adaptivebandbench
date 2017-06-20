#! /usr/bin/env python
# encoding: utf-8

from util import *

if __name__ == '__main__':

	import sys

	def naive_blast(r):
		return(score_identity(r, comp_pair = [0, 1]))
	def naive_ddiag(r):
		return(score_identity(r, comp_pair = [0, 2]))
	# aggregate(sys.argv[1], sys.argv[2], params_list = params_list, aggregator = naive_blast)
	aggregate(sys.argv[1], sys.argv[3], params_list = params_list, aggregator = naive_ddiag)
