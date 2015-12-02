#! /usr/bin/env python
# encoding: utf-8

from util import *

if __name__ == '__main__':

	params = generate_params(xs, gs, bandwidths, error_rates, lengths, count)
	for param in params:
		print('\t'.join([str(l) for l in list(param)[2:]]))


