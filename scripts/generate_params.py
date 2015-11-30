#! /usr/bin/env python
# encoding: utf-8

from util import *

lengths = [100, 300, 1000, 3000, 10000]
error_rates = [0.65, 0.75, 0.85, 0.95]
bandwidths = [16, 32, 48, 64]
xs = [-i for i in range(2, 7)]		# -1..-3
gs = [-i for i in range(2, 11)]		# -1..-5

if __name__ == '__main__':

	params = generate_params(bandwidths, xs, gs, lengths, error_rates)
	for param in params:
		print('\t'.join([str(l) for l in list(param)[2:]]))


