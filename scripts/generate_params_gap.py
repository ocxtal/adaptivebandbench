#! /usr/bin/env python
# encoding: utf-8

from generate_params import *

if __name__ == '__main__':

	params = generate_params(params_list_gap, count)
	for param in params:
		print('\t'.join([str(l) for l in list(param)[2:]]))
