#! /usr/bin/env python
# encoding: utf-8

from evaluate import *

def evaluate_gige(pbsim_path, ref_path, bin_path, gi, ge, bandwidth, error_rate, length, count):

	# m = 2
	x = -2

	return(evaluate_impl2(pbsim_path, ref_path, bin_path,
		x, gi, ge, bandwidth, error_rate, length, count,
		default_modifier, 1, 0))


if __name__ == '__main__':

	import sys

	parse_and_run(sys.argv[1], sys.argv[2], evaluate_gige)

