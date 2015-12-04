#! /usr/bin/env python
# encoding: utf-8

from evaluate import *

def randbase():
	return(['A', 'C', 'G', 'T'][random.randint(0, 3)])

def gap_insert_modifier(seq, param):
	length = len(seq)
	insert_pos = random.randint(int(0.1 * length), int(0.6 * length))
	insert_seq = ''.join([randbase() for i in range(param)])
	tail_seq = ''.join([randbase() for i in range(100)])
	# print(seq, insert_pos, insert_seq, tail_seq, length, param)
	return(seq[:insert_pos] + insert_seq + seq[insert_pos:] + tail_seq)


def evaluate_gap(pbsim_path, ref_path, bin_path,
	ref_gap_length, read_gap_length, bandwidth, error_rate, length, count):

	# m = 2
	x = -4
	gi = -4
	# ge = -2
	return(evaluate_impl(pbsim_path, ref_path, bin_path,
		gi, x, bandwidth, error_rate, length, count,
		gap_insert_modifier, ref_gap_length, read_gap_length))



if __name__ == '__main__':

	import sys

	parse_and_run(sys.argv[1], sys.argv[2], evaluate_gap)

