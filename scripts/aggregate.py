#! /usr/bin/env python
# encoding: utf-8

from util import *

# utility
def num(s):
	try: return int(s)
	except ValueError: return float(s)


def array(size):
	if len(size) == 1:
		return([0 for i in range(size[0])])
	else:
		return([array(size[1:]) for i in range(size[0])])
def set_array(arr, indices, val):
	if len(indices) == 1:
		arr[indices[0]] = val
	else:
		set_array(arr[indices[0]], indices[1:], val)
def get_array(arr, indices):
	if len(indices) == 1:
		return(arr[indices[0]])
	else:
		return(get_array(arr[indices[0]], indices[1:]))


if __name__ == '__main__':

	results_linear = array(dimensions)
	results_affine = array(dimensions)

	import sys

	lines = []
	with open(sys.argv[1]) as r: lines = r.readlines()

	for line in lines:
		p = [num(s) for s in line.split()]
		# print(p)
		indices = [l.index(q) for (l, q) in zip(params_list, p)]
		# print(indices)

		r = p[len(params_list) + 1:]

		set_array(results_linear, indices, r[0])
		set_array(results_affine, indices, r[1])

	with open(sys.argv[2], "w") as w:
		w.write(str(results_linear) + '\n')
		w.write(str(results_affine) + '\n')
