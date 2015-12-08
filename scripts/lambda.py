#! /usr/bin/env python
# encoding: utf-8

if __name__ == '__main__':

	from sys import argv
	from math import exp

	m = float(argv[1])
	x = float(argv[2])
	p = 0.25 # const

	expected_score = p * m + 3 * p * x
	(l, high, low) = (1, 2, 0)

	while (high - low) > 0.001:
		sum = 4 * p * p * exp(l * m) + 12 * p * p * exp(l * x)
		if sum > 1.0:
			high = l
			l = (l + low) / 2.0
		else:
			low = l
			l = (l + high) / 2.0

	identity = 4.0 * p * p * exp(l * m)
	h = l * m * identity + l * x * (1.0 - identity)

	print("expected score: {0}".format(expected_score))
	print("lambda: {0}".format(l))
	print("H: {0}".format(h))
	print("identity: {0}".format(identity))

