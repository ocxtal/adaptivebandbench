#! /usr/bin/env python
# encoding: utf-8


from util import *

if __name__ == '__main__':
	
	prefix = 'probe'
	length = 10
	count = 10000
	bandwidth = 32
	tail_length = 10

	pbsim(pbsim_path, ref_path, prefix,
		length,
		1.2 * count * length / ref_length,
		0.8)
	pairs = parse_maf('./{0}_0001.maf'.format(prefix))

	for pair in pairs:
		ref = default_modifier(pair[0], 0, length = tail_length)
		read = default_modifier(pair[1], 0, length = tail_length)

		result = align(
			['{0}/full-{1}'.format(bin_path, bandwidth),
			 '{0}/ddiag-{1}'.format(bin_path, bandwidth)],
			'affine', ref, read,
			2, -2, -2, -2,
			100)
		# print(result)

		if result[0][0] != result[1][0]:
			print('\t'.join([str(x) for x in [result[0][0], result[1][0], ref, read]]))

