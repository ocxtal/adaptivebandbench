#! /usr/bin/env python
# encoding: utf-8

import subprocess
import random
from params import *


# initalize rand seed
random.seed(None)

# pbsim
def pbsim(pbsim_path, ref_path, prefix, length, depth, accuracy):

	length_sd = length * 0.05
	accuracy_sd = 0.01

	ret = subprocess.call([
		'/'.join([pbsim_path, 'pbsim']),
		ref_path,
		'--prefix={0}'.format(prefix),
		'--data-type=CLR',
		'--depth={0}'.format(depth),
		'--length-mean={0}'.format(length),
		'--length-min={0}'.format(length - length_sd),
		'--length-max={0}'.format(length + length_sd),
		'--length-sd={0}'.format(length_sd),
		'--accuracy-mean={0}'.format(accuracy),
		'--accuracy-min={0}'.format(accuracy - accuracy_sd),
		'--accuracy-max={0}'.format(accuracy + accuracy_sd),
		'--accuracy-sd={0}'.format(accuracy_sd),
		'--model_qc={0}'.format('/'.join([pbsim_path, '../data/model_qc_clr']))])
	return(ret)

def cleanup_pbsim(prefix):
	subprocess.call(['rm', './{0}_0001.fastq'.format(prefix)])
	subprocess.call(['rm', './{0}_0001.maf'.format(prefix)])
	subprocess.call(['rm', './{0}_0001.ref'.format(prefix)])

def parse_maf(maf_path):
	with open(maf_path) as f:
		a = f.readline()
		while(a):
			ref = f.readline().split()[6].replace('-', '')
			read = f.readline().split()[6].replace('-', '')
			lf = f.readline()
			a = f.readline()

			yield((ref, read))

def align(align_paths, algorithms, ref, read, m, x, gi, ge, xdrop_coef):
	return([[int(subprocess.check_output([
		path, alg, ref, read, str(m), str(x), str(gi), str(ge),
		str((60 + xdrop_coef) * -gi)]).split()[0])
			for path in align_paths]
			for alg in algorithms])

# parameter generation
def generate_params(params_list, count):
	from itertools import product
	params = [[pbsim_path, ref_path] + list(p) + [count]
		for p in list(product(*params_list))]
	return(params)


def randbase():
	return(['A', 'C', 'G', 'T'][random.randint(0, 3)])

def default_modifier(seq, param):
	return(seq + ''.join([randbase() for i in range(100)]))


def evaluate_impl(pbsim_path, ref_path, bin_path,
	gi, x, bandwidth, error_rate, length, count,
	modifier, param1, param2):

	prefix = 's_{0}_{1}_{2}_{3}_{4}_{5}_{6}'.format(
		-x, -gi, bandwidth, error_rate, length, param1, param2)
	pbsim(pbsim_path, ref_path, prefix,
		length,
		1.2 * count * length / ref_length,		# depth
		error_rate)
	pairs = parse_maf('./{0}_0001.maf'.format(prefix))

	tot = 0
	acc = [0, 0]
	# fail = [0, 0]

	for pair in pairs:

		ref = modifier(pair[0], param1)
		read = modifier(pair[1], param2)

		# print(ref, read)

		scores = align(
			['{0}/blast-{1}'.format(bin_path, bandwidth), '{0}/ddiag-{1}'.format(bin_path, bandwidth)],
			['linear', 'affine'], ref, read,
			2, x, gi, -2,		# m, x, gi, ge
			max(param1, param2))
		succ = [1 if score[0] == score[1] else 0 for score in scores]

		acc = [sum(i) for i in zip(succ, acc)]
		tot = tot + 1

		# print(scores, succ, acc)
		if tot == count:
			break
	
	cleanup_pbsim(prefix)

	# print(bandwidth, m, x, gi, ge, error_rate, length, acc)
	return(acc)


def evaluate(pbsim_path, ref_path, bin_path, gi, x, bandwidth, error_rate, length, count):
	return(evaluate_impl(pbsim_path, ref_path, bin_path,
		gi, x, bandwidth, error_rate, length, count,
		default_modifier, 0, 1))



def load_results(filename):
	
	from numpy import array

	linear = []
	affine = []
	with open(filename) as f:
		linear = eval(f.readline())
		affine = eval(f.readline())

	return(array(linear), array(affine))

"""
def apply(argv): return(argv[0](*argv[1:]))

if __name__ == '__main__':

	# lengths = [100, 300, 1000, 3000, 10000]
	# error_rates = [0.65, 0.75, 0.85, 0.95]
	# bandwidths = [16, 32, 48, 64]
	# xs = [-i for i in range(2, 7)]		# -1..-3
	# gs = [-i for i in range(2, 11)]		# -1..-5

	from multiprocessing import Pool
	func_args = [(evaluate, pbsim_path, ref_path, '.', 2, x, g, -2, b, e, l, 100)
		for x in xs for g in gs
		for b in bandwidths
		for e in error_rates for l in lengths]

	p = Pool(12)
	results = p.map(apply, func_args)
	# results = func_args

	# print(results)

	# unflatten list
	from itertools import islice
	def nest(flat,levels): return _nest(flat,levels).__next__()
	def _nest(flat,levels):
		if levels:
			it = _nest(flat,levels[1:])
			while 1: yield list(islice(it,levels[0]))
		else:
			for d in flat: yield d

	results = nest(results, dimensions)

	# print(results)

	from numpy import *
	res_arr = array(results)

	res_linear = res_arr[
		0:dimensions[0], # xs
		0:dimensions[1], # gs
		0:dimensions[2], # bandwidths
		0:dimensions[3], # error_rates
		0:dimensions[4], # lengths
		0]
	res_affine = res_arr[
		0:dimensions[0], # xs
		0:dimensions[1], # gs
		0:dimensions[2], # bandwidths
		0:dimensions[3], # error_rates
		0:dimensions[4], # lengths
		1]

	print(res_linear.tolist())
	print(res_affine.tolist())
"""

