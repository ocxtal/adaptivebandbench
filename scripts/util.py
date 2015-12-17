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

def ints(s):
	try: return(int(s))
	except ValueError: return(str(s))


def align(align_paths, alg, ref, read, m, x, gi, ge, xdrop_coef):
	print(align_paths)
	return([[ints(s) for s in subprocess.check_output([
		path, alg, ref, read, str(m), str(x), str(gi), str(ge),
		str((60 + xdrop_coef) * -gi)]).split()]
			for path in align_paths])

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

	results_linear = []
	results_affine = []
	tot = 0
	for pair in pairs:

		ref = modifier(pair[0], param1)
		read = modifier(pair[1], param2)

		results_linear.append(align(
			['{0}/full-{1}'.format(bin_path, bandwidth),
			 '{0}/blast-{1}'.format(bin_path, bandwidth),
			 '{0}/ddiag-{1}'.format(bin_path, bandwidth)],
			'linear', ref, read,
			2, x, gi, -2,		# m, x, gi, ge
			max(param1, param2)))

		results_affine.append(align(
			['{0}/full-{1}'.format(bin_path, bandwidth),
			 '{0}/blast-{1}'.format(bin_path, bandwidth),
			 '{0}/ddiag-{1}'.format(bin_path, bandwidth)],
			'affine', ref, read,
			2, x, gi, -2,		# m, x, gi, ge
			max(param1, param2)))
		
		tot = tot + 1
		if tot == count:
			break
	
	cleanup_pbsim(prefix)
	return([results_linear, results_affine])

def evaluate(pbsim_path, ref_path, bin_path, gi, x, bandwidth, error_rate, length, count):
	return(evaluate_impl(pbsim_path, ref_path, bin_path,
		gi, x, bandwidth, error_rate, length, count,
		default_modifier, 0, 1))

"""
def aggregate_hist():
	indices = [clip(score[1] - score[0] + hist_size/2, 0, hist_size-1) for score in scores]
	# print(indices, scores)
	for (i, index) in zip(range(len(indices)), indices):
		hist[i][index] = hist[i][index] + 1
	tot = 0
	hist = [[0 for i in range(hist_size)], [0 for j in range(hist_size)]]
"""

# utility
def clip(x, min, max):
	if x < min:
		return min
	elif x > max:
		return max
	else:
		return x

def num(s):
	try: return(int(s))
	except ValueError: return(float(s))


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
def slice_impl(arr, index):
	if index == -1:
		return(arr)
	else:
		return(arr[index])
def slice_array(arr, indices):

	def slice_array_intl(arr, indices, fn):
		if indices[0] == -1:
			return([fn(a, indices[1:]) for a in arr])
		else:
			return(fn(arr[indices[0]], indices[1:]))

	if len(indices) == 1:
		return(slice_array_intl(arr, indices, lambda x, y: x))
	else:
		return(slice_array_intl(arr, indices, slice_array))

def aggregate(input_file, output_file, params_list):

	dimensions = [len(p) for p in params_list]
	results_linear = array(dimensions)
	results_affine = array(dimensions)

	import sys

	lines = []
	with open(input_file) as r: lines = r.readlines()

	for line in lines:
		p = [eval(s) for s in line.split('\t')]
		# print(p)
		indices = [l.index(q) for (l, q) in zip(params_list, p)]
		# print(indices)

		r = p[len(params_list) + 1:]

		set_array(results_linear, indices, r[0])
		set_array(results_affine, indices, r[1])

	with open(output_file, "w") as w:
		w.write(str(results_linear) + '\n')
		w.write(str(results_affine) + '\n')

def load_result_impl(filename):

	# from numpy import array

	linear = []
	affine = []
	with open(filename) as f:
		linear = eval(f.readline())
		affine = eval(f.readline())

	return(linear, affine)

def load_result(filename):
	return(load_result_impl(filename))

# print utilities
def tab_formatter(arr):
	return('\n'.join(['\t'.join(e) for e in arr]))

def extract_labels(params_list, indices):
	label = [params_list[0]] if indices[0] == -1 else []
	if len(indices) == 1:
		return(label)
	else:
		return(label + extract_labels(params_list[1:], indices[1:]))

def score_identity(elems, comp_pair = [0, 1]):
	try: return(sum([elem[comp_pair[0]][0] == elem[comp_pair[1]][0] for elem in elems]))
	except TypeError: return('-')

def score_difference(elems, bin_size = 30, comp_pair = [0, 1]):
	bin = [0 for i in range(bin_size)]
	try:
		for elem in elems:
			bin[clip(elem[comp_pair[0]][0] - elem[comp_pair[1]][0] + bin_size/2, 0, bin_size-1)] += 1
		return([str(e) for e in bin])
	except TypeError: return([str(e) for e in bin])
	# return([clip(elem[0][0] - e[0] + bin_size/2, 0, bin_size-1) for e in elem[1:]])

def table(arr, fn, indices, params_list = params_list, formatter = tab_formatter, comp_pair = [0, 1]):
	# indices must have two -1's
	labels = extract_labels(params_list, indices)
	table_elems = [['$'] + [str(e) for e in labels[1]]]
	table_elems += [[str(l)] + [str(fn(e, comp_pair = comp_pair)) for e in es]
		for (l, es) in zip(labels[0], slice_array(arr, indices))]
	print(formatter(table_elems))

def hist(arr, fn, indices, params_list = params_list, formatter = tab_formatter, comp_pair = [0, 1], bin_size = 30):
	# indices must have just one -1 in it
	bins = [['%'] + [str(i) for i in range(bin_size)]]
	labels = extract_labels(params_list, indices)
	for (label, result) in zip(labels[0], slice_array(arr, indices)):
		bins.append([str(label)] + fn(result, bin_size = bin_size, comp_pair = comp_pair))
	print(bins)
	print(formatter(bins))
	# print('\t'.join([str(e) for e in bin]))




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

