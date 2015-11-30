#! /usr/bin/env python
# encoding: utf-8

import subprocess

pbsim_path = '/home/suzukihajime/src/pbsim-1.0.3/src/'
ref_path = '/home/suzukihajime/oni/work/resource/NC_000913.fna'
ref_length = 4700000
align_path = './a.out'
# ./pbsim --data-type CLR --length-min 900 --length-max 1100 --accuracy-min 0.84 --accuracy-max 0.86 --model_qc ../data/model_qc_clr --length-mean 1000 --length-sd 100 --accuracy-mean 0.85 --accuracy-sd 0.01 ~/docs/oni/work/NC_000913.fna

def pbsim(pbsim_path, ref_path, prefix, length, depth, accuracy):

	length_sd = length * 0.05
	accuracy_sd = 0.01

	ret = subprocess.call([
		'/'.join([pbsim_path, 'pbsim']),
		ref_path,
		'--prefix={}'.format(prefix),
		'--data-type=CLR',
		'--depth={}'.format(depth),
		'--length-mean={}'.format(length),
		'--length-min={}'.format(length - length_sd),
		'--length-max={}'.format(length + length_sd),
		'--length-sd={}'.format(length_sd),
		'--accuracy-mean={}'.format(accuracy),
		'--accuracy-min={}'.format(accuracy - accuracy_sd),
		'--accuracy-max={}'.format(accuracy + accuracy_sd),
		'--accuracy-sd={}'.format(accuracy_sd),
		'--model_qc={}'.format('/'.join([pbsim_path, '../data/model_qc_clr']))])
	return(ret)

def cleanup_pbsim(prefix):
	subprocess.call(['rm', './{}_0001.fastq'.format(prefix)])
	subprocess.call(['rm', './{}_0001.maf'.format(prefix)])
	subprocess.call(['rm', './{}_0001.ref'.format(prefix)])

def parse_maf(maf_path):
	with open(maf_path) as f:
		a = f.readline()
		while(a):
			ref = f.readline().split()[6].replace('-', '')
			read = f.readline().split()[6].replace('-', '')
			lf = f.readline()
			a = f.readline()

			yield((ref, read))

def align(align_paths, algorithms, ref, read, m, x, gi, ge):
	return([[int(subprocess.check_output([
		path, alg, ref, read, str(m), str(x), str(gi), str(ge), str(20 * -x)]).split()[0])
			for path in align_paths]
			for alg in algorithms])


def evaluate(pbsim_path, ref_path, bandwidth, m, x, gi, ge, length, error_rate, count):

	prefix = '{}_{}_{}_{}_{}_{}'.format(m, x, gi, ge, length, error_rate)
	pbsim(pbsim_path, ref_path, prefix,
		length,
		2 * count * length / ref_length,		# depth
		error_rate)
	pairs = parse_maf('./{}_0001.maf'.format(prefix))

	tot = 0
	acc = [0, 0]
	# fail = [0, 0]

	for pair in pairs:

		ref = pair[0] + ''.join(['A' for i in range(100)])
		read = pair[1] + ''.join(['T' for i in range(100)])

		# print(ref, read)

		scores = align(
			['./blast-{}'.format(bandwidth), './ddiag-{}'.format(bandwidth)],
			['linear', 'affine'], ref, read, m, x, gi, ge)
		succ = [1 if score[0] == score[1] else 0 for score in scores]
		# print(scores, succ, acc)
		acc = [sum(x) for x in zip(succ, acc)]
		tot = tot + 1
		if tot == count:
			break
	
	cleanup_pbsim(prefix)

	print(bandwidth, m, x, gi, ge, length, error_rate, acc)
	return(acc)

def apply(argv): return(argv[0](*argv[1:]))

if __name__ == '__main__':

	lengths = [100, 300, 1000, 3000, 10000]
	error_rates = [0.65, 0.75, 0.85, 0.95]
	bandwidths = [16, 32, 48, 64]
	xs = [-i for i in range(2, 7)]		# -1..-3
	gs = [-i for i in range(2, 11)]		# -1..-5

	from multiprocessing import Pool
	func_args = [(evaluate, pbsim_path, ref_path, b, 2, x, g, -2, l, e, 100)
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

	results = nest(results,
		[len(xs), len(gs), len(bandwidths), len(error_rates), len(lengths)])

	# print(results)

	from numpy import *
	res_arr = array(results)

	res_linear = res_arr[
		0:len(xs),
		0:len(gs),
		0:len(bandwidths),
		0:len(error_rates),
		0:len(lengths),
		0]
	res_affine = res_arr[
		0:len(xs),
		0:len(gs),
		0:len(bandwidths),
		0:len(error_rates),
		0:len(lengths),
		1]

	print(res_linear.tolist())
	print(res_affine.tolist())


