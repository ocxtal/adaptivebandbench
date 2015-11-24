#! /usr/bin/env python
# encoding: utf-8

import subprocess

pbsim_path = '/home/suzukihajime/src/pbsim-1.0.3/src/'
ref_path = '/home/suzukihajime/oni/work/resource/NC_000913.fna'
align_path = './a.out'
# ./pbsim --data-type CLR --length-min 900 --length-max 1100 --accuracy-min 0.84 --accuracy-max 0.86 --model_qc ../data/model_qc_clr --length-mean 1000 --length-sd 100 --accuracy-mean 0.85 --accuracy-sd 0.01 ~/docs/oni/work/NC_000913.fna

def pbsim(pbsim_path, ref_path, length, depth, accuracy):

	length_sd = length * 0.05
	accuracy_sd = 0.01

	ret = subprocess.call([
		'/'.join([pbsim_path, 'pbsim']),
		ref_path,
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

def parse_maf(maf_path):
	with open(maf_path) as f:
		a = f.readline()
		while(a):
			ref = f.readline().split()[6].replace('-', '')
			read = f.readline().split()[6].replace('-', '')
			lf = f.readline()
			a = f.readline()

			yield((ref, read))

def align(align_paths, algorithms, ref, read):
	return([[int(subprocess.check_output([
		path, alg, ref, read, '2', '-3', '-5', '-1', '30']).split()[0])
			for path in align_paths]
			for alg in algorithms])


def evaluate(pbsim_path, ref_path, algorithm, length, error_rate):

	pbsim(pbsim_path, ref_path, length, 20, error_rate)
	pairs = parse_maf('./sd_0001.maf')

	tot = 0
	acc = [0, 0]
	fail = [0, 0]

	for pair in pairs:

		ref = pair[0] + ''.join(['A' for i in range(100)])
		read = pair[1] + ''.join(['T' for i in range(100)])
		# print(ref, read)
		scores = align(['./blast', './ddiag'], ['affine', 'linear'], ref, read)
		succ = [1 if score[0] == score[1] else 0 for score in scores]
		# print(scores, succ, acc)
		acc = [sum(x) for x in zip(succ, acc)]
		tot = tot + 1
		if tot == 1000:
			break
		

	# print(succ, fail)
	return(acc)

if __name__ == '__main__':

	lengths = [100, 200, 500, 1000, 2000, 5000, 10000]
	error_rates = [0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95]

	arg_list = [(len, err) for len in lengths for err in error_rates]

	for a in arg_list:
		r = evaluate(pbsim_path, ref_path, 'affine', a[0], a[1])
		print("{}, {}, {}, {}".format(a[0], a[1], r[0], r[1]))


