#! /usr/bin/env python

import bench


def evaluate(pbsim_path, ref_path, algorithm, length, error_rate):

	# bench.pbsim(pbsim_path, ref_path, length, 20, error_rate)
	pairs = bench.parse_maf('./sd_0001.maf')

	tot = 0
	acc = [0, 0]
	acc_blast = [0, 0]
	acc_ddiag = [0, 0]
	fail = [0, 0]

	hist = [[0 for i in range(1, 2000)] for j in [0, 1]]
	print(hist)

	for pair in pairs:

		ref = pair[0] + ''.join(['A' for i in range(100)])
		read = pair[1] + ''.join(['T' for i in range(100)])

		print(ref, read)

		scores = bench.align(['./ddiag', './ddiag-40'], ['linear', 'affine'], ref, read)
		succ = [1 if score[0] == score[1] else 0 for score in scores]
		blast = [1 if score[0] > score[1] else 0 for score in scores]
		ddiag = [1 if score[0] < score[1] else 0 for score in scores]

		print(scores, succ, acc)
		acc = [sum(x) for x in zip(succ, acc)]
		acc_blast = [sum(x) for x in zip(blast, acc_blast)]
		acc_ddiag = [sum(x) for x in zip(ddiag, acc_ddiag)]

		hist[0][scores[0][0] - scores[0][1] + 1000] += 1
		hist[1][scores[1][0] - scores[1][1] + 1000] += 1

		tot = tot + 1
		if tot == 1000:
			break
		

	# print(succ, fail)
	print(acc)
	print(acc_blast)
	print(acc_ddiag)
	print(hist)
	return(acc)

r = evaluate(bench.pbsim_path, bench.ref_path, 'affine', 100, 0.75)
