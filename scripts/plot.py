#! /usr/bin/env python
# encoding: utf-8

from params import *
from util import *
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

dashes = [
	[8, 2],
	[2, 2],
	[6, 6],
	[6, 2, 2, 2],
	[6, 2, 2, 2, 2, 2],
	[4, 2, 4, 2, 2, 2],
	[4, 2, 2, 2],
	[8, 2, 2, 2]]

linewidth = 1.5
fontsize = 20

def dump_heatmap_data(in_list, first_label, first_legend, second_label, second_legend, out_file, count = count):

	with open(out_file, "w") as o:

		o.write('\t'.join([first_label, second_label, 'Recall']) + '\n')
		for (f, l) in zip(first_legend, in_list):
			for (s, e) in zip(second_legend, l):
				o.write('\t'.join(['"' + str(x) + '"' for x in [f, s, e * 100.0 / count]]) + '\n')

result_linear_2_x_y_2 = np.array(load_result('../results/linear_2_x_y_2_ddiag_table.txt'))
result_affine_2_x_y_1 = np.array(load_result('../results/affine_2_x_y_1_ddiag_table.txt'))
result_affine_2_x_y_2 = np.array(load_result('../results/affine_2_x_y_2_ddiag_table.txt'))
result_affine_2_x_y_3 = np.array(load_result('../results/affine_2_x_y_3_ddiag_table.txt'))

result_gige = np.array(load_result('../results/gige_id_ddiag_table.txt'))

# fix x
xs_fixed = [-x for x in xs]

# fix gi
# convert gap penalty model
# from g(k) = g_i + (k - 1) g_e
# to   g(k) = g_i + k g_e
gs_fixed_linear_2_x_y_2 = [-x for x in gs]
gs_fixed_affine_2_x_y_1 = [-x - 1 for x in gs]
gs_fixed_affine_2_x_y_2 = [-x - 2 for x in gs]
gs_fixed_affine_2_x_y_3 = [-x - 3 for x in gs]

# plot bw-id
dump_heatmap_data(result_affine_2_x_y_2[2, 2, :, :, 3].transpose(), 'Identity', error_rates, 'BW', bandwidths, 'bw_id.txt')

# plot id-len
dump_heatmap_data(result_affine_2_x_y_2[2, 2, 1, :, :], 'Identity', error_rates, 'Length', lengths, 'id_len_16.txt')

# plot id-len
dump_heatmap_data(result_affine_2_x_y_2[2, 2, 3, :, :], 'Identity', error_rates, 'Length', lengths, 'id_len_32.txt')

# plot x-gi
dump_heatmap_data(result_linear_2_x_y_2[:, :, 3, 4, 3], 'Ge', gs_fixed_linear_2_x_y_2, 'X', xs_fixed, 'x_ge_linear_2_x_y_2.txt')
dump_heatmap_data(result_affine_2_x_y_1[0:len(gs)-2, :, 3, 4, 3], 'Gi', gs_fixed_affine_2_x_y_1, 'X', xs_fixed, 'x_gi_affine_2_x_y_1.txt')
dump_heatmap_data(result_affine_2_x_y_2[1:len(gs)-1, :, 3, 4, 3], 'Gi', gs_fixed_affine_2_x_y_2, 'X', xs_fixed, 'x_gi_affine_2_x_y_2.txt')
dump_heatmap_data(result_affine_2_x_y_3[2:len(gs), :, 3, 4, 3], 'Gi', gs_fixed_affine_2_x_y_3, 'X', xs_fixed, 'x_gi_affine_2_x_y_3.txt')

# plot gige
dump_heatmap_data(result_gige[:, :, 0, 0, 0], 'Gi', map(abs, gis), 'Ge', map(abs, ges), 'gige_id_ddiag_extract.txt')

def plot_gap_bench(in_file, out_file):
	# clear figure
	plt.clf()

	# ref_gaps, read_gaps = 0, bw, id = 0.75, len = 10k
	gap_affine = np.array(load_result(in_file))[:, 0, :, 4, 1]

	for r, bw, d in zip(gap_affine.transpose(), bandwidths, dashes):
		y = [x / 10.0 for x in r.tolist()]
		x = range(len(y))
		# print(x, y)
		l, = plt.plot(x, y, color = 'black', linewidth = linewidth, linestyle = '-', label = str(bw))
		l.set_dashes(d)

	plt.xlim(0, 100)
	plt.xticks(fontsize = fontsize)
	plt.xlabel('Gap insert size (bp)', fontsize = fontsize)

	plt.ylim(0, 105)
	plt.yticks(fontsize = fontsize)
	plt.ylabel('Recall rate (%)', fontsize = fontsize)

	plt.title('Gap insert size - recall rate', fontsize = fontsize)
	plt.legend(title = 'BW', loc = 'best', fontsize = fontsize * 0.6)

	plt.savefig(out_file, dpi = 72)

# plot gap bench

plot_gap_bench('../results/result_gap_linear.txt', 'gap_linear.eps')
plot_gap_bench('../results/result_gap_affine.txt', 'gap_affine.eps')

def read_csv(filename):

	with open(filename) as f:
		a = []
		l = f.readline()
		while l:
			a.append(map(eval, l.strip().split('\t')))
			l = f.readline()

	return(a)


def plot_calc_time(in_file, out_file):
	
	lengths = [100, 150, 250, 350, 500, 650, 800, 1000, 1500, 2500, 3500, 5000, 6500, 8000, 10000]
	variants = ['BLAST', 'BLAST (SIMD)', 'Adaptive (32)', 'Wavefront']

	# clear figure
	plt.clf()

	# ref_gaps, read_gaps = 0, bw, id = 0.75, len = 10k
	data = np.array(read_csv(in_file))

	for r, bw, d in zip(data.transpose(), variants, dashes):
		ys = r.tolist()
		xs = lengths

		# filt out zeros
		xs_filt = [x for (x, y) in zip(xs, ys) if y is not 0]
		ys_filt = [y for (x, y) in zip(xs, ys) if y is not 0]

		# print(x, y)
		l, = plt.plot(xs_filt, ys_filt, color = 'black', linewidth = linewidth, linestyle = '-', label = str(bw))
		l.set_dashes(d)

	plt.xlim(90, 10500)
	plt.xticks(fontsize = fontsize)
	plt.xlabel('Length (bp)', fontsize = fontsize)
	plt.xscale('log')

	dl = [[x for x in xs if x is not 0] for xs in data.tolist()]
	ymin = 0.95 * min([min(x) for x in dl])
	ymax = 1.05 * max([max(x) for x in dl])
	plt.ylim(ymin, ymax)
	plt.yticks(fontsize = fontsize)
	plt.ylabel('Calc. time (us)', fontsize = fontsize)
	plt.yscale('log')

	plt.title('Calc. time - Sequence length', fontsize = fontsize)
	plt.legend(title = 'BW', loc = 'best', fontsize = fontsize * 0.6)

	plt.savefig(out_file, dpi = 72)

plot_calc_time('../results/bench.txt', 'bench.eps')
