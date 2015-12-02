#! /usr/bin/env python
# encoding: utf-8

pbsim_path = '/home/suzukihajime/src/pbsim-1.0.3/src/'
ref_path = '/home/suzukihajime/oni/work/resource/NC_000913.fna'
ref_length = 4700000
bin_path = '../bin'
# align_path = './a.out'
# ./pbsim --data-type CLR --length-min 900 --length-max 1100 --accuracy-min 0.84 --accuracy-max 0.86 --model_qc ../data/model_qc_clr --length-mean 1000 --length-sd 100 --accuracy-mean 0.85 --accuracy-sd 0.01 ~/docs/oni/work/NC_000913.fna

lengths = [100, 200, 500, 1000, 2000, 5000, 10000]
error_rates = [0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95]
bandwidths = [16, 24, 32, 40, 48, 64]
xs = [-i for i in range(2, 9)]		# -1..-4
gs = [-i for i in range(2, 13)]		# -1..-6
count = 1000

params_list = [gs, xs, bandwidths, error_rates, lengths]
dimensions = [len(x) for x in params_list]
