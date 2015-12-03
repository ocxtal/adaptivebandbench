#! /usr/bin/env python
# encoding: utf-8

# global params
pbsim_path = '/home/suzukihajime/src/pbsim-1.0.3/src/'
ref_path = '/home/suzukihajime/oni/work/resource/NC_000913.fna'
ref_length = 4700000
bin_path = '../bin'
count = 1000

# general evaluation
lengths = [100, 200, 500, 1000, 2000, 5000, 10000]
error_rates = [0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95]
bandwidths = [16, 24, 32, 40, 48, 64]
xs = [-i for i in range(2, 9)]		# -1..-4
gs = [-i for i in range(2, 13)]		# -1..-6
params_list = [gs, xs, bandwidths, error_rates, lengths]

# gap-length evaluation
ref_gaps = [i for i in range(0, 50, 1)]
read_gaps = [0]
bandwidths_gap = [16, 24, 32, 40, 48, 56, 64]
params_list_gap = [ref_gaps, read_gaps, bandwidths_gap, error_rates, [1000, 10000]]

