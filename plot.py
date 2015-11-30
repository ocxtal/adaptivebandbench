
filename = 'bench.csv'

def plot(filename, labels, xs, ys_arr):
	
	import pylab
	pylab.figure(figsize = (15, 9), dpi = 72)

	dashes = [
		[1],
		[8, 4],
		[12, 2, 2, 2],
		[8, 2, 2, 2, 2, 2],
		[4, 2, 2, 2],
		[2, 2]]

	width = 1.5 	# linewidth
	size = 24 		# fontsize
	# gather benchmark results of 16-bit diag and diff algorithms

	for (label, ys, d) in zip(labels, ys_arr, dashes):
		line, = pylab.plot(xs, ys, color = 'black', linewidth = 1.5, linestyle = '-',  label = label)
		line.set_dashes(d)

	# x axis
	pylab.xlim(480, 20200)
	pylab.xticks(fontsize = size)
	pylab.xlabel('Query length (bp)', fontsize = size)
	pylab.xscale('log')

	pylab.yticks(fontsize = size)
	pylab.yscale('log')

	pylab.legend(loc = 'best', fontsize = size)

	# save figure
	pylab.savefig(filename, dpi = 72)
	# pylab.show()


if __name__ == '__main__':
	with open(filename) as f:
		xs = [int(a) for a in f.readline().split(',')[1:]]
		ys_arr_linear = []
		ys_arr_linear += [[int(a) for a in f.readline().split(',')[1:]]]
		ys_arr_linear += [[int(a) for a in f.readline().split(',')[1:]]]
		ys_arr_linear += [[int(a) for a in f.readline().split(',')[1:]]]
		ys_arr_linear += [[int(a) for a in f.readline().split(',')[1:]]]
		ys_arr_linear += [[int(a) for a in f.readline().split(',')[1:]]]

		ys_arr_affine = []
		ys_arr_affine += [[int(a) for a in f.readline().split(',')[1:]]]
		ys_arr_affine += [[int(a) for a in f.readline().split(',')[1:]]]
		ys_arr_affine += [[int(a) for a in f.readline().split(',')[1:]]]
		ys_arr_affine += [[int(a) for a in f.readline().split(',')[1:]]]
		ys_arr_affine += [[int(a) for a in f.readline().split(',')[1:]]]

		print(xs)
		print(ys_arr_linear)
		print(ys_arr_affine)

		labels = ['rognes', 'blast', 'simdblast', 'diag', 'ddiag']
		plot('linear.eps', labels, xs, ys_arr_linear)
		plot('affine.eps', labels, xs, ys_arr_affine)

