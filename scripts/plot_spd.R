#! /usr/bin/env Rscript
#
# rows must be:
# len blast simdblast aband wavefront parasail farrar
#
#
# cat seqs/pair.window.20.32.30kover.txt | head -20000 | ./bin/bench -a - $l
#

args <- commandArgs(trailingOnly = T) 
a = read.csv(args[1],sep='\t')

pdf(args[2])

# plot affine-gap results 

a$b40 = a$b40/1000.0
a$b80 = a$b80/1000.0
a$s40 = a$s40/1000.0
a$s80 = a$s80/1000.0
a$a32 = a$a32/1000.0
a$a64 = a$a64/1000.0
a$wave = a$wave/1000.0
a$parasail = a$parasail/1000.0
a$ssw = a$ssw/1000.0

par(pin = c(5,5))

plot(a$len,rep(1,length(a$len)),col='black',type='n',log='xy',xlim=c(0.9,10000),ylim=c(5,10000),xlab='Query sequence length (bp)',ylab='Calculation time (us)')
lines(a$len,a$b40,col='black',lty=3)
lines(a$len,a$b80,col='black',lty=2)
# lines(a$len,a$s40,col='black',lty=6)
# lines(a$len,a$s80,col='black',lty=4)
lines(a$len,a$a32,col='black',lty=5)
lines(a$len,a$a64,col='black',lty=1)
lines(a$len,a$wave,col='black',lty=6)
lines(a$len,a$parasail,col='black',lty=4)

# legend('topleft',legend=c('BLAST (X = 40)', 'BLAST (X = 80)', 'BLAST (SIMD; X = 40)', 'BLAST (SIMD; X = 80)', 'Adaptive (W = 32)', 'Adaptive (W = 64)', 'Wavefront', 'Farrar'),
	# col=c('black','black','black','black','black','black','black','black'),lty=c(3,2,6,4,5,1,5,1))
legend('topleft',legend=c('Adaptive (W = 32)', 'Adaptive (W = 64)', 'BLAST (X = 40)', 'BLAST (X = 80)', 'Wavefront', 'Farrar'),
	col=c('black','black','black','black','black','black'),lty=c(5,1,3,2,6,4))

dev.off()
