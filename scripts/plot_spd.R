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

a$blast = a$blast/10000.0
a$simdblast = a$simdblast/10000.0
a$aband = a$aband/10000.0
a$wave = a$wave/10000.0
a$parasail = a$parasail/10000.0

par(pin = c(3,3))

plot(a$len,rep(1,length(a$len)),col='black',type='n',log='xy',xlim=c(0.9,10000),ylim=c(5,10000),xlab='Query sequence length (bp)',ylab='Calculation time (us)')
lines(a$len,a$blast,col='black',lty=3)
lines(a$len,a$simdblast,col='black',lty=5)
lines(a$len,a$aband,col='black',lty=1)
lines(a$len,a$wave,col='black',lty=6)
lines(a$len,a$parasail,col='black',lty=4)

legend('topleft',legend=c('Adaptive', 'BLAST', 'BLAST (SIMD)', 'Wavefront', 'Farrar'),
	col=c('black','black','black','black','black'),lty=c(1,3,5,6,4))

dev.off()

