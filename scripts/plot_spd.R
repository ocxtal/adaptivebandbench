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
a$wavefront = a$wavefront/10000.0
a$parasail = a$parasail/10000.0
a$farrar = a$farrar/10000.0

par(pin = c(5,5))

plot(a$len,a$blast,col='black',type='n',pch=20,log='xy',xlim=c(0.9,10000),ylim=c(0.5,1000),xlab='Query sequence length (bp)',ylab='Calculation time (ms)')
points(a$len,a$simdblast,col='black',type='n',pch=20)
points(a$len,a$aband,col='black',type='n',pch=20)
points(a$len,a$wavefront,col='black',type='n',pch=20)
points(a$len,a$parasail,col='black',type='n',pch=20)
lines(a$len,a$blast,col='black',lty=2)
lines(a$len,a$simdblast,col='black',lty=4)
lines(a$len,a$aband,col='black',lty=1)
lines(a$len,a$wavefront,col='black',lty=5)
lines(a$len,a$parasail,col='black',lty=3)

legend('bottomright',legend=c('BLAST', 'BLAST (SIMD)', 'Adaptive', 'Wavefront', 'Parasail'),
	col=c('black','black','black','black','black'),lty=c(2,4,1,5,3))

dev.off()

