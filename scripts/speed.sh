#! /bin/sh

DBAND_HOME=..

for l in 1 2 5 10 20 50 100 180 300 470 680 1000 1800 3000 4700 6800 10000 18000 30000;
do
	printf "%d\t" $l;
	cat $DBAND_HOME/seqs/pair.window.20.32.30kover.txt | head -20000 | $DBAND_HOME/bin/bench -a - $l;
done
