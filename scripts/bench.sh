#! /usr/bin/env bash

for l in 100 150 250 350 500 650 800 1000 1500 2500 3500 5000 6500 8000 10000; do ../bin/bench -a $l $1; done
