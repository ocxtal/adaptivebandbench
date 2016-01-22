#! /usr/bin/env bash

for l in 100 200 500 1000 2000 5000 10000; do ../bench -a $l 1000; done
