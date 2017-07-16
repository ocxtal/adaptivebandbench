# Assessment on Adaptive-Banded Dynamic Programming algorithm for the nucleotide semi-global alignment



This repository contains benchmarking (recall benchmark and speed benchmark) programs and scripts of the adaptive-banded semi-global alignment algorithm. If you are thinking of porting (or implementing) adaptive banded algorithm to your program, you should consider using [GABA library](https://github.com/ocxtal/libgaba) that implements the algorithm with another acceleration algorithm (it is much faster and much more stable).


## Overview of the Adaptive-Banded algorithm

The adaptive-banded algorithm is a fast matrix calculation algorithm for the nucleotide sequence alignments. It is a modification of the conventional banded DP, found in the original implementation of the MOSAIK aligner, the BWA-MEM aligner, and the SeqAn library, making the advancing direction of the band determined adaptively. The advancing direction, right or down, is selected on every vector update, comparing the two cells on the upper and lower edges of the band (see Figure 1(a) below). This method keeps the two cells on the both edges balanced, resulting in capturing the cell with maximum score around the center of the band.

<img src="https://github.com/ocxtal/adaptivebandbench/blob/master/fig/adaptiveband.png">

Figure 1. (a) Vector placement of the band. (b) Vectorized update operation of the affine-gap penalty (Gotoh) algorithm.

### Relation to the semi-global DP routine in the BLAST package

The adaptive-banded algorithm is similar to the semi-global DP routine found in the NCBI BLAST+ package (the gapped alignment phase) in that the algorithm determines the area of the shrinked matrix on the way to calculate it. The difference is that the band width (W in the figure above) is fixed to some constant in the adaptive banded algorithm. The termination of the extension is detected with a X-drop-like test, comparing the center cell of the vector at the head and the maximum score of the center cells of the previous vectors.


### SIMD parallelization

In the adaptive banded algorithm, vectors are placed parallel to the anti-diagonal line of the matrix. The cells having no mutual dependencies and the width of vectors fixed to a constant, the algorithm can be parallelized with Single-Instruction-Multiple-Data (SIMD) operations. The prototype implementation in this repository adopts the SSE4.1 instructions on the x86\_64 processors, calculating eight of 16-bit cells simultaneously. The figure (b) shows the parallelized vector update operation of the affine-gap cost algorithm.



## Assessment of sensitivity

In order to confirm the algorithm reports the optimal scores (and corresponding paths, not tested this time) of the original semi-global alignment algorithm, four experiments are performed.


### Effect of band width

To assess necessary band width, a set of sequence pairs with its mean identity 0.85 and mean length 1000 are given to the algorithm (Figure 2(b)). Reported scores are compared to the result of the naive implementation of the semi-global alignment algorithm. Sum of 1000 trials are shown in the heatmap as [0, 1] normalized recall. The sequence set is generated with [PBSIM](http://bioinformatics.oxfordjournals.org/content/29/1/119.full) read simulator (Figure 2(a)).


### Effect of gap penalties

To assess the performance of the band steering mechanism of the algorithm, recall rates with various gap open and extension penalties are measured. The result shows that the low gap extension penalties, which may cause disappearance of the gradient of scores in vectors, tends to degrade the recall rates (Figure 2(c): the bottom rows of the heatmap).


### Effect of query sequence lengths

The probability of band deviation may increase when the query sequence length gets long. The effect of the length is measured in this experiment, varing the band width 16 and 32. The result (Fig. 2(d)) shows that the 32-cell wide implementation sucessfully captures the optimal scores regardless of the lengths. The 16-cell band is not wide enough to capture long sequences (but applicable to short sequences like Illumina).


### Indel tolerance

The algorithm fundamentally drops alignment paths with large insertions and deletions. The longest acceptable indel lengths are evaluated in this experiment. A set of sequence pairs with its indentity 0.85 and length 1000, and insertion of various lengths is given to the algorithm. The result (Fig. 2(e)) shows that the algortihm is able to recover optimal alignment when the insertion length is less than W - 4.

<img src="https://github.com/ocxtal/adaptivebandbench/blob/master/fig/sensitivitybench.png">

Figure 2. Results of the accuracy assessment

## Speed benchmark

Calculation time to report the maximum score is measured on the SIMD-parallelized implementation of the adaptive-banded algorithm and the BLAST gapped alignment routine. The result (Figure 3) revealed that the SSE4.1 adaptive banded implementation was 7 times faster than the BLAST DP on an Intel Ivy Bridge processor.

<img src="https://github.com/ocxtal/adaptivebandbench/blob/master/fig/speedbench.png" width="400">

Figure 3. Speed benchmark (sum of 1000 runs).

## Running scripts

### Contents

Several algorithms calculating semi-global alignment and benchmarking scripts are included.

* Naive, full-sized semi-global alignment with affine-gap penalty model.
* Adaptive banded DP with affine-gap penalty. (acceptable bandwidth is multiple of 8, determined at compile time with -DBW=32)
* Re-implementation of the semi-gapped alignment function in the NCBI BLAST+ package.
* SIMD parallelized variant of the BLAST semi-gapped alignment function.
* Myers' wavefront algorithm (with some heuristics, described in the DALIGNER paper), extracted from the [DALIGNER](https://github.com/thegenemyers/DALIGNER) repository.


### Build

`make` will generate all the binaries needed in the benchmarks with gcc compiler.


### Recall benchmarks


#### Dependencies

The benchmark script are written in python with numpy library and internally calls [PBSIM reads simulator](http://bioinformatics.oxfordjournals.org/content/29/1/119.short). Optionally it uses Sun Grid Engine to run tasks on cluster computing environments. Figure generation script depends on python with matplotlib and R (rscript command) with ggplot2.


#### General evaluation

Running `scripts/evaluate.py` generates recall rate evaluation with various sequence lengths, sequence identities, band widths (of the adaptive banded algorithm) and mismatch / gap penalties. The results of the three figures, BW-identity plot and length-identity plot with BW=16 and BW=32, are generated with this script. The script invokes naive algorithm and dynamic banded algorithm in `../bin` directory (binaries must be generated beforehand) and puts the results into a tsv-formatted file.

The parametes, a set of ranges of sequence lengths, sequence identities, band widths, and mismatch / gap penalties, are stored in `scripts/params.py`. `generate_params.py` will calculate all the combination of the parameters and print them. The list of combination of parameters is needed as an input parameter of `evaluate.py`.

The results generated by `evaluate.py` is a large tsv-formatted file, can be aggregated (and shrinked) into python-readable text file with `aggregate.py`. The output will contain multi-dimensional (dimension the same as the input parameters) list in python format. You can convert the content into python list with `eval` function.

```
$ python scripts/generate_params.py > params.tsv
$ python scripts/evaluate.py params.tsv out.tsv
$ python scripts/aggregate.py out.tsv result.txt
```

In python, `result.txt` can be converted to `numpy.array` object. `l` and `a` holds the results of linear-gap and affine-gap algorithms, respectively.

 
```
>>> from scripts.util import *
>>> (l, a) = load_result('result.txt')
```

The whole evaluation task is very heavy (mainly due to the slow matrix calculation in the naive, full-sized implementation), you can run scripts in a distributed mannar with Sun Grid Engines. The `sim.sh` script carry outs parameter generation and evaluation task distribution with `qsub` command. Multiple output files will be generated, prefixed with task numbers, can be aggregated into `result.txt` with `aggr.sh`.


```
$ cd path/to/work
$ ../path/to/scripts/sim.sh
$ ../path/to/scripts/aggr.sh
```

The `results/` directory contains aggregated results of the benchmark with gap extension penalty Ge = 0 (linear), 1, 2, and 3.


#### Gap penalty evaluation

The figure 3(b), recall rate with different gap panalties, was generated by other scripts, `generate_params_gige.py`, `evaluate_gige.py` and `aggregate_id_gige.py`. The script runs similar evaluation tasks to the general evaluation, different in that the mismatch penalty is fixed to 2 and both gap open and gap extension penalties are varied. The results are placed in `results/`, named `result_gap_linear.txt` and `result_gap_affine.txt`.


#### Indel tolerance evaluation

The figure 3(d), recall rate with different gap insertion size, was generated with `generate_params_gap.py`, `evaluate_gap.py` and `aggregate_gap.py`. The scripts generates simulated reference-read pairs with a gap inserted on the reference sequence and recall rate evaluation with parameters similar to the general evaluation. The result is placed in `results/` with the name `gige_id_ddiag_table.txt`.



### Speed benchmark

`scripts/speed.sh` invokes `../bin/bench` with sequences in `../seq` with lengths varied within [1, 10k].


### Figure generation


The results used in the paper is contained in results/ directory and plotted to eps files with  `plot.py` and `plot_heatmap.R`. The `plot.py` scripts open the result files in `../results` directory and perform some preprocesses required in `plot_heatmap.R`. Another script `plot_spd.R` generates a length-calc.time plotting from a result of the `speed.sh`

Plotting heatmap:

```
$ cd results
$ python ../scripts/plot.py
$ Rscript ../scripts/plot.R
```

Plotting graph:

```
$ Rscript plot_spd.R result.tsv result.pdf
```


## License

Apache v2.

Copyright (c) 2015-2016 Hajime Suzuki






