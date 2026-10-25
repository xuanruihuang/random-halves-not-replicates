#!/bin/bash
# Paths: set WORK_ROOT (GPU server work folder) and STORE_ROOT (storage server folder) before running.
# cell-eval ceilings for CD4 donors 2-4, resting (run_celleval_ceiling.py, as donor 1); subsets from the A100, deleted after use
# Usage: run_celleval_d234.sh <d2|d3|d4> [threads]
d=$1; T=${2:-8}; B=${WORK_ROOT}/cross_context; O=$B/celleval_$d; mkdir -p $O
A100=${A100_SSH:?set A100_SSH to the ssh command that reaches the storage server}
f=cd4_${d}_rest_subset.npz; src=${STORE_ROOT}/review_tests/celleval/$f
[ -s $O/$f ] || $A100 "cat $src" > $O/$f
a=$($A100 "md5sum $src" | cut -d' ' -f1); b=$(md5sum $O/$f | cut -d' ' -f1)
[ "$a" = "$b" ] || { echo "$(date) $d md5 mismatch" >> $B/celleval_d234.log; exit 1; }
export OMP_NUM_THREADS=$T OPENBLAS_NUM_THREADS=$T MKL_NUM_THREADS=$T POLARS_MAX_THREADS=$T NUMPY_MADVISE_HUGEPAGE=0
$B/env_celleval/bin/python $B/code/run_celleval_ceiling.py $O/$f $O $T > $O/run.log 2>&1 && { rm -f $O/$f; echo "$(date) DONE $d" >> $B/celleval_d234.log; } || echo "$(date) FAILED $d" >> $B/celleval_d234.log
