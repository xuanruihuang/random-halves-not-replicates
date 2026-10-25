#!/bin/bash
# Paths: set WORK_ROOT (GPU server work folder) and STORE_ROOT (storage server folder) before running.
# A100: all-genes copies of the donor 2-4 pilot subsets (same cells; cd4_pilot_extract.py with PILOT_ALL_GENES=1), for State Embedding
cd ${STORE_ROOT}; export PYTHONPATH=${STORE_ROOT}/pylib OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PILOT_ALL_GENES=1
ME=${STORE_ROOT}/model_eval
for d in D2 D3 D4; do
  PILOT_DONOR=$d python3 $ME/code/cd4_pilot_extract.py $ME/sc_$d $ME/select_$d.npz $ME/pilot_$d 24 > $ME/pilot_$d/allgenes.log 2>&1 && echo "DONE $d" >> $ME/allgenes_status.log || echo "FAILED $d" >> $ME/allgenes_status.log
done
echo ALL_DONE >> $ME/allgenes_status.log
