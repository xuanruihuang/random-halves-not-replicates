#!/bin/bash
# Paths: set WORK_ROOT (GPU server work folder) and STORE_ROOT (storage server folder) before running.
# Model evaluation, the full State recipe in donors 2-4 (5090): SE-600M embeddings of every pilot cell
# (pilot_se_embed.py, from the all-genes copy of the donor's pilot subset written by run_allgenes_d234.sh on the A100),
# then STATE on the embeddings with the settings of code/pilot/run_state_se.sh (README example, 40,000 steps, batch 8,
# embed_key X_state, output_space gene) and the donor's STATE few-shot split (state_<D>/cd4_fewshot.toml, the dataset
# path pointed at the embedded h5ad); predict, score (pilot_state_score.py).
# Usage: run_5090_se_donor.sh <donor D2|D3|D4> <gpu>. Needs data_<D> (pilot subset) on the 5090. The all-genes files
# (about 14 GB per donor) are deleted once embedded; the embedded h5ad and predictions at the end.
D=$1; G=$2
R=${WORK_ROOT}/model_eval; DD=$R/data_$D; SD=$R/state_$D; ED=$R/state_se_$D; O=$R/out_$D; C=$R/code
PY=${WORK_ROOT}/vcc2026/stateenv/bin/python; ST=${WORK_ROOT}/vcc2026/stateenv/bin/state
A100=${A100_SSH:?set A100_SSH to the ssh command that reaches the storage server}
export PILOT_DONOR=$D CUDA_VISIBLE_DEVICES=$G NUMPY_MADVISE_HUGEPAGE=0 OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8 WANDB_MODE=disabled PYTHONPATH=$C
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True
mkdir -p $ED $O
if [ ! -s $ED/data/cd4_d1_se.h5ad ]; then
  for st in ${D}_Rest ${D}_Stim48hr; do
    f=pilot_${st}_allgenes.npz
    [ -s $ED/emb_$st.npy ] && continue   # already embedded (pilot_se_embed.py keeps emb_<state>.npy)
    [ -s $DD/$f ] || $A100 "cat ${STORE_ROOT}/model_eval/pilot_$D/$f" > $DD/$f
    a=$($A100 "md5sum ${STORE_ROOT}/model_eval/pilot_$D/$f" | cut -d' ' -f1); b=$(md5sum $DD/$f | cut -d' ' -f1)
    [ "$a" = "$b" ] || { echo "$(date) $D $f md5 mismatch"; exit 1; }
  done
  echo "$(date) $D all-genes data copied"
  $PY $C/pilot_se_embed.py $DD ${WORK_ROOT}/pilot/se600m $ED 32 > $ED/embed.log 2>&1 || { echo "$(date) $D embed failed"; exit 1; }
  rm -f $DD/pilot_${D}_*_allgenes.npz
  echo "$(date) $D embedded"
fi
sed "s#$SD/data/#$ED/data/#" $SD/cd4_fewshot.toml > $ED/cd4_fewshot_se.toml
grep -q "$ED/data/" $ED/cd4_fewshot_se.toml || { echo "$(date) $D toml path not replaced"; exit 1; }
cd $ED
$ST tx train \
  data.kwargs.toml_config_path=$ED/cd4_fewshot_se.toml \
  data.kwargs.embed_key=X_state \
  data.kwargs.output_space=gene \
  data.kwargs.pert_col=gene \
  data.kwargs.cell_type_key=cell_type \
  data.kwargs.batch_col=gem_group \
  data.kwargs.control_pert=non-targeting \
  data.kwargs.num_workers=8 \
  training.max_steps=40000 \
  training.batch_size=8 \
  use_wandb=false \
  model=state \
  output_dir=$ED/runs \
  name=cd4_${D}_se > $ED/train.log 2>&1 || { echo "$(date) $D SE train failed"; exit 1; }
$ST tx predict --output-dir runs/cd4_${D}_se --checkpoint final.ckpt --predict-only > $ED/predict.log 2>&1 || { echo "$(date) $D SE predict failed"; exit 1; }
E=$ED/runs/cd4_${D}_se/eval_final.ckpt
$PY $C/pilot_state_score.py $DD $E/adata_pred.h5ad $O/state_se_scores.json > $ED/score.log 2>&1 || { echo "$(date) $D SE score failed"; exit 1; }
rm -f $E/adata_pred.h5ad $E/adata_real.h5ad $ED/data/cd4_d1_se.h5ad
echo "$(date) ${D}_SE_DONE"
