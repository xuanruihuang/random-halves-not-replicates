#!/bin/bash
# Paths: set WORK_ROOT (GPU server work folder) and STORE_ROOT (storage server folder) before running.
# Model evaluation, published models in donors 2-4 (5090): Arc pseudobulk, Arc perturb_mean and PIE, with the donor-1
# pilot's code and settings unchanged (code/pilot/run_state_model.sh, run_pie.sh and the PIE prep commands of the pilot),
# on the donor's pilot subset (run_a100_donors.sh) and the STATE few-shot split of run_5090_donor.sh (same TOML; the
# STATE input h5ad is rebuilt by pilot_state_prep.py, whose validation guides are checked against the saved prep.json).
# Usage: run_5090_pub_donor.sh <donor D2|D3|D4> <gpu pseudobulk> <gpu perturb_mean> <gpu PIE> [only: pie|arc].
# Outputs: out_<D>/arc_pseudobulk_scores.json, arc_perturb_mean_scores.json, pie_scores.json. Large intermediates are
# deleted at the end (the 5090 disk is nearly full).
D=$1; G1=$2; G2=$3; G3=$4; ONLY=${5:-all}
R=${WORK_ROOT}/model_eval; DD=$R/data_$D; SD=$R/state_$D; O=$R/out_$D; C=$R/code; PD=$R/pie_$D
PY=${WORK_ROOT}/vcc2026/stateenv/bin/python; ST=${WORK_ROOT}/vcc2026/stateenv/bin/state; PIE=${WORK_ROOT}/pilot/pieenv/bin/pie
A100=${A100_SSH:?set A100_SSH to the ssh command that reaches the storage server}
export PILOT_DONOR=$D NUMPY_MADVISE_HUGEPAGE=0 OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4 WANDB_MODE=disabled PYTHONPATH=$C
mkdir -p $DD $SD $O $PD/cd4
for f in pilot_${D}_Rest.npz pilot_${D}_Stim48hr.npz pilot_design.json; do
  [ -s $DD/$f ] || $A100 "cat ${STORE_ROOT}/model_eval/pilot_$D/$f" > $DD/$f
  a=$($A100 "md5sum ${STORE_ROOT}/model_eval/pilot_$D/$f" | cut -d' ' -f1); b=$(md5sum $DD/$f | cut -d' ' -f1)
  [ "$a" = "$b" ] || { echo "$(date) $f md5 mismatch"; exit 1; }
done
echo "$(date) $D data copied"
if [ $ONLY != pie ]; then   # PIE needs only prep.json (its validation guides); the STATE input is left alone
  cp $SD/prep.json $SD/prep_saved.json
  $PY $C/pilot_state_prep.py $DD $SD || { echo "$(date) state prep failed"; exit 1; }
  $PY -c "import json,sys; a=json.load(open('$SD/prep.json'))['val']; b=json.load(open('$SD/prep_saved.json'))['val']; sys.exit(0 if a==b else 1)" \
    || { echo "$(date) validation guides differ from the STATE run"; exit 1; }
  echo "$(date) $D state input rebuilt"
fi

arc() {   # <model> <gpu>: train as run_state_model.sh, predict, score
  M=$1; cd $SD
  CUDA_VISIBLE_DEVICES=$2 $ST tx train \
    data.kwargs.toml_config_path=$SD/cd4_fewshot.toml \
    data.kwargs.embed_key=X_hvg \
    data.kwargs.output_space=gene \
    data.kwargs.pert_col=gene \
    data.kwargs.cell_type_key=cell_type \
    data.kwargs.batch_col=gem_group \
    data.kwargs.control_pert=non-targeting \
    data.kwargs.num_workers=8 \
    training.max_steps=40000 \
    training.batch_size=8 \
    use_wandb=false \
    model=$M \
    output_dir=$SD/runs \
    name=cd4_${D}_$M > $SD/train_$M.log 2>&1 || { echo "$(date) $M train failed"; return 1; }
  CUDA_VISIBLE_DEVICES=$2 $ST tx predict --output-dir runs/cd4_${D}_$M --checkpoint final.ckpt --predict-only > $SD/predict_$M.log 2>&1 || { echo "$(date) $M predict failed"; return 1; }
  E=$SD/runs/cd4_${D}_$M/eval_final.ckpt
  $PY $C/pilot_state_score.py $DD $E/adata_pred.h5ad $O/arc_${M}_scores.json > $SD/score_$M.log 2>&1 || { echo "$(date) $M score failed"; return 1; }
  rm -f $E/adata_pred.h5ad $E/adata_real.h5ad
  echo "$(date) $D $M done"
}

pie() {   # PIE as in the donor-1 pilot: prep inputs, DE (pie process), labels (pie prep, then again without PuroR), train, eval, score
  export PIE_ENV_FILE=$PD/common.sh CUDA_VISIBLE_DEVICES=$G3
  printf 'WANDB_ENTITY=none\nWANDB_PROJECT=none\nPIE_DATA_ROOT=%s\nPIE_RUNS_ROOT=%s/runs\nPIE_CACHE_DIR=%s/cache\n' $PD $PD $PD > $PD/common.sh
  S=${WORK_ROOT}/pilot/pie_data/sources_pub; P=$PD/cd4
  rm -rf $P/de $P/preprocessed $P/preprocessed_withpuro $P/counts $P/expression   # leftovers of an earlier attempt
  $PY $C/pilot_pie_prep.py $DD $SD/prep.json $P > $P/inputs.log 2>&1 || { echo "$(date) PIE inputs failed"; return 1; }
  cd $PD
  OMP_NUM_THREADS=8 $PIE process input=cd4/counts/cd4.h5ad de.enabled=true de.groupby=target_gene de.reference=non-targeting \
    de.context_column=cell_line de.output_dir=cd4/de > cd4/process.log 2>&1 || { echo "$(date) PIE process failed"; return 1; }
  OMP_NUM_THREADS=8 $PIE prep name=cd4 label_format=pie_process "labels=cd4/de/*.parquet" h5ad=cd4/expression/cd4.h5ad \
    output_dir=cd4/preprocessed_withpuro obs.context_col=cell_line obs.pert_col=target_gene obs.control_label=non-targeting > cd4/prep.log 2>&1 \
    || { echo "$(date) PIE prep failed"; return 1; }
  ${WORK_ROOT}/pilot/pieenv/bin/python -c "import json; g=json.load(open('cd4/preprocessed_withpuro/meta.json'))['genes']; open('cd4/genes_nopuro.txt','w').write('\n'.join(x for x in g if x!='PuroR')+'\n')"
  OMP_NUM_THREADS=8 $PIE prep name=cd4 label_format=pie_process "labels=cd4/de/*.parquet" h5ad=cd4/expression/cd4.h5ad \
    output_dir=cd4/preprocessed genes=cd4/genes_nopuro.txt obs.context_col=cell_line obs.pert_col=target_gene obs.control_label=non-targeting > cd4/prep2.log 2>&1 \
    || { echo "$(date) PIE prep2 failed"; return 1; }
  rm -rf cd4/preprocessed_withpuro cd4/counts cd4/expression
  $PIE train experiment_name=cd4_$D \
    "data.preprocessed_dirs=[$P/preprocessed]" data.dataset_weights=null data.split_dir=$P/splits \
    "data.source_dirs={esm2:$S/esm2,string_space:$S/string_space,depmap_gene_effect:$S/depmap_gene_effect,context_text:$S/context_text,perturbation_text:$S/perturbation_text}" \
    data.gene_text_dir=$S/gene_text data.aliases_path=$P/aliases.yaml \
    data.delta_p.bin_width_fold_change=1.1 data.evidence.seed=42 data.evidence.chunk=256 \
    trainer.devices=1 trainer.num_nodes=1 trainer.max_steps=5000 trainer.accumulate_grad_batches=8 \
    logger.enabled=false > $P/train.log 2>&1 || { echo "$(date) PIE train failed"; return 1; }
  $PIE eval experiment_name=cd4_$D ckpt=best_auprc split_path=$P/splits/test.json row_set=test save_predictions=true > $P/eval.log 2>&1 \
    || { echo "$(date) PIE eval failed"; return 1; }
  ${WORK_ROOT}/pilot/pieenv/bin/python $C/pilot_pie_score.py $DD $PD/runs/cd4_$D/eval/test/predictions_best_auprc.parquet $O/pie_scores.json > $P/score.log 2>&1 \
    || { echo "$(date) PIE score failed"; return 1; }
  echo "$(date) $D PIE done"
}

r1=0; r2=0; r3=0
[ $ONLY != pie ] && { arc pseudobulk $G1 & p1=$!; arc perturb_mean $G2 & p2=$!; }
[ $ONLY != arc ] && { rm -rf $PD/runs/cd4_$D; pie & p3=$!; }
[ $ONLY != pie ] && { wait $p1; r1=$?; wait $p2; r2=$?; }
[ $ONLY != arc ] && { wait $p3; r3=$?; }
[ $ONLY != pie ] && rm -f $SD/data/cd4_d1.h5ad
rm -rf $PD/cache $PD/runs/cd4_$D/checkpoints
[ $r1 -eq 0 ] && [ $r2 -eq 0 ] && [ $r3 -eq 0 ] && echo "$(date) ${D}_PUB_DONE" || echo "$(date) ${D}_PUB_FAILED $r1 $r2 $r3"   # data_<D> is kept for scGen and State
