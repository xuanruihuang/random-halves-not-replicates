#!/bin/bash
# Paths: set WORK_ROOT (GPU server work folder) and STORE_ROOT (storage server folder) before running.
# Model evaluation on the A100 (moved from the 5090, whose shared /data disk filled repeatedly): the same code and
# settings as code/model_eval/run_5090_pub_donor_v3.sh (Arc pseudobulk, PIE) and run_5090_se_donor.sh (full State),
# with the conda environments copied by pack_env_to_a100.sh and the pilot subsets already here (model_eval/pilot_<D>).
# Usage: run_a100_models.sh <donor D2|D3|D4> <pseudobulk|pie|se> [gpu]. Outputs: model_eval/out_<D>/.
D=$1; JOB=$2; export CUDA_VISIBLE_DEVICES=${3:-1}
R=${STORE_ROOT}/model_eval; DD=$R/pilot_$D; SD=$R/state_$D; O=$R/out_$D; C=$R/code_models; PD=$R/pie_$D; ED=$R/state_se_$D
PY=${STORE_ROOT}/envs/stateenv/bin/python; PIE=${STORE_ROOT}/envs/pieenv/bin/pie
ST="$PY -m state"   # the env was copied as a directory, so its console scripts still point at the 5090 path
export PILOT_DONOR=$D NUMPY_MADVISE_HUGEPAGE=0 OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4 WANDB_MODE=disabled PYTHONPATH=$C
mkdir -p $SD $O

prep_state() {   # STATE input and few-shot TOML; the validation guides must equal those of the 5090 STATE run (prep_5090.json)
  $PY $C/pilot_state_prep.py $DD $SD > $SD/prep.log 2>&1 || { echo "$(date) $D state prep failed"; return 1; }
  $PY -c "import json,sys; a=json.load(open('$SD/prep.json'))['val']; b=json.load(open('$SD/prep_5090.json'))['val']; sys.exit(0 if a==b else 1)" \
    || { echo "$(date) $D validation guides differ from the STATE run"; return 1; }
}

case $JOB in
pseudobulk)   # own folder, so that a concurrent full-State run of the same donor cannot delete its input
  mkdir -p ${SD}_arc && cp $SD/prep_5090.json ${SD}_arc/ && SD=${SD}_arc
  prep_state || exit 1
  cd $SD
  $ST tx train data.kwargs.toml_config_path=$SD/cd4_fewshot.toml data.kwargs.embed_key=X_hvg data.kwargs.output_space=gene \
    data.kwargs.pert_col=gene data.kwargs.cell_type_key=cell_type data.kwargs.batch_col=gem_group data.kwargs.control_pert=non-targeting \
    data.kwargs.num_workers=8 training.max_steps=40000 training.batch_size=8 use_wandb=false model=pseudobulk \
    output_dir=$SD/runs name=cd4_${D}_pseudobulk > $SD/train_pseudobulk.log 2>&1 || { echo "$(date) $D pseudobulk train failed"; exit 1; }
  $ST tx predict --output-dir runs/cd4_${D}_pseudobulk --checkpoint final.ckpt --predict-only > $SD/predict_pseudobulk.log 2>&1 || { echo "$(date) $D pseudobulk predict failed"; exit 1; }
  E=$SD/runs/cd4_${D}_pseudobulk/eval_final.ckpt
  $PY $C/pilot_state_score.py $DD $E/adata_pred.h5ad $O/arc_pseudobulk_scores.json > $SD/score_pseudobulk.log 2>&1 || { echo "$(date) $D pseudobulk score failed"; exit 1; }
  rm -f $E/adata_pred.h5ad $E/adata_real.h5ad $SD/data/cd4_d1.h5ad
  echo "$(date) ${D}_PSEUDOBULK_DONE" ;;
pie)
  cp $SD/prep_5090.json $SD/prep.json
  export PIE_ENV_FILE=$PD/common.sh; mkdir -p $PD/cd4
  printf 'WANDB_ENTITY=none\nWANDB_PROJECT=none\nPIE_DATA_ROOT=%s\nPIE_RUNS_ROOT=%s/runs\nPIE_CACHE_DIR=%s/cache\n' $PD $PD $PD > $PD/common.sh
  S=${STORE_ROOT}/models/pie_sources_pub; P=$PD/cd4
  rm -rf $P/de $P/preprocessed $P/preprocessed_withpuro $P/counts $P/expression $PD/runs/cd4_$D
  $PY $C/pilot_pie_prep.py $DD $SD/prep.json $P > $P/inputs.log 2>&1 || { echo "$(date) $D PIE inputs failed"; exit 1; }
  cd $PD
  OMP_NUM_THREADS=8 $PIE process input=cd4/counts/cd4.h5ad de.enabled=true de.groupby=target_gene de.reference=non-targeting \
    de.context_column=cell_line de.output_dir=cd4/de > cd4/process.log 2>&1 || { echo "$(date) $D PIE process failed"; exit 1; }
  OMP_NUM_THREADS=8 $PIE prep name=cd4 label_format=pie_process "labels=cd4/de/*.parquet" h5ad=cd4/expression/cd4.h5ad \
    output_dir=cd4/preprocessed_withpuro obs.context_col=cell_line obs.pert_col=target_gene obs.control_label=non-targeting > cd4/prep.log 2>&1 \
    || { echo "$(date) $D PIE prep failed"; exit 1; }
  ${STORE_ROOT}/envs/pieenv/bin/python -c "import json; g=json.load(open('cd4/preprocessed_withpuro/meta.json'))['genes']; open('cd4/genes_nopuro.txt','w').write('\n'.join(x for x in g if x!='PuroR')+'\n')"
  OMP_NUM_THREADS=8 $PIE prep name=cd4 label_format=pie_process "labels=cd4/de/*.parquet" h5ad=cd4/expression/cd4.h5ad \
    output_dir=cd4/preprocessed genes=cd4/genes_nopuro.txt obs.context_col=cell_line obs.pert_col=target_gene obs.control_label=non-targeting > cd4/prep2.log 2>&1 \
    || { echo "$(date) $D PIE prep2 failed"; exit 1; }
  rm -rf cd4/preprocessed_withpuro cd4/counts cd4/expression
  $PIE train experiment_name=cd4_$D \
    "data.preprocessed_dirs=[$P/preprocessed]" data.dataset_weights=null data.split_dir=$P/splits \
    "data.source_dirs={esm2:$S/esm2,string_space:$S/string_space,depmap_gene_effect:$S/depmap_gene_effect,context_text:$S/context_text,perturbation_text:$S/perturbation_text}" \
    data.gene_text_dir=$S/gene_text data.aliases_path=$P/aliases.yaml \
    data.delta_p.bin_width_fold_change=1.1 data.evidence.seed=42 data.evidence.chunk=256 \
    trainer.devices=1 trainer.num_nodes=1 trainer.max_steps=5000 trainer.accumulate_grad_batches=8 \
    logger.enabled=false > $P/train.log 2>&1 || { echo "$(date) $D PIE train failed"; exit 1; }
  $PIE eval experiment_name=cd4_$D ckpt=best_auprc split_path=$P/splits/test.json row_set=test save_predictions=true > $P/eval.log 2>&1 \
    || { echo "$(date) $D PIE eval failed"; exit 1; }
  ${STORE_ROOT}/envs/pieenv/bin/python $C/pilot_pie_score.py $DD $PD/runs/cd4_$D/eval/test/predictions_best_auprc.parquet $O/pie_scores.json > $P/score.log 2>&1 \
    || { echo "$(date) $D PIE score failed"; exit 1; }
  rm -rf $PD/cache $PD/runs/cd4_$D/checkpoints
  echo "$(date) ${D}_PIE_DONE" ;;
se)
  prep_state || exit 1
  mkdir -p $ED; export OMP_NUM_THREADS=8 PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True
  $PY $C/pilot_se_embed.py $DD ${STORE_ROOT}/models/se600m $ED 32 > $ED/embed.log 2>&1 || { echo "$(date) $D embed failed"; exit 1; }
  sed "s#$SD/data/#$ED/data/#" $SD/cd4_fewshot.toml > $ED/cd4_fewshot_se.toml
  grep -q "$ED/data/" $ED/cd4_fewshot_se.toml || { echo "$(date) $D toml path not replaced"; exit 1; }
  cd $ED
  $ST tx train data.kwargs.toml_config_path=$ED/cd4_fewshot_se.toml data.kwargs.embed_key=X_state data.kwargs.output_space=gene \
    data.kwargs.pert_col=gene data.kwargs.cell_type_key=cell_type data.kwargs.batch_col=gem_group data.kwargs.control_pert=non-targeting \
    data.kwargs.num_workers=8 training.max_steps=40000 training.batch_size=8 use_wandb=false model=state \
    output_dir=$ED/runs name=cd4_${D}_se > $ED/train.log 2>&1 || { echo "$(date) $D SE train failed"; exit 1; }
  $ST tx predict --output-dir runs/cd4_${D}_se --checkpoint final.ckpt --predict-only > $ED/predict.log 2>&1 || { echo "$(date) $D SE predict failed"; exit 1; }
  E=$ED/runs/cd4_${D}_se/eval_final.ckpt
  $PY $C/pilot_state_score.py $DD $E/adata_pred.h5ad $O/state_se_scores.json > $ED/score.log 2>&1 || { echo "$(date) $D SE score failed"; exit 1; }
  rm -f $E/adata_pred.h5ad $E/adata_real.h5ad $ED/data/cd4_d1_se.h5ad $SD/data/cd4_d1.h5ad
  echo "$(date) ${D}_SE_DONE" ;;
esac
