#!/bin/bash
set -e

#sparklespray requirements
sparkles_config=$1
job_name=$2
task_filelist=$3
task_params=$4

#ensemble requirements
model_config=$5
targets=$6
folds=$7
confounders=$8
file_info=$9
task_mode=${10}

#compile requirements
feat_suffix="features"
pred_suffix="predictions"
data_name=${11}
target=${12}
output_dir=${13}

{
    read
    while IFS=',', read -r start end model
    do
			echo "Start : $start"
			echo "End : $end"
			echo "Model : $model"
			python3 'src/ensemble_prediction_pipeline/RunEnsemble.py' \
				--model-config $model_config \
				--task-mode $task_mode \
				--targets $targets \
				--nfolds $folds \
				--confounders $confounders \
				--feature-dir . \
				--feature-info $file_info \
				--model $model \
				--start-col $start \
				--end-col $end \
				--feat-suffix $feat_suffix \
				--pred-suffix $pred_suffix \
        --output-dir $output_dir
    done
} < $task_params

#Move results from proper folder
#mv data/temp data/processed/$data_name-$target

#compile tasks per model files
Rscript src/ensemble_prediction_pipeline/compile_ensemble_tasks.R $task_params $feat_suffix $pred_suffix data/processed/$data_name-$target data/processed/$data_name $target
