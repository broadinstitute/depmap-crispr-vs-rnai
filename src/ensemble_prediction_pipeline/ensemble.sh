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
feat_suffix="features.csv"
pred_suffix="predictions.csv"
data_name=${11}


#Submit job
/usr/local/envs/sparkles/bin/sparkles --config $sparkles_config \
	sub --no-wait \
	--clean \
	-n $job_name \
	-u '@'$task_filelist \
	--params $task_params \
	python '^src/ensemble_prediction_pipeline/RunEnsemble.py' \
		--model-config $model_config \
		--task-mode $task_mode \
		--targets $targets \
		--nfolds $folds \
		--confounders $confounders \
		--feature-dir . \
		--feature-info $file_info \
		--model '{model}' \
		--start-col '{start}' \
		--end-col '{end}' \
		--feat-suffix $feat_suffix \
		--pred-suffix $pred_suffix

#Increases nodes
/usr/local/envs/sparkles/bin/sparkles addnodes $job_name 80
/usr/local/envs/sparkles/bin/sparkles status --wait $job_name

#fetch results from cloud
mkdir ${data_name}-tasks
gsutil -m cp gs://tda/${job_name}/*/*.csv ${data_name}-tasks/

#compile tasks per model files
Rscript src/ensemble_prediction_pipeline/compile_ensemble_tasks.R $task_params $feat_suffix $pred_suffix ${data_name}-tasks $data_name

