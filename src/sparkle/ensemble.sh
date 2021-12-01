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

#compile requirements
feat_suffix="features.csv"
pred_suffix="predictions.csv"
data_name=${10}

#gather requirements
out_file=${11}

#Submit job
/usr/local/envs/sparkles/bin/sparkles --config $sparkles_config \
	sub --no-wait \
	--clean \
	-n $job_name \
	-u '@'$task_filelist \
	--params $task_params \
	python '^src/Python/RunEnsemble.py' \
		--model-config $model_config \
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

# #Increases nodes
/usr/local/envs/sparkles/bin/sparkles addnodes $job_name 60
/usr/local/envs/sparkles/bin/sparkles status --wait $job_name

#fetch results from cloud
mkdir ${data_name}-ensemble-tasks
gsutil -m cp gs://tda/${job_name}/*/*.csv ${data_name}-ensemble-tasks/

#compile tasks per model files
Rscript src/R/compile_ensemble_tasks.R $task_params $feat_suffix $pred_suffix ${data_name}-ensemble-tasks $data_name

#gather ensemble results and compute goodness of fit
Rscript src/R/gather_ensemble.R $task_params $data_name $feat_suffix $pred_suffix $targets $out_file
