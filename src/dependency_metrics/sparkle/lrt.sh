#!/bin/bash
set -e

#sparklespray requirements
sparkles_config=$1 
job_name=$2
task_filelist=$3
task_params=$4

#LRT requirements
targets=$5

#gather requirements
data_name=$6
out_file=$7

#Submit job
/usr/local/envs/sparkles/bin/sparkles --config $sparkles_config \
	sub --no-wait \
	--clean \
	-n $job_name \
	-u '@'$task_filelist \
	--params $task_params \
	Rscript '^src/dependency_metrics/LRT.R' $targets '{start}' '{end}' $out_file

#Increases nodes
/usr/local/envs/sparkles/bin/sparkles addnodes $job_name 20
/usr/local/envs/sparkles/bin/sparkles status --wait $job_name

#fetch results from cloud
mkdir ${data_name}-LRT-tasks
gsutil -m cp gs://tda/${job_name}/*/*.csv ${data_name}-LRT-tasks/

#gather LRT results by rbind individual tasks
Rscript src/dependency_metrics/compile_LRT_tasks.R $task_params ${data_name}-LRT-tasks $out_file
