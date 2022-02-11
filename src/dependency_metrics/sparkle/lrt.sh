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

{
    read
    while IFS=',', read -r start end
    do
			echo "Start : $start"
			echo "End : $end"
			Rscript 'src/dependency_metrics/LRT.R' $targets $start $end $out_file
		done
} < $task_params


#fetch results from cloud
mkdir -p data/processed/${data_name}-LRT-tasks
cp data/processed/*_LRT_res.csv data/processed/${data_name}-LRT-tasks/

#gather LRT results by rbind individual tasks
Rscript src/dependency_metrics/compile_LRT_tasks.R $task_params data/processed/${data_name}-LRT-tasks $out_file
