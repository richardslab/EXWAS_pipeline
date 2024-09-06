#!/bin/bash
#SBATCH --job-name=exwas_masks  # Specify the job name
#SBATCH --nodes=1         # Specify the number of nodes
#SBATCH --ntasks=1        # Specify the number of tasks (typically used for MPI jobs)
#SBATCH --cpus-per-task=15 # Specify the number of CPUs per task
#SBATCH --mem=64G          # Specify the memory limit per node
#SBATCH --time=10:00:00   # Specify the maximum walltime for the job
#SBATCH --output=/home/richards/ethan.kreuzer/projects/richards/ethan.kreuzer/EXWAS/job-%A-%a.log

ANNOT_HELPER=$1
INPUT=$2

#cd /home/richards/ethan.kreuzer/projects/richards/ethan.kreuzer/EXWAS
cd  /scratch/richards/ethan.kreuzer

#CREATE MASK FILE

echo -e "M1\t$INPUT" > mask.def.regeneron.txt

#MUST filter to the 5 algorithms we want

{
    read -r header

    while IFS=$'\t' read -r id gene consequence plugins num_deleterious CADD_value; do
        
	declare -A plugin_consequence_array=();

        IFS=',' read -r -a plugin_list <<< "$plugins"

	IFS=',' read -r -a consequence_list <<< "$consequence"

        for plugin in "${plugin_list[@]}"; do
            plugin_consequence_array["$plugin"]=1
        done

	for consequence in "${consequence_list[@]}"; do
            plugin_consequence_array["$consequence"]=1
        done

    	if [[ -n "${plugin_consequence_array["$INPUT"]}" ]]; then

            echo -e "$id\t$gene\t$INPUT" >> regenie.annotation.regeneron.txt
	fi




    done


} < "${ANNOT_HELPER}"

