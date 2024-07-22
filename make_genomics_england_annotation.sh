#!/bin/bash
#SBATCH --job-name=exwas_masks  # Specify the job name
#SBATCH --nodes=1         # Specify the number of nodes
#SBATCH --ntasks=1        # Specify the number of tasks (typically used for MPI jobs)
#SBATCH --cpus-per-task=15 # Specify the number of CPUs per task
#SBATCH --mem=64G          # Specify the memory limit per node
#SBATCH --time=12:00:00   # Specify the maximum walltime for the job
#SBATCH --output=/home/richards/ethan.kreuzer/projects/richards/ethan.kreuzer/EXWAS/job-%A-%a.log

ANNOT_HELPER=$1

cd /home/richards/ethan.kreuzer/projects/richards/ethan.kreuzer/EXWAS


#CREATE MASK FILE

echo -e "M1\tLoF" > mask.def.genomics_england.txt
echo -e "M2\tLoF,CADD>20" >> mask.def.genomics_england.txt
echo -e "M3\tLoF,CADD>20,CADD>10" >> mask.def.genomics_england.txt




{
    read -r header

    while IFS=$'\t' read -r id gene consequence plugins num_deleterious CADD_value; do


	declare -A plugin_consequence_array=();

        IFS=',' read -r -a plugin_list <<< "$plugins"

	IFS=',' read -r -a consequence_list <<< "$consequence"

	for consequence in "${consequence_list[@]}"; do
            plugin_consequence_array["$consequence"]=1
        done


        for plugin in "${plugin_list[@]}"; do
            if [[ "$plugin" == "LoF" ]] || [[ "$plugin" == "CADD_PHRED" ]]; then
                plugin_consequence_array["$plugin"]=1
            fi
        done

        if [[ -n "${plugin_consequence_array[LoF]}" ]]; then
            echo -e "$id\t$gene\tLoF" >> regenie.annotation.genomics_england.txt
        elif awk -v val="$CADD_value" 'BEGIN {exit !(val > 20)}'; then

		if [[ -n "${plugin_consequence_array[missense_variant]}" ]]; then

            		echo -e "$id\t$gene\tCADD>20" >> regenie.annotation.genomics_england.txt

		fi
        elif awk -v val="$CADD_value" 'BEGIN {exit !(val > 10)}'; then

		if [[ -n "${plugin_consequence_array[missense_variant]}" ]]; then

                        echo -e "$id\t$gene\tCADD>10" >> regenie.annotation.genomics_england.txt

                fi
        fi

    done
} < "${ANNOT_HELPER}"

