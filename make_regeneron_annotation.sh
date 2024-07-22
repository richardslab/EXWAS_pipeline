#!/bin/bash
#SBATCH --job-name=exwas_masks  # Specify the job name
#SBATCH --nodes=1         # Specify the number of nodes
#SBATCH --ntasks=1        # Specify the number of tasks (typically used for MPI jobs)
#SBATCH --cpus-per-task=15 # Specify the number of CPUs per task
#SBATCH --mem=64G          # Specify the memory limit per node
#SBATCH --time=10:00:00   # Specify the maximum walltime for the job
#SBATCH --output=/home/richards/ethan.kreuzer/projects/richards/ethan.kreuzer/EXWAS/job-%A-%a.log

ANNOT_HELPER=$1

#cd /home/richards/ethan.kreuzer/projects/richards/ethan.kreuzer/EXWAS
cd  /scratch/richards/ethan.kreuzer

#CREATE MASK FILE

echo -e "M1\tpLoF" > mask.def.regeneron.txt
echo -e "M2\tpLoF,missense.0in5,missense.1in5,missense.5in5" >> mask.def.regeneron.txt
echo -e "M3\tpLoF,missense.5in5" >> mask.def.regeneron.txt
echo -e "M4\tpLoF,missense.1in5,missense.5in5" >> mask.def.regeneron.txt


echo -e "VARIANT_ID\tGENE\tCONSEQUENCE\tDELETERIOUS\t#_DELETERIOUS" > "${ANNOT_HELPER}.filtered"

{
    read -r header

    while IFS=$'\t' read -r id gene consequence plugins num_deleterious CADD_value; do
        IFS=',' read -r -a plugin_list <<< "$plugins"

        filtered_names=()
        count=0

        for plugin in "${plugin_list[@]}"; do
            if [[ "$plugin" == "SIFT_pred" ]] || [[ "$plugin" == "LRT_pred" ]] || [[ "$plugin" == "MutationTaster_pred" ]] || [[ "$plugin" == "Polyphen2_HVAR_pred" ]] || [[ "$plugin" == "Polyphen2_HDIV_pred" ]]; then
                count=$((count + 1))
                filtered_names+=("$plugin")
            else 

	    	if [[ "$plugin" == "IMPACT" ]]; then

	    		filtered_names+=("$plugin")

	    	fi
	    fi


        done

        plugin_names_str=$(IFS=','; echo "${filtered_names[*]}")
        echo -e "$id\t$gene\t$consequence\t$plugin_names_str\t$count" >> "${ANNOT_HELPER}.filtered"

    done
} < "$ANNOT_HELPER"


#MUST filter to the 5 algorithms we want

> regenie.annotation.regeneron.txt.tmp
> regenie.annotation.regeneron.txt

{

    while IFS=$'\t' read -r id gene consequence plugins num_deleterious; do

        declare -A plugin_consequence_array=()

        IFS=',' read -r -a plugin_list <<< "$plugins"
        IFS=',' read -r -a consequence_list <<< "$consequence"
        IFS=$' \t\n'

        for plugin in "${plugin_list[@]}"; do
            plugin_consequence_array["$plugin"]=1
        done

        for consequence in "${consequence_list[@]}"; do
            plugin_consequence_array["$consequence"]=1
        done


        if [[ -n "${plugin_consequence_array[IMPACT]}" ]]; then
            echo -e "$id\t$gene\tpLoF" >> regenie.annotation.regeneron.txt.tmp
        else
            if [[ -n "${plugin_consequence_array[missense_variant]}" ]] && [[ "$num_deleterious" -eq 5 ]]; then
                echo -e "$id\t$gene\tmissense.5in5" >> regenie.annotation.regeneron.txt.tmp
            elif [[ -n "${plugin_consequence_array[missense_variant]}" ]] && [[ "$num_deleterious" -ge 1 ]]; then
                echo -e "$id\t$gene\tmissense.1in5" >> regenie.annotation.regeneron.txt.tmp
            elif [[ -n "${plugin_consequence_array[missense_variant]}" ]]; then
                echo -e "$id\t$gene\tmissense.0in5" >> regenie.annotation.regeneron.txt.tmp
            fi
        fi

    done

} < "${ANNOT_HELPER}.filtered"

awk '!seen[$1]++' regenie.annotation.regeneron.txt.tmp > regenie.annotation.regeneron.txt
rm regenie.annotation.regeneron.txt.tmp
rm "${ANNOT_HELPER}.filtered"

