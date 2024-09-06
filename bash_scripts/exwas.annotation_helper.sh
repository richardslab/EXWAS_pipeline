#!/bin/bash
#SBATCH --job-name=exwas_masks  # Specify the job name
#SBATCH --nodes=1         # Specify the number of nodes
#SBATCH --ntasks=1        # Specify the number of tasks (typically used for MPI jobs)
#SBATCH --cpus-per-task=10 # Specify the number of CPUs per task
#SBATCH --mem=75G          # Specify the memory limit per node
#SBATCH --time=10:00:00   # Specify the maximum walltime for the job
#SBATCH --output=/scratch/richards/ethan.kreuzer/job-%A-%a.log

INPUT_VCF=$1

#cd /project/richards/ethan.kreuzer/EXWAS
cd /scratch/richards/ethan.kreuzer

is_deleterious() {
    local plugin=$1
    local value=$2

    case $plugin in

    	IMPACT)

	    [[ "$value" == "HIGH" ]] && return 0
            ;;
        LoF)
            [[ "$value" == "HC" ]] && return 0
            ;;
        CADD_PHRED)
            return 0
            ;;
        AlphaMissense_pred|EVE_Class25_pred)
            [[ "$value" == *P* ]] && return 0
            ;;
        MutationTaster_pred)
            [[ "$value" == *D* || "$value" == *A* ]] && return 0
            ;;
        LRT_pred|Polyphen2_HDIV_pred|Polyphen2_HVAR_pred|SIFT4G_pred|SIFT_pred)
            [[ "$value" == *D* ]] && return 0
            ;;
    esac
    return 1
}


#Preprocess the file to filter out lines, select desired columns, and save to a new file

grep -v '^##' "${INPUT_VCF}".finalAnnot.filter_vep.txt | awk '{print $1, $4, $7, $14}' | uniq | \

awk '
function getImpactValue(impact) {
    if (impact == "HIGH") return 4;
    if (impact == "MODERATE") return 3;
    if (impact == "LOW") return 2;
    if (impact == "MODIFIER") return 2;
    return 1;
}

{
    split($4, a, ";");
    split(a[1], b, "=");
    impactValue = getImpactValue(b[2]);

    if ($1 != prev_id) {
        if (prev_id != "") print best_line;
        best_line = $0;
        best_impact_value = impactValue;
        prev_id = $1;
    } else {
        if (impactValue > best_impact_value) {
            best_line = $0;
            best_impact_value = impactValue;
        }
    }
}
END {
    if (prev_id != "") print best_line;
}' > "${INPUT_VCF}.finalAnnot.filter_vep.txt.tmp"

# Process the annotated vcf file


echo -e "VARIANT_ID\tGENE\tCONSEQUENCE\tDELETERIOUS\t#_DELETERIOUS\tCADD_PHRED" > regenie.annotation.helper.txt


{
    read -r header

    while IFS=' ' read -r id gene consequence extra; do

        IFS=';' read -r -a plugins <<< "$extra"
        plugin_names=()  # declare array to store plugin names
        counter=0        # declare counter to count number of plugins that report deleterious
        CADD_value=-1    # declare CADD_value outside the loop to reset it for each line

        for plugin in "${plugins[@]}"; do
            plugin_name="${plugin%%=*}"
            plugin_value="${plugin#*=}"
	    

            if is_deleterious "$plugin_name" "$plugin_value"; then
                if [[ "$plugin_name" == "IMPACT" ]]; then
                    plugin_names+=("$plugin_name")

                elif [[ "$plugin_name" == "CADD_PHRED" ]]; then
                    CADD_value="$plugin_value"
                    plugin_names+=("$plugin_name")
                else
                    counter=$((counter + 1))  # increment counter
                    plugin_names+=("$plugin_name")
                fi
            fi
        done

        plugin_names_str=$(IFS=','; echo "${plugin_names[*]}")

        if [[ $(awk "BEGIN {print ($CADD_value == -1)}") -eq 1 ]]; then

	    echo -e "$id\t$gene\t$consequence\t$plugin_names_str\t$counter" >> regenie.annotation.helper.txt
        else

	    echo -e "$id\t$gene\t$consequence\t$plugin_names_str\t$counter\t$CADD_value" >> regenie.annotation.helper.txt
        fi
    done
} < "${INPUT_VCF}.finalAnnot.filter_vep.txt.tmp"

rm "${INPUT_VCF}.finalAnnot.filter_vep.txt.tmp"
