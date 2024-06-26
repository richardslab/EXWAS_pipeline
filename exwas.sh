#!/bin/bash
#SBATCH --job-name=exwas_annotation  # Specify the job name
#SBATCH --nodes=1         # Specify the number of nodes
#SBATCH --ntasks=1        # Specify the number of tasks (typically used for MPI jobs)
#SBATCH --cpus-per-task=15 # Specify the number of CPUs per task
#SBATCH --mem=64G          # Specify the memory limit per node
#SBATCH --time=24:00:00   # Specify the maximum walltime for the job
#SBATCH --output=/home/richards/ethan.kreuzer/projects/richards/ethan.kreuzer/EXWAS/job-%A-%a.log

INPUT_VCF=$1
mask_file=$2
LOF_threshold=$3
CADD_threshold=$4

cd /home/richards/ethan.kreuzer/projects/richards/ethan.kreuzer/EXWAS
module load tabix
module load bcftools
module load apptainer/1.1.8
sif=/scratch/richards/guillaume.butler-laporte/WGS/COVID19_GenOMICC_AVT_analysis/vep_v105.sif


# Normalize, left-align and subset ancestry of the input VCF
bcftools norm -m -any --check-ref w -f /scratch/richards/ethan.kreuzer/vep/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna "${INPUT_VCF}" -Ou | bcftools annotate --set-id '%CHROM:%POS:%REF:%FIRST_ALT' -Ou | bcftools view --drop-genotypes -Oz  > "${INPUT_VCF}".set_id.no_genotypes

tabix -p vcf "${INPUT_VCF}".set_id.no_genotypes


chromosome=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X)
chr="${chromosome[SLURM_ARRAY_TASK_ID]}"

apptainer run --bind ${PWD}:${PWD} ${sif} vep -i "${INPUT_VCF}".set_id.no_genotypes \
         --assembly GRCh38 \
         --format vcf \
         --cache \
         --dir_cache /scratch/richards/ethan.kreuzer/vep \
         -o "${INPUT_VCF}".finalAnnot.txt \
         --plugin LoF,loftee_path:/opt/micromamba/share/ensembl-vep-105.0-1,human_ancestor_fa:/scratch/richards/ethan.kreuzer/vep/human_ancestor.fa.gz,conservation_file:/scratch/richards/ethan.kreuzer/vep/loftee.sql,gerp_bigwig:/scratch/richards/ethan.kreuzer/vep/gerp_conservation_scores.homo_sapiens.GRCh38.bw \
         --plugin CADD,/scratch/richards/ethan.kreuzer/vep/whole_genome_SNVs.tsv.gz,/scratch/richards/ethan.kreuzer/vep/gnomad.genomes.r3.0.indel.tsv.gz \
         --plugin dbNSFP,/scratch/richards/ethan.kreuzer/vep/dbNSFP4.8a_grch38.gz,Ensembl_transcriptid,EVE_Class25_pred,VEP_canonical,LRT_pred,SIFT_pred,SIFT4G_pred,MutationTaster_pred,Polyphen2_HDIV_pred,Polyphen2_HVAR_pred,AlphaMissense_pred \
	 --everything \
         --force_overwrite \
         --offline \
         --fork 1 \
         --quiet

#Make set-list input for regenie

transpose() {
    awk '
    BEGIN {
        FS = "\n";
        max_x = 0;
        max_y = 0;
    }

    {
        max_y++;
        for (i = 1; i <= NF; i++) {
            if (i > max_x) max_x = i;
            A[i, max_y] = $i;
        }
    }

    END {
        for (x = 1; x <= max_x; x++) {
            for (y = 1; y <= max_y; y++) {
                if ((x, y) in A) printf "%s", A[x, y];
                if (y != max_y) printf " ";
            }
            printf "\n";
        }
    }
    '
}

awk '{ print $4, $1 }' "${INPUT_VCF}".finalAnnot.txt | sort -u -k 1.5,1 -k 2.6,2 > canon.chr2.var.gene.txt

awk '{if(!seen[$1]++) { match($2, /^chr([0-9]+)/, arr); print $1, arr[1], ++line }}' canon.chr2.var.gene.txt > canon.chr2.gene.txt

# Read genes into an array
genes=($(awk '{print $1}' canon.chr2.gene.txt | tr '\n' ' '))

# Iterate through each gene
for gene in "${genes[@]}"; do
  if [ -n "$gene" ]; then
    variants=$(awk -v a="$gene" '$1 == a {print $2}' canon.chr2.var.gene.txt | uniq | transpose | tr -s '[:blank:]' "," | awk 'BEGIN{FS=",";OFS=","} { $1=$1; print $0 }')
    first_variant=$(echo $variants | awk -F, '{print $1}')
    position=$(echo $first_variant | awk -F: '{print $2}')
    chromosome=$(echo $first_variant | awk -F: '{print $1}')
    echo "$gene $chromosome $position $variants" >> regenie.set.list.txt.tmp
  fi
done

sed -E 's/([^ ]+ [^ ]+ )[0-9]+chr/\1chr/' regenie.set.list.txt.tmp  > regenie.set.list.txt

rm regenie.set.list.txt.tmp


#Make annotation file for regenie input


# Read the mask file into an array
mapfile -t mask_selection < "$mask_file"

# Function to check if a plugin value is deleterious
is_deleterious() {
    local plugin=$1
    local value=$2

    case $plugin in
        LoF)
            [[ "$value" == ${LOF_threshold} ]] && return 0
            ;;
        CADD_PHRED)
            (( $(echo "$value >= ${CADD_threshold}" | bc -l) )) && return 0
            ;;
        AlphaMissense_pred|EVE_Class25_pred)
            [[ "$value" == *P* ]] && return 0
            ;;
        LRT_pred)
            [[ "$value" == *D* ]] && return 0
            ;;
        MutationTaster_pred)
            [[ "$value" == *D* || "$value" == *A* ]] && return 0
            ;;
        Polyphen2_HDIV_pred|Polyphen2_HVAR_pred|SIFT4G_pred|SIFT_pred)
            [[ "$value" == *D* ]] && return 0
            ;;
    esac
    return 1
}

# Function to check if a plugin or consequence is in the mask selection
in_mask_selection() {
    local item=$1
    for mask in "${mask_selection[@]}"; do
        if [[ "$mask" == "$item" ]]; then
            return 0
        fi
    done
    return 1
}

# Process the input file
{
    # Skip the header line
    read -r header

    # Process each line
    while IFS=$'\t' read -r id location allele gene feature feature_type consequence cDNA_position CDS_position protein_position amino_acids codons existing_variation extra; do
        # Skip lines where gene is "-"
        if [[ "$gene" == "-" ]]; then
            continue
        fi


        # Split the consequences and print each one if it is in the mask selection
        IFS=',' read -r -a consequences_array <<< "$consequence"
        for single_consequence in "${consequences_array[@]}"; do
            if in_mask_selection "$single_consequence"; then
                echo -e "$id\t$gene\t$single_consequence" >> regenie.annotation.txt
            fi
        done

        # Iterate over each plugin in Extra field
        IFS=';' read -r -a plugins <<< "$extra"
        for plugin in "${plugins[@]}"; do
            plugin_name=$(echo "$plugin" | cut -d'=' -f1)
            plugin_value=$(echo "$plugin" | cut -d'=' -f2)


            # Check if the plugin value is deleterious and if the plugin name is in the mask selection
            if is_deleterious "$plugin_name" "$plugin_value" && in_mask_selection "$plugin_name"; then
                echo -e "$gene\t$id\t$plugin_name" >> regenie.annotation.txt
            fi
        done
    done
} < reheadered.chr2.vcf.gz.finalAnnot.txt
