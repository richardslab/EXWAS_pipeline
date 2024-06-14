#!/bin/bash
#SBATCH --job-name=exwas_annotation  # Specify the job name
#SBATCH --nodes=1         # Specify the number of nodes
#SBATCH --ntasks=1        # Specify the number of tasks (typically used for MPI jobs)
#SBATCH --cpus-per-task=15 # Specify the number of CPUs per task
#SBATCH --mem=64G          # Specify the memory limit per node
#SBATCH --time=24:00:00   # Specify the maximum walltime for the job
#SBATCH --output=/home/richards/ethan.kreuzer/projects/richards/ethan.kreuzer/EXWAS/job-%A-%a.log

INPUT_VCF=$1

#PHENOTYPE=$2

ANCESTRY=$2

#MASKS=$4


#assume for now that VEP and LOFTEE work

module load tabix

module load bcftools

cd /home/richards/ethan.kreuzer/projects/richards/ethan.kreuzer/EXWAS

# Normalize the input VCF
bcftools view --drop-genotypes "${INPUT_VCF}" -Ou | bcftools norm -m -any --check-ref w -f /scratch/richards/ethan.kreuzer/vep/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz -Ou | bcftools annotate --set-id '%CHROM:%POS:%REF:%FIRST_ALT' -Oz > "${INPUT_VCF}".set_id

# Index the output VCF
tabix -p vcf "${INPUT_VCF}".set_id


#Subset to specified Ancestry

bcftools view -S /home/richards/ethan.kreuzer/scratch/1000G_ancestries/"${ANCESTRY}"_sample_ids.txt "${INPUT_VCF}".set_id -Oz > "${INPUT_VCF}".set_id."${ANCESTRY}"

tabix -p vcf "${INPUT_VCF}".set_id."${ANCESTRY}"


/home/richards/ethan.kreuzer/projects/richards/ethan.kreuzer/EXWAS/ensembl-vep/vep -i "${INPUT_VCF}".set_id."${ANCESTRY}" \
        --plugin LoFtool \
        -everything \
        --assembly GRCh38 \
        --offline \
        --dir_cache /home/richards/ethan.kreuzer/.vep \
        --cache -o "${INPUT_VCF}".set_id."${ANCESTRY}".finalAnnot.txt ;


