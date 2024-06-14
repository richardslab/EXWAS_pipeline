#!/bin/bash
#SBATCH --job-name=exwas_annotation  # Specify the job name
#SBATCH --nodes=1         # Specify the number of nodes
#SBATCH --ntasks=1        # Specify the number of tasks (typically used for MPI jobs)
#SBATCH --cpus-per-task=15 # Specify the number of CPUs per task
#SBATCH --mem=64G          # Specify the memory limit per node
#SBATCH --time=24:00:00   # Specify the maximum walltime for the job
#SBATCH --output=/home/richards/ethan.kreuzer/projects/richards/ethan.kreuzer/EXWAS/job-%A-%a.log

INPUT_VCF=$1

ANCESTRY=$2


module load tabix
module load bcftools
module load apptainer/1.1.8
sif=/scratch/richards/guillaume.butler-laporte/WGS/COVID19_GenOMICC_AVT_analysis/vep_v105.sif

cd /home/richards/ethan.kreuzer/projects/richards/ethan.kreuzer/EXWAS


# Normalize, left-align and drop genotypes of the input VCF
bcftools view --drop-genotypes "${INPUT_VCF}" -Ou | bcftools norm -m -any --check-ref w -f /scratch/richards/ethan.kreuzer/vep/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz -Ou | bcftools annotate --set-id '%CHROM:%POS:%REF:%FIRST_ALT' -Oz > "${INPUT_VCF}".set_id
tabix -p vcf "${INPUT_VCF}".set_id


#Subset to specified Ancestry
bcftools view -S /home/richards/ethan.kreuzer/scratch/1000G_ancestries/"${ANCESTRY}"_sample_ids.txt "${INPUT_VCF}".set_id -Oz > "${INPUT_VCF}".set_id."${ANCESTRY}"
tabix -p vcf "${INPUT_VCF}".set_id."${ANCESTRY}"


chromosome=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X)
chr="${chromosome[SLURM_ARRAY_TASK_ID]}"

apptainer run --bind ${PWD}:${PWD} ${sif} vep -i "${INPUT_VCF}".set_id."${ANCESTRY}" \
         --assembly GRCh38 \
         --vcf \
         --format vcf \
         --cache \
         --dir_cache /scratch/richards/ethan.kreuzer/vep \
         -o ${inputvcf}.finalAnnot_loftee.vcf \
         --plugin LoF,loftee_path:/opt/micromamba/share/ensembl-vep-105.0-1,human_ancestor_fa:vep_data/human_ancestor.fa.gz,conservation_file:vep_data/loftee.sql,gerp_bigwig:vep_data/gerp_conservation_scores.homo_sapiens.GRCh38.bw \
         --plugin CADD,/scratch/richards/ethan.kreuzer/vep/whole_genome_SNVs.tsv.gz,/scratch/richards/ethan.kreuzer/vep/gnomad.genomes.r3.0.indel.tsv.gz \
         --plugin dbNSFP,/scratch/richards/ethan.kreuzer/vep/dbNSFP4.8a_grch38.gz,Ensembl_transcriptid,Uniprot_acc,VEP_canonical,LRT_pred,SIFT_pred,MutationTaster_pred,Polyphen2_HDIV_pred,Polyphen2_HVAR_pred \
	 --everything \
         --force_overwrite \
         --offline \
         --fork 1 \
         --quiet

#bgzip and index output
bgzip -f ${inputvcf}.finalAnnot_loftee.vcf
tabix -f ${inputvcf}.finalAnnot_loftee.vcf

