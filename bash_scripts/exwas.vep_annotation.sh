#!/bin/bash
#SBATCH --job-name=exwas_vep_annotation  # Specify the job name
#SBATCH --nodes=1         # Specify the number of nodes
#SBATCH --ntasks=1        # Specify the number of tasks (typically used for MPI jobs)
#SBATCH --cpus-per-task=15 # Specify the number of CPUs per task
#SBATCH --mem=64G          # Specify the memory limit per node
#SBATCH --time=14:00:00   # Specify the maximum walltime for the job
#SBATCH --output=/home/richards/ethan.kreuzer/projects/richards/ethan.kreuzer/EXWAS/job-%A-%a.log

INPUT_VCF=$1

cd /scratch/richards/ethan.kreuzer
module load tabix
module load bcftools
module load apptainer/1.1.8
sif=/scratch/richards/guillaume.butler-laporte/WGS/COVID19_GenOMICC_AVT_analysis/vep_v105.sif

# Normalize, left-align and subset ancestry of the input VCF

bcftools view --drop-genotypes -Ou "${INPUT_VCF}" | bcftools norm -m -any --check-ref w -f /scratch/richards/ethan.kreuzer/vep/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna -Ou | bcftools annotate --set-id '%CHROM:%POS:%REF:%FIRST_ALT' -Oz -o "${INPUT_VCF}.set_id.no_genotypes"

tabix -p vcf "${INPUT_VCF}".set_id.no_genotypes


apptainer run --bind ${PWD}:${PWD} ${sif} vep -i "${INPUT_VCF}".set_id.no_genotypes \
         --assembly GRCh38 \
         --format vcf \
         --cache \
	 -o "${INPUT_VCF}".finalAnnot.filter_vep.txt \
         --dir_cache /scratch/richards/ethan.kreuzer/vep \
         --plugin LoF,loftee_path:/opt/micromamba/share/ensembl-vep-105.0-1,human_ancestor_fa:/scratch/richards/ethan.kreuzer/vep/human_ancestor.fa.gz,conservation_file:/scratch/richards/ethan.kreuzer/vep/loftee.sql,gerp_bigwig:/scratch/richards/ethan.kreuzer/vep/gerp_conservation_scores.homo_sapiens.GRCh38.bw \
         --plugin CADD,/scratch/richards/ethan.kreuzer/vep/whole_genome_SNVs.tsv.gz,/scratch/richards/ethan.kreuzer/vep/gnomad.genomes.r3.0.indel.tsv.gz \
         --plugin dbNSFP,/scratch/richards/ethan.kreuzer/vep/dbNSFP4.8a_grch38.gz,gnomAD_genomes_AFR_AF,gnomAD_genomes_AMI_AF,gnomAD_genomes_AMR_AF,gnomAD_genomes_ASJ_AF,gnomAD_genomes_EAS_AF,gnomAD_genomes_FIN_AF,gnomAD_genomes_MID_AF,gnomAD_genomes_NFE_AF,gnomAD_genomes_SAS_AF,EVE_Class25_pred,VEP_canonical,LRT_pred,SIFT_pred,SIFT4G_pred,MutationTaster_pred,Polyphen2_HDIV_pred,Polyphen2_HVAR_pred,AlphaMissense_pred \
         --force_overwrite \
         --offline \
	 --symbol \
	 --coding_only \
	 --no_stats \
         --fork 1 \
         --quiet 



