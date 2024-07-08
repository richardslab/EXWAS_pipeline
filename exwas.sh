#!/bin/bash
#SBATCH --job-name=exwas_pipeline  # Specify the job name
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

#IN the final pipeline I will combine the drop genotypes command with this above and assume the user provided the ancestry selected plink files themselves

#/scratch/richards/yiheng.chen/Plink1.9/plink --vcf "${INPUT_VCF}".set_id --make-bed --out "${INPUT_VCF}".set_id.plk

#/scratch/richards/yiheng.chen/Plink1.9/plink --bfile "${INPUT_VCF}".set_id.plk  --hwe 1E-15 midp  --maf 0.01  --geno 0.1  --indep-pairwise 50 5 0.05  --out pruned_variants.txt


chromosome=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X)
chr="${chromosome[SLURM_ARRAY_TASK_ID]}"

apptainer run --bind ${PWD}:${PWD} ${sif} vep -i "${INPUT_VCF}".set_id.no_genotypes \
         --assembly GRCh38 \
         --format vcf \
         --cache \
	 -o "${INPUT_VCF}".finalAnnot.txt \
         --dir_cache /scratch/richards/ethan.kreuzer/vep \
         --plugin LoF,loftee_path:/opt/micromamba/share/ensembl-vep-105.0-1,human_ancestor_fa:/scratch/richards/ethan.kreuzer/vep/human_ancestor.fa.gz,conservation_file:/scratch/richards/ethan.kreuzer/vep/loftee.sql,gerp_bigwig:/scratch/richards/ethan.kreuzer/vep/gerp_conservation_scores.homo_sapiens.GRCh38.bw \
         --plugin CADD,/scratch/richards/ethan.kreuzer/vep/whole_genome_SNVs.tsv.gz,/scratch/richards/ethan.kreuzer/vep/gnomad.genomes.r3.0.indel.tsv.gz \
         --plugin dbNSFP,/scratch/richards/ethan.kreuzer/vep/dbNSFP4.8a_grch38.gz,gnomAD_genomes_AFR_AF,gnomAD_genomes_AMI_AF,gnomAD_genomes_AMR_AF,gnomAD_genomes_ASJ_AF,gnomAD_genomes_EAS_AF,gnomAD_genomes_FIN_AF,gnomAD_genomes_MID_AF,gnomAD_genomes_NFE_AF,gnomAD_genomes_SAS_AF,EVE_Class25_pred,VEP_canonical,LRT_pred,SIFT_pred,SIFT4G_pred,MutationTaster_pred,Polyphen2_HDIV_pred,Polyphen2_HVAR_pred,AlphaMissense_pred \
         --force_overwrite \
         --offline \
         --fork 1 \ #DO NOT INCREASE THIS OR LOFTEE WILL FAIL
         --quiet 

/home/richards/ethan.kreuzer/projects/richards/ethan.kreuzer/EXWAS/ensembl-vep/filter_vep --format tab -i "${INPUT_VCF}".finalAnnot.txt --filter "Gene" | /home/richards/ethan.kreuzer/projects/richards/ethan.kreuzer/EXWAS/ensembl-vep/filter_vep --filter "LRT_pred match .*D.* or CADD_PHRED >= ${CADD_threshold} or AlphaMissense_pred match .*P.* or EVE_Class25_pred match .*P.* or LoF is HC or LoF is LC or MutationTaster_pred match .*D.* or MutationTaster_pred match .*A.* or Polyphen2_HDIV_pred match .*D.* or Polyphen2_HVAR_pred match .*D.* or SIFT4G_pred match .*D.* or SIFT_pred match .*P.*" -o "${INPUT_VCF}".finalAnnot.filter_vep.txt


# The following code is to create the annotation file needed for regenie according to the user specs

declare -A mask_selection

declare -A vep_consequences_map

while IFS= read -r mask; do
    mask_selection["$mask"]=1
done < "${mask_file}"


vep_consequences_str="transcript_ablation,splice_acceptor_variant,splice_donor_variant,stop_gained,frameshift_variant,stop_lost,start_lost,transcript_amplification,feature_elongation,feature_truncation,inframe_insertion,inframe_deletion,missense_variant,protein_altering_variant,splice_donor_5th_base_variant,splice_region_variant,splice_donor_region_variant,splice_polypyrimidine_tract_variant,incomplete_terminal_codon_variant,start_retained_variant,stop_retained_variant,synonymous_variant,coding_sequence_variant,mature_miRNA_variant,5_prime_UTR_variant,3_prime_UTR_variant,non_coding_transcript_exon_variant,intron_variant,NMD_transcript_variant,non_coding_transcript_variant,coding_transcript_variant,upstream_gene_variant,downstream_gene_variant,TFBS_ablation,TFBS_amplification,TF_binding_site_variant,regulatory_region_ablation,regulatory_region_amplification,regulatory_region_variant,intergenic_variant,sequence_variant"

# Convert the comma-separated string into an array
IFS=',' read -r -a vep_consequences <<< "$vep_consequences_str"

for vep_consequence in "${vep_consequences[@]}"; do
    vep_consequences_map["$vep_consequence"]=1
done

# Function to check if a plugin value is deleterious

is_deleterious() {
    local plugin=$1
    local value=$2

    case $plugin in
        LoF)
            if [[ "$LOF_threshold" == "HC" ]]; then
                [[ "$value" == "HC" ]] && return 0
            else
                return 0
            fi
            ;;
        CADD_PHRED)
            (( $(awk "BEGIN {print ($value >= ${CADD_threshold})}") )) && return 0
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

in_mask_selection() {
    local item=$1
    [[ -n "${mask_selection[$item]}" ]] && return 0
    return 1
}

is_combined_deleterious() {
    local combined=$1
    local consequences=$2
    local extra=$3
    IFS='&&' read -r -a criteria <<< "$combined"
    IFS=';' read -r -a plugins <<< "$extra"
    IFS=',' read -r -a consequence_array <<< "$consequences"
    
    for criterion in "${criteria[@]}"; do
    	#echo "Criterion: $criterion"
	local found=1    
        
	if [[ -z "$criterion" ]]; then
            continue
        fi


	# Check if criterion is a VEP consequence
	if [[ ${vep_consequences_map["$criterion"]+_} ]]; then
		for single_consequence in "${consequence_array[@]}"; do
                
			if [[ "$single_consequence" == "$criterion" ]]; then
              			found=0
				break
                	fi
     		done                
		
		if [[ found -eq 1 ]]; then
			
			#echo "Criterion $criterion not found in consequences"
                	return 1
		fi
	fi

	if [[ found -eq 0 ]]; then

		continue
	fi

	found=1

	for plugin in "${plugins[@]}"; do
            
	    name="${plugin%%=*}"
            value="${plugin#*=}"            
	    
	    if [[ "$name" == "$criterion" ]]; then

		if is_deleterious "$name" "$value"; then
			#echo "Plugin $name is deleterious"
			found=0
			break
		fi
            fi
        done

	if [[ found -eq 1 ]]; then
		#echo "Criterion $criterion not found as deleterious plugin"
		return 1
	fi
    done
    return 0
}


# Process the annotated vcf file
{
    read -r header
    output=""

    # Process each line
    while IFS=$'\t' read -r id location allele gene feature feature_type consequence cDNA_position CDS_position protein_position amino_acids codons existing_variation extra; do
        #echo $id
	#echo $consequence
        # Split the consequences and print each one if it is in the mask selection
        IFS=',' read -r -a consequences_array <<< "$consequence"
        for single_consequence in "${consequences_array[@]}"; do
            if in_mask_selection "$single_consequence"; then
            	
		output+="$gene\t$id\t$single_consequence\n"

	    fi
        done

        # Iterate over each plugin in Extra field
        IFS=';' read -r -a plugins <<< "$extra"
        for plugin in "${plugins[@]}"; do
            
	    plugin_name="${plugin%%=*}"
	    plugin_value="${plugin#*=}"

	    #echo $plugin_name
	    #echo $plugin_value
            # Check if the plugin value is deleterious and if the plugin name is in the mask selection
            if is_deleterious "$plugin_name" "$plugin_value" && in_mask_selection "$plugin_name"; then
            
	    	output+="$gene\t$id\t$plugin_name\n"

            fi
        done

        #Check for combined criteria in mask selection
        for mask in "${!mask_selection[@]}"; do
            if [[ $mask == *&&* ]]; then
	    	#echo $mask
                if is_combined_deleterious "$mask" "$consequence" "$extra"; then
                    
		    output+="$gene\t$id\t$mask\n"
                fi
            fi
        done
    done
} < "${INPUT_VCF}".finalAnnot.filter_vep.txt



echo -e "$output" > regenie.annotation.txt

# The next section of code is to create the set-list file for regenie

# Define the transpose function using awk
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

awk '{ print $4, $1 }' "${INPUT_VCF}".finalAnnot.filter_vep.txt | sort -u -k 1.5,1 -k 2.6,2 > canon.chr2.var.gene.txt

awk '{if(!seen[$1]++) { match($2, /^chr([0-9]+)/, arr); print $1, arr[1], ++line }}' canon.chr2.var.gene.txt > canon.chr2.gene.txt

# Read genes into an array
genes=($(awk '{print $1}' canon.chr2.gene.txt | tr '\n' ' '))

# Create the output file
output_file="regenie.set.list.txt.tmp"
> "$output_file"

# Iterate through each gene

for gene in "${genes[@]}"; do
  if [ -n "$gene" ]; then
  	awk -v a="$gene" '$1 == a {print $2}' canon.chr2.var.gene.txt | uniq | transpose | tr -s '[:blank:]' "," | \
	awk 'BEGIN{FS=",";OFS=","} { $1=$1; print $0 }' | \
	awk -F, -v gene="$gene" '{split($1, a, ":"); $1 = gene OFS a[1] OFS a[2] OFS $0; gsub(/^[ \t]+/, "", $1); print $1}' >> "$output_file"
  fi
done

sed -E 's/([^ ]+ [^ ]+ )[0-9]+chr/\1chr/' "$output_file" > regenie.set.list.txt


