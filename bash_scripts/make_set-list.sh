#!/bin/bash
#SBATCH --job-name=exwas_set-list  # Specify the job name
#SBATCH --nodes=1         # Specify the number of nodes
#SBATCH --ntasks=1        # Specify the number of tasks (typically used for MPI jobs)
#SBATCH --cpus-per-task=15 # Specify the number of CPUs per task
#SBATCH --mem=64G          # Specify the memory limit per node
#SBATCH --time=12:00:00   # Specify the maximum walltime for the job
#SBATCH --output=/scratch/richards/ethan.kreuzer/job-%A-%a.log

INPUT_VCF=$1

cd /scratch/richards/ethan.kreuzer

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

awk '!/^##/ && !/^#/' "${INPUT_VCF}".finalAnnot.filter_vep.txt > "${INPUT_VCF}".finalAnnot.filter_vep.txt.no_hash

awk '{ print $4, $1 }' "${INPUT_VCF}".finalAnnot.filter_vep.txt.no_hash | sort -u -k 1.5,1 -k 2.6,2 > canon.chr2.var.gene.txt

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

sed -E 's/([^ ]+ [^ ]+ )[0-9]+chr/\1chr/' "$output_file" > regenie.set-list.txt
rm "$output_file"
rm canon.chr2.gene.txt canon.chr2.var.gene.txt "${INPUT_VCF}".finalAnnot.filter_vep.txt.no_hash
