#!/bin/bash
#SBATCH --job-name=ULTreeNG
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=100G
#SBATCH --account=r00324
#SBATCH --qos=allocated
#SBATCH --chdir=/N/project/Lennon_Sequences/metabolake
#SBATCH --output=/N/project/Lennon_Sequences/metabolake/ULTreeNG.%j.out

module load raxmlng/1.2.0

raxml-ng --all \
  --msa UL.afa \
  --model GTR+G \
  --threads 4 \
  --seed 12345 \
  --bs-trees autoMRE \
  --outgroup Methanosarcina \
  --prefix UL.ml
