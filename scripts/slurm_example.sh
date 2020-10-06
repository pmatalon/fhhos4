#!/bin/bash

# 1. Copy/paste this example file.
# 2. Edit the copy according to the platform.
# 3. Run the script with a command like the following:
#    sbatch --partition=prod --time=02:00:00 slurm.sh -geo square -n 128


# Options SBATCH
#SBATCH --job-name=DGHHO
#SBATCH --nodes=1
#SBATCH --cpus-per-task=36

#SBATCH --mail-type=ALL
#SBATCH --mail-user=matalon@cerfacs.fr

# Load required modules
#module purge

# Concatenate all the arguments to build the output filename
arguments="$*"
outfilename=${arguments//[[:blank:]]/}
outfilename=${outfilename:4}

# Call the program
/scratch/algo/matalon/dghho "$@" > /home/algo/matalon/dghho/out/${outfilename}.log