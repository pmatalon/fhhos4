#!/bin/bash
#PBS -S /bin/bash

# 1. Copy/paste this example file.
# 2. Edit the copy according to the platform.
# 3. Run the script with a command like the following:
#    qsub -q <queue> -l ppn=36 -l walltime=02:00:00 torque.sh -geo square -n 128


# Job name
#PBS -N fhhos4

# Nodes and number of processsors by node
#PBS -l nodes=1

# Sends email if the job is (a) aborted, when it (b) begins, and when it (e) ends
#PBS -m abe
#PBS -M pierre.matalon@polimi.it

# Joins standard output and standard error
#PBS -j oe

# Load required modules
module load gcc-glibc/9

# Start the job in the current directory
cd ${PBS_O_WORKDIR}


# Concatenate all the arguments to build the output filename
arguments="$*"
outfilename=${arguments//[[:blank:]]/}
outfilename=${outfilename:4}

# Call the program
/u/matalon/fhhos4/build/bin/fhhos4 "$@" > /u/matalon/fhhos4/out/${outfilename}.log 2>&1