#!/bin/bash -l

#$ -P frgeeeph       # Specify the SCC project name you want to use
#$ -N fputAGPNormAvg        # Give job a name
#$ -l cpu_type=Gold-6242
#$ -l h_rt=72:00:00   # Specify the hard time limit for the job
#$ -j y               # Merge the error and output streams into a single file
#$ -m e               #Send email after job is done

./_fputAGPNormAvg.out 64 1 1 0.$SGE_TASK_ID 1 $SGE_TASK_ID
