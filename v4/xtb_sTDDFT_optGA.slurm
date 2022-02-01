#!/bin/bash
#SBATCH --output=test.out
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --time=0-05:00:00

# allows cron to run bash
source /ihome/crc/install/lmod/lmod/init/bash

# Load modules
module purge

export XTBHOME=/ihome/ghutchison/geoffh/xtb

inp=${SLURM_JOB_NAME}

# Copy files to SLURM_SCRATCH
cp $SLURM_JOB_NAME $SLURM_SCRATCH

# cd to the SCRATCH space
cd $SLURM_SCRATCH

base=${inp%%.xyz}
echo ${base}
# just getting the NFA file name
molname=${base##*/}
echo ${molname}

# run GFN2-xTB geometry optimization
${XTBHOME}/xtb $SLURM_JOB_NAME -opt > /ihome/ghutchison/blp62/GA/running_GA/opt_GA/GFN2_output/${molname}.out

# copy incremental progress
#cp ${base}*.out ${base}_xtb.out
cp xtbopt.xyz /ihome/ghutchison/blp62/GA/running_GA/opt_GA/opt_xyz/${molname}.xyz
rm -f *restart 

# Load modules
module purge
module load intel openbabel

# create orca input file for sTD-DFT
obabel xtbopt.xyz -O .inp -o orcainp -m -xf /ihome/ghutchison/blp62/GA/running_GA/cam_b3lyp_def2_SVP_pp.txt

module purge
module load openmpi/3.1.4
module load orca/4.2.0


# run sTDDFT, $(which orca) is necessary
$(which orca) xtbopt.inp > ${base}.out

# Copy outputs to submit directory
#cp $SLURM_SCRATCH/*.{inp,prop} $SLURM_SUBMIT_DIR

cp $SLURM_SCRATCH/xtbopt.inp ${base}.inp
cp $SLURM_SCRATCH/xtbopt.prop ${base}.prop
cp ${base}.out /ihome/ghutchison/blp62/GA/running_GA/sequence/sTDDFT_output/${molname}.out
cp ${base}.prop /ihome/ghutchison/blp62/GA/running_GA/sequence/sTDDFT_output/${molname}.prop






