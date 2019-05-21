#!/bin/tcsh
#SBATCH -A gsd-hpcs
#SBATCH -q debug
###SBATCH -A sena
###SBATCH -q fge
#SBATCH -o gtest

## NOTE: 20 FGE cores, 8 GPUs
#SBATCH --ntasks-per-node=16
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --time=00:10:00

setenv USE_GPUS 0

cd ${SLURM_SUBMIT_DIR}

module load intel/18.1.163
module load impi
module load cuda


#setenv LD_LIBRARY_PATH /apps/cuda/cuda-8.0/lib64

setenv I_MPI_PIN         disable

# Number of CPU cores
setenv CPU_CORES_PER_RANK 1

setenv CUDA_DEVICE_MAX_CONNECTIONS 12
setenv CUDA_COPY_SPLIT_THRESHOLD_MB 1

# Using impi, must set the following to
# see the filesystem:
setenv I_MPI_EXTRA_FILESYSTEM on
setenv I_MPI_EXTRA_FILESYSTEM_LIST lustre:panfs

# KMP_AFFINITY controls affinity with impi:
setenv KMP_AFFINITY scatter

setenv NP   $SLURM_NTASKS
setenv NPPN $SLURM_NTASKS_PER_NODE
#setenv OMP_NUM_THREADS $SLURM_CPUS_PER_TASK

if ( $NPPN == 2 ) then
  setenv OMP_NUM_THREADS 10
else if ( $NPPN == 4 ) then
  setenv OMP_NUM_THREADS 4
else if ( $NPPN == 8 ) then
  setenv OMP_NUM_THREADS 2
else
##@ OMP_NUM_THREADS = 20 / $NPPN
setenv OMP_NUM_THREADS 20
endif



# Following used for use on FGA system due to
# use different IB cards than on reg theia nodes:
setenv I_MPI_FABRICS shm:ofa


# Set GPU ids, if necessary:
if ( $USE_GPUS == 1 ) then

  if ( $NP <= 4 ) then
    setenv CUDA_VISIBLE_DEVICES "0,1,2,3"
  else
    setenv CUDA_VISIBLE_DEVICES "0,1,2,3,4,5,6,7"
  endif


  if ( $NPPN > 8 ) then
    echo "Must have 1 MPI task bound to 1 GPU!"
  endif

endif

echo "NP=" $NP
echo "NPPN=" $NPPN
echo "OMP_NUM_THREADS=" $OMP_NUM_THREADS


module list
scontrol show hostname $SLURM_NODELIST#  exit

#srun -np $NP  affin_theia.sh
srun -n $NP  ../../bin/gtest_burgers
