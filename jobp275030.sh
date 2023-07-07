#!/bin/bash
#SBATCH --time=7-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --job-name=script_p275030
#SBATCH --mem=64G
#SBATCH --mail-user=n.bhushan@rug.nl
#SBATCH --mail-type=BEGIN,END,FAIL

pwd
module load R/3.4.2-foss-2016a-X11-20160819
export LD_LIBRARY_PATH=$HOME/v8/usr/lib64/:$LD_LIBRARY_PATH
cd $HOME/CausalDiscoverySimulation/
Rscript cdSimulate.R
