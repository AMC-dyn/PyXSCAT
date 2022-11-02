#!/bin/sh
# Grid Engine options (lines prefixed with #$)
#$ -N lif-ac
#$ -cwd
#$ -l h_rt=48:00:00
#$ -l h_vmem=8G
#$ -pe sharedmem 1
#  These options are:
#  job name: -N
#  use the current working directory: -cwd
#  runtime limit of 5 minutes: -l h_rt
#  memory limit of 2 Gbyte: -l h_vmem
#  parallel environment and no. of slots: -pe

input=lif_ac_m.inp

# Initialise the environment modules
. /etc/profile.d/modules.sh

# Export environment variables

# Run the program
# Defaults are $TMPDIR for tmp scratch and $HOME/wfu for wfu
/exports/applications/apps/community/chem/Molpro/molpros_2012_1_Linux_x86_64_i8/bin/molpro $input
