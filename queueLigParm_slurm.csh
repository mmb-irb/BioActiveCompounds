#!/bin/tcsh 

# Usage: queueLigParm.csh <PDB_Code> <pdb> <top> <crd> <type>
if ( ($#argv < 4) || ($#argv > 12) ) then
        echo "USAGE: $0 <pdbCode> <ligCode> <numReps> <numProcsPerRep> [<MDlength>] [<ini_temp>] [<end_temp>] [<boxSize>] [<boxType>] [<replex>] [<ensemble>] [<cluster>]"
        exit
endif

set pdb = $1
set lig = $2
set nreps = $3
set nProcsPerRep = $4
set mdlength = $5
set iniTemp = $6
set endTemp = $7
set boxsize = $8
set boxtype = $9
set replex = $10
set ensemble = $11
set cluster = $12
set days = 5  
set nprocs = 8

set tmpFile = parm${pdb}-${lig}.cluster.csh

echo Queueing parm${pdb}-${lig}...

set dir = `pwd`

echo "#\!/bin/bash " > $tmpFile
echo "#SBATCH -t $days-00:00:00" >> $tmpFile
echo "#SBATCH -J g$pdb-$lig-3" >> $tmpFile
echo "#SBATCH -e g$pdb-$lig-3.e" >> $tmpFile
echo "#SBATCH -o g$pdb-$lig-3.o" >> $tmpFile
echo "#SBATCH -p MPI" >> $tmpFile
echo "#SBATCH -D $dir" >> $tmpFile
echo "#SBATCH --ntasks=$nprocs" >> $tmpFile

#source /etc/profile.d/modules.sh

echo " " >> $tmpFile

#echo "module load CUDA" >> $tmpFile
#echo "module load gaussian/09" >> $tmpFile
echo "module gromacs-plumed-hrex/4.6.7" >> $tmpFile
echo "module unload openmpi/1.8.1" >> $tmpFile
echo "module load openmpi-slurm" >> $tmpFile
echo "module load OpenBabel" >> $tmpFile

cat <<EOT >> $tmpFile

perl /orozco/netapp/hospital/Sanja/Scripts/ligParmHrexSmilesSlurmMPIRUN_newEqNVT_gmx4.pl -pdb $pdb -lig $lig -reps $nreps -procs $nProcsPerRep -len $mdlength -tlow $iniTemp -thigh $endTemp -boxsize $boxsize -boxtype $boxtype -replex $replex -ensemble $ensemble -cluster $cluster 

EOT

cat $tmpFile

sbatch $tmpFile

echo Queued...

