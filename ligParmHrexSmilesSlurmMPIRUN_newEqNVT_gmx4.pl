#!/usr/bin/perl -w

# REQUIRED ENVIRONMENT:
# module load OpenBabel
# module load gromacs
# module load amber-14

use Getopt::Long;
use POSIX;
use strict;

# Input Parameters.
my $long=@ARGV;
if ($long == 0) {
	&usage();
	exit;
}

#######################  Check input parameters  #######################

my %limits;
$limits{'maxLength'} = 100;	# Max MD time: 100ns
$limits{'maxReps'} = 50;	# Max replicas: 50
$limits{'minTemp'} = 100;	# Min temperature: 100K
#$limits{'maxTemp'} = 500;	# Max temperature: 500K
$limits{'maxTemp'} = 1000;	# Max temperature: 1000K

$limits{'maxClInc'} = 0.1;      # Max Increment in iterative clustering: 1 Angstrom (0.1 nm)
$limits{'maxClModels'} = 50;    # Max number of cluster conformations: 50

my %input;

&init(\%input,%limits);

#######################  Auxiliar variables  #######################

# Home dir of the project (will change)
my $homeDir = "/orozco/netapp/hospital/Sanja/Scripts";

# Work dir of the project (might change)
my $workDir = `pwd`; 
chomp($workDir);

if ($input{'name'}){
	my $projectName = $input{'name'};
	mkdir("$projectName") if(! -s "$projectName");
	$workDir .="/$projectName";
}

# Path to the PDB mirror
#my $mirror = "/orozco/raid7/mirror/pdb/data/structures/all/pdb";
#my $mirror = "/orozco/projects/pdb/pdb/data/structures/all/pdb";
my $mirror = ".";

# ACPYPE python script
# https://code.google.com/p/acpype/
my $acpype = "$homeDir/ACPYPE/acpype.py";

# Gromacs minimization config file
my $gmx_em_file = "minimizeBoxSolvent.mdp";

# Gromacs REMD equilibration config files
my $gmx_eq_file = "equilibrateBoxSolvent.gmx4.mdp";
my @gmx_eq_files;
$gmx_eq_files[1] = "equilibrateBoxSolvent_step1.NVT.gmx4.mdp";
$gmx_eq_files[2] = "equilibrateBoxSolvent_step2.NVT.gmx4.mdp";

# Gromacs REMD production config files
#my $gmx_run_file = "runMDBoxSolventNPT.mdp";
#my $gmx_run_file = "runMDBoxSolvent$input{'ensemble'}.mdp";
my $gmx_run_file = "runMDBoxSolvent$input{'ensemble'}.gmx4.mdp";
#my $gmx_run_file = "runMDBoxSolvent$input{'ensemble'}.shake.gmx4.mdp";

# Molecule charge
my $charge = '';

#######################  Input Parameters  #######################

my $pdbfile;

if ($input{'ligFile'} && !$input{'lig'}) {

	my $res = '';
	open LIG,"$input{'ligFile'}";
	while (<LIG>){
		chomp;
		next if !(/^ATOM|^HETATM/);
		# HETATM   14  N14 E20 A2001
		$input{'lig'} = substr($_,17,3);
		last;
	}
	close LIG;

	$pdbfile = $input{'ligFile'};
	`cp $pdbfile $workDir`;
}

if ($input{'pdbFile'}) {
	
	$pdbfile = $input{'pdbFile'};
	`cp $pdbfile $workDir`;
}

if ($input{'smilesFile'}){

	# Converting SMILES file to PDB
	`babel -ismi $input{'smilesFile'} -opdb $workDir/LIG.pdb --gen3d`;

	$pdbfile = "LIG.pdb";
	$input{'lig'} = "LIG";
	`cp $pdbfile $workDir`;
} 

chdir("$workDir");

if ($input{'pdb'}) {

	# Pdb file with mirror format (e.g. pdb1a32.ent.gz)
	$pdbfile = "pdb".lc($input{'pdb'}).".ent";

	# Getting PDB file from MMB PDB API (through a REST call)
	`wget http://mmb.irbbarcelona.org/api/pdb/$input{'pdb'} -O $pdbfile -o REST.log`;
}

if ($input{'smiles'}){

	# SMILES string to file
	open SMI,">smiles.tmp";
	print SMI $input{'smiles'};
	close SMI;

	# Converting SMILES file to PDB
	`babel -ismi smiles.tmp -opdb LIG.pdb --gen3d`;

	$pdbfile = "LIG.pdb";
	$input{'lig'} = "LIG";
} 

if (!$input{'pdb'}){
	#my $date = localtime();
	#my ($day,$mon,$num,$hour,$year) = split ' ',$date;
	#$input{'pdb'} = "$num$mon$year";
	
	$input{'pdb'} = "UNKN";
}

if (!$input{'lig'}) {
	
	#my $cont = 1;
	#my $codeLig = sprintf("%03d",$cont);
	#my $proposDir = "$workDir/$input{'pdb'}-$codeLig";
	#while (-s "$proposDir"){
	#	$cont++;
	#	$codeLig = sprintf("%03d",$cont);
	#	$proposDir = "$workDir/$input{'pdb'}-$codeLig";
	#}
	#$input{'lig'} = $codeLig;
	#print "CodeLig: $codeLig \n";
	$input{'lig'} = "LIG";
}

# Creating Working dir for the particular input
my $workDirProject = "$workDir/$input{'pdb'}-$input{'lig'}";

if($input{'force'}){
	`rm -r $workDirProject` if(-s "$workDirProject");
}

mkdir("$workDirProject") if(! -s "$workDirProject");

# Getting rid of path folders and storing the name of the PDB file to work with it.
if ($pdbfile =~ /\//){
	$pdbfile=~/\/([^\/]+)$/;
	$pdbfile = $1;
}
`cp $pdbfile $workDirProject` if (! -s "$workDirProject/$pdbfile");

chdir("$workDirProject");

#######################  Ligand Extraction  #######################

if (! $input{'onlyREMD'} and ! $input{'onlyCluster'} and ! $input{'onlyClusterMin'} and ! $input{'onlyGaussian'}) {

	print "We are going to try to parameterize ligand $input{'lig'} from pdb $input{'pdb'} running a REMD with $input{'numReps'} replicas of $input{'length'} ns from $input{'ini_temp'} to $input{'end_temp'} K\n\n";

	print "Extracting ligand $input{'lig'} from pdb $input{'pdb'}... \n";

	# Converting lig code to uppercase (as it is in the pdb file)
	$input{'lig'} = uc($input{'lig'});

	# Extracting ligand from pdb file
	my $pdbLig = extractLig($input{'lig'},$mirror,$pdbfile);

	# DEBUG: Printing resulting Lig pdb
	print "LIG:\n$pdbLig\n" if ($input{'verbose'});

	# Writing Ligand PDB
	open LIG, ">lig$input{'lig'}.pdb";
	print LIG $pdbLig;
	close LIG;

	print "Done!\n";

	exit if ($input{'onlyLig'});

	#######################  Ligand Parameterization  #######################

	print "Computing parameters for ligand $input{'lig'}...  \n";

	# Running OpenBabel to add protons at a particular pH
	print "Running OpenBabel (Addding protons and minimizing the structure)...\n";
	my $ph = 7.4;
	`babel lig$input{'lig'}.pdb lig$input{'lig'}.H.pdb -h -p$ph`;

	# Computing ligand charge just summing up number of protons added by OpenBabel.
	#my $charge = `grep "H  " lig$input{'lig'}.H.pdb | wc -l`;
	$charge = getCharge("lig$input{'lig'}.H.pdb");
	print "Ligand Charge: $charge\n" if ($input{'verbose'});

	# Minimizing the energy of the ligand with protons a little bit
	`obminimize -sd -c 1e-10 -ipdb lig$input{'lig'}.H.pdb -opdb lig$input{'lig'}.H.min.pdb >& lig$input{'lig'}.H.min.log`;
	`egrep "ATOM|HETATM|CONECT|MASTER|COMPND|AUTHOR" lig$input{'lig'}.H.min.log > lig$input{'lig'}.H.min.pdb`;

        # Generating minimized pdb for ligand without hydrogens (to be used in future RMSd calculations)
        &removeH("lig$input{'lig'}.H.min.pdb","lig$input{'lig'}.min.pdb");

	# Running ACPYPE
	print "Running ACPYPE...\n";
	# DEBUG: Printing ACPYPE cmd line
	print "python $acpype -i lig$input{'lig'}.H.min.pdb -b $input{'pdb'}-$input{'lig'} -n $charge\n" if ($input{'verbose'});
	my $acpypeLOG = `python $acpype -i lig$input{'lig'}.H.min.pdb -b $input{'pdb'}-$input{'lig'} -n $charge`;

	# DEBUG: Printing ACPYPE log
	print "$acpypeLOG\n" if ($input{'verbose'});

	# Check output
	my $outAC = "$input{'pdb'}-$input{'lig'}.acpype/sqm.out";
	if (! -s "$outAC"){
		die ("ERROR in Parameterization step (Antechamber -- ACPYPE). Please check log files.\n");
	}
	else{
		my $complete = `grep "Calculation Completed" $outAC`;
		if (! $complete) {
			die ("ERROR in Parameterization step (Antechamber -- ACPYPE). Please check log files.\n");
		}
	}

	print "Done!\n";

	exit if ($input{'onlyParam'});
}

#########################  Gromacs REMD  #########################

if (! $input{'onlyCluster'} and ! $input{'onlyClusterMin'} and ! $input{'onlyGaussian'}){

	# Gromacs REMD preparation
	mkdir("$workDirProject/GMX") if (! -e "$workDirProject/GMX");
	chdir("$workDirProject/GMX");

	print "\nPreparing Gromacs files generated with ACPYPE... ";

	my $eq_top_file = "$input{'pdb'}-$input{'lig'}_GMX.mod.top";
	my $eq_gro_file = "$input{'pdb'}-$input{'lig'}_wat.gro";

	&prepare_gmx_files(\%input,$eq_top_file,$eq_gro_file);

	print "\nDone!\n";

	print "Now equilibrating all replicas...\n";
	my $eqSteps = 1;
	for (my $eqStep = 1; $eqStep <= $eqSteps; $eqStep++){
		&equilibrate_REMD(\%input,$eqStep,$eq_top_file,$eq_gro_file);		
	}

	print "\nDone!\n";

	exit if ($input{'onlyEquil'});

	print "Now running production REMD...\n";
	my $check = &production_REMD(\%input,$eqSteps,$eq_top_file,$eq_gro_file);

	if ($check) {
		print "Done!\n";
	}
	else{
		print "ERROR: Final Trajectory is smaller than expected. Please check output and log files. Exiting...\n";
		exit;
	}

	exit if ($input{'onlyREMD'});

}

#########################  Gromacs Clustering  #########################

chdir("$workDirProject/GMX"); # If just computing clustering...

if (! $input{'onlyClusterMin'} and ! $input{'onlyGaussian'}){

	print "\nComputing cluster in REMD...  ";

        my $cl_cutoff = $input{'cluster'};
        my $cl_nmodels = $input{'clModels'};
        my $cl_percentage = $input{'clPercentage'};
        my $cl_inc = $input{'clIncrement'};

        # Iterative clustering: increasing clustering cutoff by $cl_inc until number of models <= $cl_nmodels
        my $done = 0;
        while (!$done) {

		# Computing cluster in traj with lowest temperature
		my $options = "";
		$options .= "-method $input{'clMethod'} "; 
		$options .= "-cutoff $cl_cutoff "; 
		$options .= "-dist clusters.$input{'clMethod'}.$cl_cutoff.rmsdDist ";
		$options .= "-sz clusters.$input{'clMethod'}.$cl_cutoff.sizes ";

		print "srun -n 1 g_cluster_mpi -s sim1/topol.tpr -f sim1/$input{'pdb'}-$input{'lig'}.imaged.rot.xtc -cl clusters.pdb -dista $options\n";
		`echo "1 1" | srun -n 1 g_cluster_mpi -s sim1/topol.tpr -f sim1/$input{'pdb'}-$input{'lig'}.imaged.rot.xtc -cl clusters.pdb -dista $options`;

                if (! -s "clusters.pdb"){
                        print "Ups, clusters.pdb has not been generated by gmx cluster. Please check it and try again!\n";
                        exit;
                }

                my $nmodels = `grep MODEL clusters.pdb | wc -l`;

                my $cl_ok = &checkClusterPercentage("cluster.log",$cl_nmodels,$cl_percentage);

                #if ($nmodels > $cl_nmodels) {
                if (!$cl_ok){
                        $cl_cutoff = $cl_cutoff + $cl_inc;
                        unlink "clusters.pdb";
                        unlink "rmsd-clust.xpm";
                        $done = 0;
                }
                else {
                        $done = 1;
                }
	}
}

if (! $input{'onlyGaussian'}){

	# Getting number of clusters and pdb files corresponding to them 

	my $titles = `grep ^TITLE clusters.pdb`;

	mkdir("$workDirProject/CLUSTERS/") if (! -e "$workDirProject/CLUSTERS");
	chdir("$workDirProject/CLUSTERS");

	my $i = 1;
	foreach my $title (split '\n', $titles){

		mkdir("$workDirProject/CLUSTERS/CL$i") if (! -e "$workDirProject/CLUSTERS/CL$i");
		chdir("$workDirProject/CLUSTERS/CL$i");

		# TITLE     frame t= 6598.000
		my ($tag,$frame,$t,$dump) = split ' ',$title;
		print "Dump: $dump\n";

		print " srun -n 1 trjconv_mpi -f $workDirProject/GMX/sim1/$input{'pdb'}-$input{'lig'}.imaged.rot.xtc -s $workDirProject/GMX/sim1/topol.tpr -o cluster$i.gro -dump $dump >& trjconv_dump.log";
		`echo 0 |  srun -n 1 trjconv_mpi -f $workDirProject/GMX/sim1/$input{'pdb'}-$input{'lig'}.imaged.rot.xtc -s $workDirProject/GMX/sim1/topol.tpr -o cluster$i.gro -dump $dump >& trjconv_dump.log`;

		print "Trying to energetically minimize cluster number $i...\n";

		# Building system box
		print " echo 0 |  srun -n 1 editconf_mpi -bt $input{'boxtype'} -f cluster$i.gro -o cluster$i.box.gro -d $input{'boxsize'} -c -princ >& cluster${i}_editconf.log >& editconf.log\n";
		`echo 0 |  srun -n 1 editconf_mpi -bt $input{'boxtype'} -f cluster$i.gro -o cluster$i.box.gro -d $input{'boxsize'} -c -princ >& cluster${i}_editconf.log >& editconf.log`;

		print " srun -n 1 grompp_mpi -f $homeDir/GMX_files/$gmx_em_file -c cluster$i.box.gro -p $workDirProject/GMX/$input{'pdb'}-$input{'lig'}_GMX.mod.top -o em.tpr >& grompp.log\n";
		` srun -n 1 grompp_mpi -f $homeDir/GMX_files/$gmx_em_file -c cluster$i.box.gro -p $workDirProject/GMX/$input{'pdb'}-$input{'lig'}_GMX.mod.top -o em.tpr >& grompp.log`;

		print "mpirun -np $input{'nprocs'} mdrun_mpi -s em.tpr -c cluster$i.min.pdb >& mdrun.log\n";
		`mpirun -np $input{'nprocs'} mdrun_mpi -s em.tpr -c cluster$i.min.pdb >& mdrun.log`;

		print " srun -n 1 trjconv_mpi -s $workDirProject/GMX/sim1/topol.tpr -f cluster$i.min.pdb -o cluster$i.min.imaged.pdb -pbc mol -center -ur compact >& trjconv.rot.log\n";
		#`echo "1 1" | trjconv_mpi -s $workDirProject/GMX/sim1/topol.tpr -f cluster$i.min.pdb -o cluster$i.min.imaged.pdb -pbc mol -center -ur compact`;
		`echo "1 1" |  srun -n 1 trjconv_mpi -s $workDirProject/GMX/sim1/topol.tpr -f cluster$i.box.gro -o cluster$i.min.imaged.pdb -pbc mol -center -ur compact >& trjconv.rot.log`;

		$i++;
	}
}

exit if ($input{'onlyCluster'});

#########################  Gaussian Energy Optimization  #########################

# Gaussian Energy optimization 
mkdir("$workDirProject/Gaussian") if (! -e "$workDirProject/Gaussian");
chdir("$workDirProject/Gaussian");

print "\nComputing Gaussian Optimization...  ";

if ($input{'onlyGaussian'}){
	$charge = getCharge("../lig$input{'lig'}.H.pdb");
}
print "Ligand Charge: $charge\n" if ($input{'verbose'});

my $cl_ok = $input{'clModels'};
if (-s "$workDirProject/GMX/cluster.log") {                
	$cl_ok = &checkClusterPercentage("$workDirProject/GMX/cluster.log",$input{'clModels'},$input{'clPercentage'});
	print "Working with $cl_ok clusters (from cluster.log)\n";
}
else{
	print "Output file cluster.log not found, working with $cl_ok clusters (input clModels)\n";
}

foreach my $snapCl (`ls --color=never $workDirProject/CLUSTERS`){
	chomp($snapCl);
	print "SNP: $snapCl\n";

	$snapCl =~ /CL(\d+)/;
	my $numCluster = $1;
	
	next if ($numCluster > $cl_ok);

	mkdir("$workDirProject/Gaussian/$snapCl") if (! -e "$workDirProject/Gaussian/$snapCl");
	chdir("$workDirProject/Gaussian/$snapCl");

	next if (-s "gaussianFinal.pdb");

	# Gaussian Scratch folder
	my $gaussTmp = "$workDirProject/Gaussian/$snapCl/gaussTmp";
	mkdir("$gaussTmp") if (! -e "$gaussTmp");

	# Gaussian Scratch File Prefix
	my $dirGaussianTmp = "$workDirProject/Gaussian/$snapCl/gaussTmp/gaussTmp";

	my $gaussian_zmatrix = "cluster$numCluster.HF.com";
	my $gaussian_zmatrix2 = "cluster$numCluster.B3LYP.com";

	my $inputPDB = "$workDirProject/CLUSTERS/$snapCl/cluster$numCluster.min.imaged.pdb";

	# Fixing Chloride problems in PDB format, so that babel recognizes it correctly.
	my $clCont = `grep Cl $workDirProject/lig$input{'lig'}.H.pdb`;
	if($clCont){
		print "Fixing Chloride problems in PDB format...\n";
		print "$clCont\n";
		print "sed -r 's/  CL(\\w*)/ CL\\1 /g' $inputPDB > pdbFixedCl.pdb\n";
		 `sed -r 's/  CL(\\w*)/ CL\\1 /g' $inputPDB > pdbFixedCl.pdb`; 
		$inputPDB = "pdbFixedCl.pdb";
	}

        # Fixing Bromide problems in PDB format, so that babel recognizes it correctly.
        my $brCont = `grep Br $workDirProject/lig$input{'lig'}.H.pdb`;
        if($brCont){
                print "Fixing Bromide problems in PDB format...\n";
                print "$brCont\n";
                print "sed -r 's/  BR(\\w*)/ BR\\1 /g' $inputPDB > pdbFixedBr.pdb\n";
                 `sed -r 's/  BR(\\w*)/ BR\\1 /g' $inputPDB > pdbFixedBr.pdb`;
                $inputPDB = "pdbFixedBr.pdb";
        }

	# Building Z-matrix from cluster structure energetically minimized
	print "babel -i pdb $inputPDB -o fh $gaussian_zmatrix\n";
	`babel -i pdb $inputPDB -o fh $gaussian_zmatrix\n`;

	print "\nPreparing Gaussian configuration files for step 1: HF level of theory... ";
	my $theory = "HF";
	my $gPrefix = "gaussian.1.$theory";
	&run_gaussian(\%input,$charge,$gaussian_zmatrix,$theory,$gPrefix,$dirGaussianTmp);

	# Checking output. If desired energy has not been reached, run a HF minimization again.
	my $gHF = `grep "Normal termination" $workDirProject/Gaussian/$snapCl/$gPrefix.out`;

	if (!$gHF){

		# Removing core file
		`rm core.*`;

		my $gStep = 2;
		$gaussian_zmatrix = "cluster$numCluster.HF.$gStep.com";

		# Building Z-matrix from output of previous HF run
		print "babel -i log $workDirProject/Gaussian/$snapCl/$gPrefix.out -o com $gaussian_zmatrix\n";
		`babel -i log $workDirProject/Gaussian/$snapCl/$gPrefix.out -o com $gaussian_zmatrix\n`;

		$gPrefix = "gaussian.$gStep.$theory";

		print "\nPreparing Gaussian configuration files for step 1 (repeat $gStep): HF level of theory... ";
		my $theory = "HF";
		&run_gaussian(\%input,$charge,$gaussian_zmatrix,$theory,$gPrefix,$dirGaussianTmp);
	}

	# Building Z-matrix from HF gaussian output
	print "babel -i log $workDirProject/Gaussian/$snapCl/$gPrefix.out -o com $gaussian_zmatrix2\n";
	`babel -i log $workDirProject/Gaussian/$snapCl/$gPrefix.out -o com $gaussian_zmatrix2\n`;

	print "\nPreparing Gaussian configuration files for step 2: B3LYP level of theory... ";
	$theory = "B3LYP";
	$gPrefix = "gaussian.1.$theory";
	&run_gaussian(\%input,$charge,$gaussian_zmatrix2,$theory,$gPrefix,$dirGaussianTmp);

	# Checking output. If desired energy has not been reached, run a HF minimization again.
	$gHF = `grep "Normal termination" $workDirProject/Gaussian/$snapCl/$gPrefix.out`;

	if (!$gHF){

		# Removing core file
		`rm core.*`;
		print "Sorry, something went wrong with B3LYP gaussian run...\n";
	}
	else{
		# Extracting resulting pdb file from Gaussian Energy Optimization
		print "babel -i log $workDirProject/Gaussian/$snapCl/$gPrefix.out -o pdb gaussianFinal.pdb\n";
		`babel -i log $workDirProject/Gaussian/$snapCl/$gPrefix.out -o pdb gaussianFinal.pdb\n`;
	}

	print "\nDone!\n";

	chdir("..");
}

print "Done!\n";

#########################  Subroutines  #########################  

# Command line options processing
sub init(){

my ($input,%limits) = @_;

GetOptions(
	"name=s"=>\$input{'name'},
	"smiles=s"=>\$input{'smiles'},
	"smilesFile=s"=>\$input{'smilesFile'},
	"pdb=s"=>\$input{'pdb'},
	"pdbFile=s"=>\$input{'pdbFile'},
	"lig=s"=>\$input{'lig'},
	"ligFile=s"=>\$input{'ligFile'},
	"reps:i"=> \$input{'numReps'},
	"procs:i"=> \$input{'numProcsPerRep'},
	"len:i"=> \$input{'length'},
	"tlow:f"=> \$input{'ini_temp'},
	"thigh:f"=> \$input{'end_temp'},
	"boxsize:f"=> \$input{'boxsize'},
	"boxtype:s"=>\$input{'boxtype'},
	"replex:i"=> \$input{'replex'},
	"ensemble:s"=>\$input{'ensemble'},
	"cluster:f"=> \$input{'cluster'},
	"clMethod:s"=> \$input{'clMethod'},
        "clModels:s"=> \$input{'clModels'},
        "clPercentage:s"=> \$input{'clPercentage'},
        "clIncrement:s"=> \$input{'clIncrement'},
	"onlyLig"=>\$input{'onlyLig'},
	"onlyParam"=>\$input{'onlyParam'},
	"onlyEquil"=>\$input{'onlyEquil'},
	"onlyREMD"=>\$input{'onlyREMD'},
	"onlyCluster"=>\$input{'onlyCluster'},
	"onlyClusterMin"=>\$input{'onlyClusterMin'},
	"onlyGaussian"=>\$input{'onlyGaussian'},
	"force"=>\$input{'force'},
	"v!"=>\$input{'verbose'},
	"h|help"=> sub { usage() }
) or die (usage());

if ((!$input{'pdb'} || !$input{'lig'}) && (!$input{'smiles'}) && (!$input{'smilesFile'}) && (!$input{'pdbFile'}) && (!$input{'ligFile'})) {
	print "Missing Required Arguments: -pdb + -lig or -pdbFile or -ligFile or -smiles or -smilesFile\n";
	exit;
}

if ($input{'pdbFile'} and ! -s "$input{'pdbFile'}"){
	print "Ups, couldn't find pdbFile $input{'pdbFile'}. Please check it and try again!\n";
	exit;
}

if ($input{'pdbFile'} and !$input{'lig'}){
	print "Missing required argument -lig. Please check it and try again!\n";
	exit;
}

if ($input{'ligFile'} and ! -s "$input{'ligFile'}"){
	print "Ups, couldn't find ligFile $input{'ligFile'}. Please check it and try again!\n";
	exit;
}

$input{'numReps'} = 8 if (!$input{'numReps'});			# 8 replicas
$input{'numProcsPerRep'} = 1 if (!$input{'numProcsPerRep'});	# 1 single proc per replica
$input{'ensemble'} = 'NPT' if(!$input{'ensemble'});		# in NPT ensemble
$input{'ini_temp'} = 298 if(!$input{'ini_temp'});		# from 298K
$input{'end_temp'} = 398 if(!$input{'end_temp'});		# to 398K
$input{'length'} = 10 if (!$input{'length'});			# of 10ns-length
$input{'boxsize'} = 0.8 if (!$input{'boxsize'});		# with 8 Angstroms water box
$input{'boxtype'} = 'octahedron' if (!$input{'boxtype'});	# in a octahedral box
$input{'replex'} = 500 if (!$input{'replex'});			# and 500 MD integration steps between exchange attempts (1ps --> 500 * 2fs time step)
$input{'nprocs'} = $input{'numReps'} * $input{'numProcsPerRep'};
$input{'cluster'} = 0.05 if (!$input{'cluster'});		# Clustering with 0.5 Angstroms RMSd 
$input{'clMethod'} = 'gromos' if (!$input{'clMethod'});		# Clustering Method: gromos
$input{'clModels'} = 10 if (!$input{'clModels'});               # Clustering: Max number of cluster representatives
$input{'clPercentage'} = 95 if (!$input{'clPercentage'});       # Clustering: percentage of trajectory snapshots represented by clusters
$input{'clIncrement'} = 0.005 if (!$input{'clIncrement'});      # Clustering: Cutoff increment for iterative clustering

print "WARNING: Unprocessed input parameter\n" if $ARGV[0];

if (!isdigit($input{'numReps'})) {
	print "Error, number of replicas must be a number ($input{'numReps'})\n";
	die(usage());
}
if (!isdigit($input{'length'})) {
	print "Error, MD simulation time must be a number ($input{'length'})\n";
	die(usage());
}
if (!isdigit($input{'ini_temp'})) {
	print "Error, initial temperature must be a number ($input{'ini_temp'})\n";
	die(usage());
}
if (!isdigit($input{'end_temp'})) {
	print "Error, final temperature must be a number ($input{'end_temp'})\n";
	die(usage());
}
if ( $input{'end_temp'} <= $input{'ini_temp'}) {
	print "Error, final temperature must be greather than initial temperature ($input{'ini_temp'} -- $input{'end_temp'})\n";
	die(usage());
}
if ($input{'length'} > $limits{'maxLength'}) {
	print "Error, maximum allowed REMD simulation time $limits{'maxLength'} ns exceeded ($input{'length'})\n";
	die(usage());
}
if ($input{'numReps'} > $limits{'maxReps'}) {
	print "Error, maximum allowed number of replicas limits{'maxReps'} exceeded ($input{'numReps'})\n";
	die(usage());
}
if ( ($input{'ini_temp'} > $limits{'maxTemp'}) or ($input{'end_temp'} > $limits{'maxTemp'}) ){
	print "Error, maximum allowed REMD temperature $limits{'maxTemp'} K exceeded ($input{'ini_temp'} -- $input{'end_temp'})\n";
	die(usage());
}
if ( ($input{'ini_temp'} < $limits{'minTemp'}) or ($input{'end_temp'} < $limits{'minTemp'}) ){
	print "Error, minimum allowed REMD temperature $limits{'minTemp'} K exceeded ($input{'ini_temp'} -- $input{'end_temp'})\n";
	die(usage());
}
if ($input{'clModels'} > $limits{'maxClModels'}) {
        print "Error, maximum allowed number of cluster conformations $limits{'maxClModels'} exceeded ($input{'clModels'})\n";
        die(usage());
}
if ($input{'clPercentage'} > 100) {
        print "Error, cluster percentage must be 100 or lower ($input{'clPercentage'})\n";
        die(usage());
}
if ($input{'clIncrement'} > $limits{'maxClInc'}) {
        print "Error, maximum allowed increment in iterative clustering $limits{'maxClInc'} exceeded ($input{'clIncrement'})\n";
        die(usage());
}

foreach (@ARGV) {
  print "$_\n";
}
	#my $opt_string = 'hvdf:';
	#getopts( "$opt_string", \%opt ) or usage();
	#usage() if $opt{h};
}

# Help and usage method
sub usage(){

	print STDERR << "EOF";

  $0 Perl script running a pipeline to obtain a ligand conformational ensemble using Gromacs REMD/ Plumed HREX 

  usage: $0 -pdb -lig -smiles [-smilesFile] [-pdbFile] [-ligFile] [-name]
			[-reps] [-procs] [-len] [-replex] [-tlow] [-thigh] [-boxsize] [-boxtype] [-ensemble] [-cluster] [-clMethod]
			[-clModels] [-clPercentage] [-clIncrement] [-only[Lig|Param|Equil|REMD|Cluster|ClusterMin|Gaussian]] -force -v -h

	-name    	  	: project name
	-smiles    	  	: input SMILES code (e.g. CCCCOOO)
	-smilesFile    	  	: input SMILES File 
	-pdb    	  	: input pdbCode (e.g. 1d4l)
	-pdbFile    	  	: input PDB File (e.g. 1d4l.pdb)
	-lig		      	: input ligand 3-letter code (e.g. ATP)
	-ligFile    	  	: input PDB LIG File (e.g. FD2.pdb)

	-reps (8)		: number of MD replicas in Gromacs REMD
	-procs (1)		: number of processors per MD replica
	-len (10)		: length of the MD simulations in REMD (ns)
	-replex	(500)		: number of steps between REMD exchanges
	-tlow (298)		: initial (low) temperature in REMD (K)
	-thigh (398)		: final (high) temperature in REMD (K)
	-boxsize (0.8)		: box size in REMD simulations (nm)
	-boxtype (octahedron)	: box type in REMD simulations (triclinic, cubic, dodecahedron, octahedron)
	-ensemble (NPT)      	: REMD ensemble (NPT,NVT)

        -cluster (0.05)         : Clustering RMSd (nm)
        -clMethod (gromos)      : Clustering Method (gromos, linkage, jarvis-patrick, monte-carlo, diagonalization)
        -clModels (10)          : Clustering: Max number of cluster representatives
        -clPercentage (95)      : Clustering: percentage of trajectory snapshots represented by clusters (%)
        -clIncrement (0.005)    : Clustering: Cutoff increment for iterative clustering (nm)
	
	-onlyLig		: only running ligand extraction
	-onlyParam		: only running ligand parameterization
	-onlyEquil		: only running ligand parameterization + REMD Equilibration
	-onlyREMD		: only running ligand parameterization + REMD
	-onlyCluster		: only running ligand clustering (after REMD)
	-onlyClusterMin		: only running ligand clustering minimization (after REMD)
	-onlyGaussian		: only running Gaussian Optimization

	-force			: force re-computing of all steps

	-h        		: this (help) message
	-v        		: verbose output

	Example 1: $0 -pdb 1efy -lig BZC -reps 8 -procs 2 -ensemble NPT
	Example 2: $0 -ligFile ADP.pdb -lig ADP -reps 8 -procs 2 -ensemble NPT
	Example 3: $0 -pdbFile ../../PDBS/pdb1efy.ent -lig BZC -reps 8 -procs 2 -ensemble NPT
	Example 4: $0 -smiles "O=C(N1CCC(C(=O)OC)CC1)C(NS(=O)(=O)c3cc2ccccc2cc3)Cc4cc(C(=[N\@H])N)ccc4" -reps 8 -procs 2
	Example 5: $0 -smilesFile ../SMILES/FD2.smi -reps 8 -procs 2 -ensemble NPT

EOF
    exit;
}

# Extracting ligand from pdb
sub extractLig{

	my ($lig,$mirror,$pdbfile) = @_;

	# Checking pdb existence
	#my $file = "$mirror/$pdbfile";
	#if (! -s "$file"){
	#	print "Ups, couldn't find pdb $pdbfile in $mirror. Please check it and try again!\n";
	#	exit;
	#}

	# Extracting ligand from pdb
	my $glig = `cat $pdbfile | grep "^HETATM" | grep $lig`;

	my $pdbLig = '';
	if ($glig) {

		# DEBUG: Printing ligand from grep command
		print "$glig\n" if ($input{'verbose'});

		# Getting ALL ligands with residue name $input{'lig'} in the pdb
		my %ligs = &extract_ligs($glig);

		# Do we have more than one ligand in the pdb?
		my $nligs = keys(%ligs);
		if ($nligs < 1) {
			print "Ups, couldn't find lig $lig. Please check it and try again!\n";
			exit;
		}elsif ($nligs > 1) {

			# Case we have more than one ligand. We just get the first one. 
			# We can give the user the possibility to choose which one, maybe in version 2.0 :)
			print "WARNING: We've found more than one ligand ($nligs) with code $lig in the pdb file. Getting the first one!\n";

			# Getting minimum residue number for the different ligands found in the pdb
			#my $n = (sort {$a <=> $b} keys %ligs)[0];
			my $n = (sort keys %ligs)[0];

			# DEBUG: Printing minimum residue number for the different ligands found in the pdb
			print "N: $n\n" if ($input{'verbose'});
			$pdbLig = $ligs{$n};
		}
		else{
			# Case we only have one ligand. Good!
			$pdbLig = $glig;
		}
	}
	else{
		print "Ups, couldn't find lig $input{'lig'} in pdb. Please check it and try again!\n";
		exit;
	}

	return $pdbLig;
}

# Extracting ligs from pdb (if there are more than one)
sub extract_ligs {

	my $lig = shift;

	my %ligs;
	foreach my $line (split '\n',$lig){
		# Line should be like this:
		# HETATM  750  O   HOH A 121      16.241   1.058   6.102  1.00 49.83           O

		# Getting residue number and converting it from text to number
		#my $nres = substr($line,22,4);
		#$nres += 0;
		my $nres = substr($line,21,5);
		$nres =~ s/ //g;

		# Storing ligands in a hash structure
		$ligs{$nres}.= "$line\n";
	}

	# DEBUG: printing all ligands with residue number
	# foreach my $f (sort keys %ligs){
	#	print "$f $input{'lig'}s{$f}";
	# }

	return %ligs;
}

# Finding out ligand charge from OpenBable "clues"
sub getCharge {
	# HETATM   14  N14 E20 A2001       3.836  67.789  66.275  1.00  0.00           N1+

	my $lig = shift;

	my $charge = 0;
	open LIG,"$lig";
	while(<LIG>){
		chomp;
		next if (!/^ATOM/ and !/^HETATM/);
		my @arr = split ' ';
		my $elem = $arr[$#arr];
		$charge ++ if($elem =~ /\+/ and $elem !~ /C1+/);
		$charge -- if($elem =~ /\-/);
	}
	close LIG;
	
	return $charge;
}

# Preparing Gromacs files
sub prepare_gmx_files{

	my ($input,$top,$gro) = @_;

	my $eq_gro_file = "$input{'pdb'}-$input{'lig'}_GMX.gro";

	# Copying gro and itp files generated by ACPYPE
	`cp ../$input{'pdb'}-$input{'lig'}.acpype/$input{'pdb'}-$input{'lig'}_GMX.gro .`;
	`cp ../$input{'pdb'}-$input{'lig'}.acpype/$input{'pdb'}-$input{'lig'}_GMX.itp .`;

	# Building appropiate Gromacs topology file
	# Adding Amber ff and water model
	open TOP,">$input{'pdb'}-$input{'lig'}_GMX.mod.top";
	print TOP "; Gromacs topology file automatically generated by ligParm.pl script\n\n";
	print TOP "; Include forcefield parameters\n";
	print TOP "\#include \"amber99sb-ildn.ff/forcefield.itp\"\n\n";
	print TOP "; Include $input{'pdb'}-$input{'lig'}_GMX.itp topology\n";
	print TOP "\#include \"$input{'pdb'}-$input{'lig'}_GMX.itp\"\n\n";
	print TOP "; Include water topology\n";
	print TOP "\#include \"amber99sb-ildn.ff/tip3p.itp\"\n\n";
	print TOP "[ system ]\n";
	print TOP " $input{'pdb'}-$input{'lig'}\n";
	print TOP "[ molecules ]\n";
	print TOP "; Compound        nmols \n";
	print TOP " $input{'pdb'}-$input{'lig'}         1     \n";
	close TOP;

	print "Done!\nConfiguring Gromacs REMD, $input{'numReps'} reps, from $input{'ini_temp'} to $input{'end_temp'}... \n";

	# Building system box
	`echo 0 | srun -n 1 editconf_mpi -bt $input{'boxtype'} -f $eq_gro_file -o $input{'pdb'}-$input{'lig'}_conf.gro -d $input{'boxsize'} -c -princ >& $input{'pdb'}-$input{'lig'}_editconf.log`;

	# Solvating system
	my $watBox = "spc216.gro";	# TIP3P / SPCE Water model
	#`srun -n 1 gmx_mpi solvate -cp $input{'pdb'}-$input{'lig'}_conf.gro -cs $watBox -o $gro -p $top >& $input{'pdb'}-$input{'lig'}_genbox.log`;
	`srun -n 1 genbox_mpi -cp $input{'pdb'}-$input{'lig'}_conf.gro -cs $watBox -o $gro -p $top >& $input{'pdb'}-$input{'lig'}_genbox.log`;

}

sub equilibrate_REMD {

	my ($input,$eqStep,$top,$gro) = @_;

	my $tempFixed = $input{'ini_temp'};

	my $range = $input{'end_temp'} - $input{'ini_temp'};
	my $step = $range / ($input{'numReps'} - 1);
	my $temp = $input{'ini_temp'};
	for (my $nrep = 1; $nrep <= $input{'numReps'}; $nrep++){

		my $t = int($temp);

		# DEBUG: Printing number of replica and associated temperature	
		print "Rep: $nrep, Temperature: $t\n" if ($input{'verbose'});

		# Creating new dirs for replica $nrep
		my $dirEq = "equil$nrep"."_$eqStep";
		mkdir("$dirEq") if (! -e "$dirEq");

		# Building Gromacs mdp equilibration config files (with appropiate temperature)
		`cp $homeDir/GMX_files/$gmx_eq_files[$eqStep] $dirEq`;	# Copying equilibration config file template
		#`sed -i s/REMD_temp/$t/g $dirEq/$gmx_eq_files[$eqStep]`;	# Changing Temperature in template
		`sed -i s/REMD_temp/$tempFixed/g $dirEq/$gmx_eq_files[$eqStep]`;	# Changing Temperature in template

		# Preparing Gromacs tpr files to equilibrate systems with REMD
		my $eq_tpr_file = "topol.tpr";
		my $maxWarn = 1;
		` srun -n 1 grompp_mpi -f $dirEq/$gmx_eq_files[$eqStep] -c $gro -p $top -o $dirEq/$eq_tpr_file -po $dirEq/mdout.mdp -maxwarn $maxWarn >& $dirEq/eq_grompp.log`;
		print " srun -n 1 grompp_mpi -f $dirEq/$gmx_eq_files[$eqStep] -c $gro -p $top -o $dirEq/$eq_tpr_file -po $dirEq/mdout.mdp -maxwarn $maxWarn >& $dirEq/eq_grompp.log\n";

		$temp+=$step;
	}

	print "Done!\nEquilibrating systems in REMD... \n";

	# Run equilibrations with REMD command
	my $eq_out_file = "$input{'pdb'}-$input{'lig'}-eqOut";
	my $remd_range = "";
	foreach my $f (1..$input{'numReps'}){
		$remd_range .= "equil${f}_$eqStep ";
	}

	print "mpirun -np $input{'nprocs'} mdrun_mpi -v -multidir $remd_range -o $eq_out_file\n";
	`mpirun -np $input{'nprocs'} mdrun_mpi -v -multidir $remd_range -o $eq_out_file`;

	# Check output
	my $dirEq1 = "equil1_$eqStep";
	#if (! -s "$dirEq1/traj_comp.xtc"){
	if (! -s "$dirEq1/$eq_out_file.trr"){
		die ("ERROR in equilibration step $eqStep. Please check log files.\n");
	}
}

sub equilibrate_REMD_oneGRO {

	my ($input,$eqStep,$top,$gro) = @_;

	my $tempFixed = $input{'ini_temp'};

	# Creating new dirs for replica $nrep
	my $dirEq = "equil_$eqStep";
	mkdir("$dirEq") if (! -e "$dirEq");

	# Building Gromacs mdp equilibration config files (with appropiate temperature)
	`cp $homeDir/GMX_files/$gmx_eq_files[$eqStep] $dirEq`;	# Copying equilibration config file template
	#`sed -i s/REMD_temp/$t/g $dirEq/$gmx_eq_files[$eqStep]`;	# Changing Temperature in template
	`sed -i s/REMD_temp/$tempFixed/g $dirEq/$gmx_eq_files[$eqStep]`;	# Changing Temperature in template

	# Preparing Gromacs tpr files to equilibrate systems with REMD
	my $eq_tpr_file = "topol.tpr";
	my $maxWarn = 1;
	`srun -n 1 grompp_mpi -f $dirEq/$gmx_eq_files[$eqStep] -c $gro -p $top -o $dirEq/$eq_tpr_file -po $dirEq/mdout.mdp -maxwarn $maxWarn >& $dirEq/eq_grompp.log`;
	print "srun -n 1 grompp_mpi -f $dirEq/$gmx_eq_files[$eqStep] -c $gro -p $top -o $dirEq/$eq_tpr_file -po $dirEq/mdout.mdp -maxwarn $maxWarn >& $dirEq/eq_grompp.log\n";

	# Run equilibrations with REMD command
	my $eq_out_file = "$input{'pdb'}-$input{'lig'}-eqOut";

	#print "srun -n $input{'nprocs'} gmx_mpi mdrun -s $dirEq/$eq_tpr_file -o $dirEq/$eq_out_file -c $dirEq/confout.gro\n";
	#`srun -n $input{'nprocs'} gmx_mpi mdrun -s $dirEq/$eq_tpr_file -o $dirEq/$eq_out_file -c $dirEq/confout.gro`;
	print "srun -n 1 mdrun_mpi -s $dirEq/$eq_tpr_file -o $dirEq/$eq_out_file -c $dirEq/confout.gro\n";
	`srun -n 1 mdrun_mpi -s $dirEq/$eq_tpr_file -o $dirEq/$eq_out_file -c $dirEq/confout.gro`;

	# Check output
	if (! -s "$dirEq/$eq_out_file.trr"){
		die ("ERROR in equilibration step $eqStep. Please check log files.\n");
	}
}

sub production_REMD {

	my ($input,$eqStep,$top,$gro) = @_;
	
	my $sim_gro_file = "confout.gro";

	my $tempFixed = $input{'ini_temp'};

	# Temperature Replica Exchange
	my $range = $input{'end_temp'} - $input{'ini_temp'};
	my $step = $range / ($input{'numReps'} - 1);
	my $temp = $input{'ini_temp'};
	my $tempEnd = $input{'end_temp'};
	my $t = int($temp);
	my $iniTemp = $t;
	my $te = int($tempEnd);
	my $endTemp = $te;

	# Hamiltonian Replica Exchange 
	#set tmin=300
	#set tmax=1000

	# build geometric progression
	#set list = `awk -v n=$nrep -v tmin=$tmin -v tmax=$tmax 'BEGIN{for(i=0;i<n;i++){t=tmin*exp(i*log(tmax/tmin)/(n-1)); printf(t); if(i<n-1)printf(","); }}'`
	#echo $list

	# Preparing Gromacs REMD Production Runs
	for (my $rep = 1; $rep <= $input{'numReps'}; $rep++){

		#my $t = int($temp);
		my $t = $iniTemp*exp(($rep-1)*log($endTemp/$iniTemp)/($input{'numReps'}-1));
		print "HREX geometrical progression: Rep $rep, Tini: $iniTemp, Tend: $endTemp, T: $t\n";
		print "HREX geometrical progression eq: $iniTemp*exp(($rep-1)*log($endTemp/$iniTemp)/($input{'numReps'}-1)) = $t\n";

		# Creating new dirs for replica $rep
		my $dirEq = "equil${rep}_$eqStep";
		my $dirRun = "sim$rep";
		mkdir("$dirRun") if (! -e "$dirRun");

		# Plumed.dat file required (even if it is empty)
		`touch $dirRun/plumed.dat`;

		# Preparing Gromacs TOP Hamiltonian REMD Production Topologies

		my $include = `grep GMX $top | grep "^#"`;
		#include "1if7-SBR_GMX.itp"

		my ($tag,$fileInc) = split ' ',$include;
		$fileInc =~s/\"//g;
		my $catInc = `cat $fileInc`;

		my $catTop = `cat $top`;
		$catTop =~s/$include/$catInc/g;

		open TOP,">$dirRun/top.hrex$rep.top"; 
		print TOP $catTop;
		close TOP;

		#my $lambda = `echo $list | awk 'BEGIN{FS=",";}{print $1/$'$j';}'`
		my $lambda = $iniTemp / $t;
		#my $lambda = 1.0;

		# Process topology: 1.- Choose "HOT" region/atoms adding "_" to them in the [ atoms ] section of the topology.
		# Adding the "_" to all atoms of the peptide is equivalent to REST2 (the entire solute is scaled).

		my $tophrex = "top.hrex$rep.lambda.top";
		my $tophrexHot = "top.hrex$rep.lambda.hot.top";
		
		&prep_top_hrex("$dirRun/top.hrex$rep.top","$dirRun/$tophrex");
 		#&prep_top_hrex("$dirRun/$tophrex","$dirRun/$tophrexHot");

		# Process topology: 2.- Scaling Hamiltonian terms using plumed partial_tempering script.

		print "srun -n 1 plumed partial_tempering $lambda < $dirRun/$tophrex > $dirRun/$tophrexHot \n";
		`srun -n 1 plumed partial_tempering $lambda < $dirRun/$tophrex > $dirRun/$tophrexHot`;

		# prepare tpr file
		# -maxwarn is often needed because box could be charged
		#print "grompp_mpi_d  -maxwarn 1 -o topol$i.tpr -f grompp$i.mdp -p topol$i.top\n";
		#`grompp_mpi_d  -maxwarn 1 -o topol$i.tpr -f grompp$i.mdp -p topol$i.top`;

		# Building Gromacs mdp production config files (with appropiate temperature)
		`cp $homeDir/GMX_files/$gmx_run_file $dirRun`;		# Copying production config file template
		#`sed -i s/REMD_temp/$t/g $dirRun/$gmx_run_file`;	# Changing Temperature in template
		`sed -i s/REMD_temp/$tempFixed/g $dirRun/$gmx_run_file`;	# Changing Temperature in template

		my $mdl = $input{'length'} * 1000 / 0.002;		# ns to steps ( length to ps / 2ps timestep)
		`sed -i s/md_time/$mdl/g $dirRun/$gmx_run_file`;	# Changing MD length in template 

		# Preparing Gromacs tpr files to run MD productions with REMD
		my $sim_tpr_file = "topol.tpr";
		my $maxWarn = 1;
		#print "grompp_mpi -f $dirRun/$gmx_run_file -c $dirEq/$sim_gro_file -p $top -o $dirRun/$sim_tpr_file -po $dirRun/mdout.mdp -maxwarn $maxWarn >& $dirRun/sim_grompp.log\n";
		#`grompp_mpi -f $dirRun/$gmx_run_file -c $dirEq/$sim_gro_file -p $top -o $dirRun/$sim_tpr_file -po $dirRun/mdout.mdp -maxwarn $maxWarn >& $dirRun/sim_grompp.log`;
		print " srun -n 1 grompp_mpi -f $dirRun/$gmx_run_file -c $dirEq/$sim_gro_file -p $dirRun/$tophrexHot -o $dirRun/$sim_tpr_file -po $dirRun/mdout.mdp -maxwarn $maxWarn >& $dirRun/sim_grompp.log\n";
		` srun -n 1 grompp_mpi -f $dirRun/$gmx_run_file -c $dirEq/$sim_gro_file -p $dirRun/$tophrexHot -o $dirRun/$sim_tpr_file -po $dirRun/mdout.mdp -maxwarn $maxWarn >& $dirRun/sim_grompp.log`;

		#$temp+=$step;
	}

	print "\nRunning production MDs in REMD...  ";

	# Run Gromacs REMD production simulations
	my $sim_out_file = "$input{'pdb'}-$input{'lig'}-simOut";
	my $sim_out_xtc_file_raw = "$input{'pdb'}-$input{'lig'}.xtc";
	my $remd_range = "";
	foreach my $f (1..$input{'numReps'}){
		$remd_range .= "sim$f ";
	}

	print "mpirun -np $input{'nprocs'} mdrun_mpi -hrex -plumed plumed.dat -v -multidir $remd_range -o $sim_out_file -replex $input{'replex'} -x $sim_out_xtc_file_raw\n";
	`mpirun -np $input{'nprocs'} mdrun_mpi -hrex -v -plumed plumed.dat -multidir $remd_range -o $sim_out_file -replex $input{'replex'} -x $sim_out_xtc_file_raw`;

	#print "srun mdrun_mpi -hrex -plumed plumed.dat -v -multidir $remd_range -o $sim_out_file -replex $input{'replex'} -x $sim_out_xtc_file_raw\n";
	#`srun mdrun_mpi -hrex -v -plumed plumed.dat -multidir $remd_range -o $sim_out_file -replex $input{'replex'} -x $sim_out_xtc_file_raw`;

	#print "srun -n $input{'nprocs'} gmx_mpi mdrun -hrex -v -multidir $remd_range -o $sim_out_file -replex $input{'replex'}\n";
	#`srun -n $input{'nprocs'} gmx_mpi mdrun -hrex -v -multidir $remd_range -o $sim_out_file -replex $input{'replex'}`;

	#print "mpirun -np $input{'nprocs'} mdrun_mpi -v -multidir sim$remd_range -o $sim_out_file -replex $input{'replex'} -pin on\n";
	#`mpirun -np $input{'nprocs'} mdrun_mpi -v -multidir sim$remd_range -o $sim_out_file -replex $input{'replex'} -pin on`;

	#print "srun -n $input{'nprocs'} gmx_mpi mdrun -hrex -plumed plumed.dat -v -multidir $remd_range -o $sim_out_file -replex $input{'replex'}\n";
	#`srun -n $input{'nprocs'} gmx_mpi mdrun -hrex -v -plumed plumed.dat -multidir $remd_range -o $sim_out_file -replex $input{'replex'}`;

	##print "mpirun -np $input{'nprocs'} gmx_mpi mdrun -hrex -plumed plumed.dat -v -multidir $remd_range -o $sim_out_file -replex $input{'replex'}\n";
	##`mpirun -np $input{'nprocs'} gmx_mpi mdrun -hrex -v -plumed plumed.dat -multidir $remd_range -o $sim_out_file -replex $input{'replex'}`;

	#`mpirun -np $input{'nprocs'} mdrun_mpi -multidir $remd_range -o $sim_out_file -replex $input{'replex'}`;
	#`mpirun -np $input{'nprocs'} mdrun_mpi -v -multidir sim$remd_range -o $sim_out_file `;

	# Post-processing Gromacs REMD trajectories generated

	# Extracting ligand pdb from MD system
	`echo 0 |  srun -n 1 trjconv_mpi -f $gro -s sim1/topol.tpr -o $input{'pdb'}-$input{'lig'}.ret.pdb -ur compact -pbc atom >& trjconv.retPdb.log`;

	# Plotting REMD temperature exchange between replicas
	#mkdir("DEMUX") if(! -s "DEMUX");
	#chdir("DEMUX");
	chdir("sim1");

	# De-multiplexing REMD - HREX trajectory
	`demux.pl md.log`;
	#`demux.pl ../sim1/md.log`;
	#my $simTrajs = '';
	#foreach my $f (1..$input{'numReps'}){
	#	$simTrajs .= "../sim$f/traj.xtc ";
	#}

	#`trjcat_mpi -f $simTrajs -demux replica_index.xvg -o traj.xtc >& trjcat_demux.log`;

	# Imaging trajectory sim1 
	my $sim_tpr_file = "topol.tpr";
	my $sim_out_xtc_file = "$input{'pdb'}-$input{'lig'}.imaged.xtc";
	my $sim_out_xtc_rot_file = "$input{'pdb'}-$input{'lig'}.imaged.rot.xtc";

	# Image trajectory (moving ligand to center of the water box)
	#`echo "1 0" | trjconv_mpi -s $sim_tpr_file -f traj_comp.xtc -o $sim_out_xtc_file -pbc mol -center -ur compact >& trjconv_image.log`;
	#`echo "1 0" | trjconv_mpi -s $sim_tpr_file -f 0_traj.xtc -o $sim_out_xtc_file -pbc mol -center -ur compact >& trjconv_image.log`;
	`echo "1 0" |  srun -n 1 trjconv_mpi -s $sim_tpr_file -f $sim_out_xtc_file_raw -o $sim_out_xtc_file -pbc mol -center -ur compact >& trjconv_image.log`;

	# Removing rotation
	`echo "1 0" |  srun -n 1 trjconv_mpi -s $sim_tpr_file -f $sim_out_xtc_file -o $sim_out_xtc_rot_file -fit rot+trans >& trjconv_rot.log`;

	chdir("..");

	# Imaging all trajectories (Do we really need it???)
	#for (my $rep = 1; $rep <= $input{'numReps'}; $rep++){

	#	my $dirRun = "sim$rep";
	#	chdir($dirRun);

	#	my $sim_tpr_file = "topol.tpr";
	#	my $sim_out_xtc_file = "$input{'pdb'}-$input{'lig'}.imaged.xtc";
	#	my $sim_out_xtc_rot_file = "$input{'pdb'}-$input{'lig'}.imaged.rot.xtc";

		# Image trajectory (moving ligand to center of the water box)
	#	`echo "1 1" | trjconv_mpi -s $sim_tpr_file -f traj.xtc -o $sim_out_xtc_file -pbc mol -center -ur compact >& trjconv_image.log`;

		# Removing rotation
	#	`echo "1 1" | trjconv_mpi -s $sim_tpr_file -f $sim_out_xtc_file -o $sim_out_xtc_rot_file -fit rot+trans >& trjconv_rot.log`;

		# Reducing traj??
		# $TCONV -f pre-dimer-2.10-20.imaged.rot.1.xtc -o pre-dimer-2.10-20.imaged.rot.10.xtc -dt 10 >& pre-dimer-2.10-20.trjconv_img_rot_10.log

	#	chdir("..");
	#}

	my $number_frames = $input{'length'} * 1000; # ns to ps
	my $check_frames = `grep "frame   $number_frames" sim1/trjconv_image.log`;

	if ($check_frames) {
		return 1;
	}
	else{
		return 0;
	}
}

sub run_gaussian {
	my ($input,$charge,$gaussian_zmatrix,$theory,$gauPrefix,$dirGaussianTmp) = @_;

	my $nprocs = $input{'numReps'} * $input{'numProcsPerRep'};

	# Input Z-matrix, we need to remove the first 2 lines of the file (blank line + nº of atoms)
	my $zmat;
	$zmat = `tail -n+6 $gaussian_zmatrix` if ($theory eq "HF" and $gauPrefix ne "gaussian.1.HF");
	$zmat = `tail -n+3 $gaussian_zmatrix` if ($theory eq "HF" and $gauPrefix eq "gaussian.1.HF");

	if ($theory eq "B3LYP"){
		$zmat = `tail -n+6 $gaussian_zmatrix`;	# Removing first 6 lines (just interested in coordinates)
		my @lines = split /\n/, $zmat;		# Removing last 2 lines (blank lines)
		$zmat = join "\n", @lines;
		$zmat .= "\n";				# Adding 1 blank line at the end (previous 2 commands remove ALL blank lines)
	}

	my $gauInput = "$gauPrefix.in";
	my $gauOutput = "$gauPrefix.out";
	my $gauLog = "$gauPrefix.log";
	#my $gauOutput = "gaussian.$theory.out";
	#my $gauLog = "gaussian.$theory.log";

	my $theoryLevel = '';
	if($theory eq "B3LYP"){
		$theoryLevel = "p B3LYP/6-31G* OPT=(Z-matrix,loose,Maxcycles=1000,MaxStep=60) SCRF=(IEFPCM,Read,solvent=water)";
	}
	elsif($theory eq "HF"){
		$theoryLevel = "p HF/3-21G OPT=(Z-matrix,loose,Maxcycles=1000,MaxStep=60) SCRF=(IEFPCM,Read,solvent=water)";
	}
	else{
		print "Sorry, level of theory in gaussian energy optimization is not defined, exiting...\n";
		exit;
	}

	# Building gaussian configuration file.
	open GAU, ">$gauInput";
	print GAU "\%NProcShared=$nprocs\n";
	print GAU "\%rwf=$dirGaussianTmp\n";
	print GAU "\%chk=$dirGaussianTmp\n";
	print GAU "\%int=$dirGaussianTmp\n";
	print GAU "\%d2e=$dirGaussianTmp\n";
	print GAU "\%scr=$dirGaussianTmp\n";
	#print GAU "\%inp=$dirGaussianTmp\n";
	print GAU "\%NoSave\n";
	print GAU "\%Mem=10GB\n";
	print GAU "\#$theoryLevel\n";
	print GAU "\n";
	print GAU " Gaussian ligand cluster HF optimization (generated by ligParm)\n";
	print GAU "\n";
	print GAU "$charge 1 \n";
	print GAU "$zmat\n";
	print GAU "RADII=PAULING\n";
	print GAU "OFAC=0.8\n";
	print GAU "RMIN=0.5\n\n";
	print GAU "MSTModel=DFT cav g03defaults\n";
	print GAU "\n\n";
	close GAU;

	# Building gaussian queue csh.
	#open GAUCSH, ">$gauInput.csh";
	#print GAUCSH "\#\$ -cwd\n";
	#print GAUCSH "\#\$ -pe onemachine $nprocs\n";
	#print GAUCSH "\#\$ -R y\n";
	#print GAUCSH "\#\$ -q cpu.q\n";
	#print GAUCSH "\#\$ -N g$theory-$input{'pdb'}-$input{'lig'}\n\n";
	#print GAUCSH "module load gaussian/09\n\n";
	#print GAUCSH "g09 $gauInput $gauOutput >& $gauLog\n";
	#close GAUCSH;

	# Running Gaussian.
	`g09 $gauInput $gauOutput >& $gauLog`;
}

sub prep_top_hrex {
	my ($in,$out) = @_;

	my $now = 0;

	open OUT,">$out";
	open IN,"$in";
	while(<IN>){	
		#[ atoms ]
		#;   nr  type
		if($now){
			next if (/^;/);
			if (/^\s$/){
				print OUT $_;
				next;
			}
			if(/^\[/){
				print OUT $_;
				$now = 0;
				next;
			}

			#;   nr  type  resi  res  atom  cgnr     charge      mass       ; qtot   bond_type
			#     1   c3     1   SBR   C01    1     0.049300     12.01000

			my $line = $_;
			my ($nr,$type,$resi,$res,$atom,$cgnr,$charge,$mass,$rest) = split ' ';
			chomp($rest);
			print OUT "$nr ${type}_ $resi $res $atom $cgnr $charge $mass $rest ; Plumed hrex HOT atom\n";
		}
		else {
			print OUT "$_";
		}

		if(/\[ atoms \]/){
			#print $_;
			$now = 1;
		}
	}
	close IN;
	close OUT;
}

# Removing hydrogen atoms from a pdb file.
sub removeH {

        my ($fin,$fout) = @_;

        open (OUT,">$fout");

        open(PDB,$fin);
        while (<PDB>){
                next if(!(/^ATOM/ or /^HETATM/ or /^TER/));
                my $at=substr($_,12,4);
                $at=~s/ //g;
                if($at!~/^H/ and $at!~/^\dH/){
                        if($at!~/OXT/){
                                print OUT $_;
                        }
                }
        }
        close(PDB);
        close(OUT);
}

# Check cluster representation percentage throughout the trajectory
sub checkClusterPercentage {

        my ($file,$nclusters,$cl_perc) = @_;

        my $totalSnapshots = $input{'length'} * 1000;
        my $numTotalCl = 0;
        my $numCl = 0;
        my $sum = 0;
        open FILE,"$file";
        while  (<FILE>){
                chomp;

                # Found 2632 clusters
                if (/Found (\d+) clusters/){
                        $numTotalCl = $1;
                        print "NUM CLUSTERS: $numTotalCl\n";
                        next;
                }

                if (/^\s+\d+/){
                        #print "$_\n";

                        #  2 | 114  0.055 |   5479 .045 |    200    219    221    223    229    245    258
                        my @array = split '\|';
                        my $n = $array[1];
                        my @array2 = split ' ',$n;
                        my $numC = $array2[0];
                        $sum+=$numC;
                        my $p = ($sum / $totalSnapshots) * 100;

                        $numCl++;

                        print "Perc: $p, NumCl: $numCl, Sum: $sum, NumC: $numC\n"; #if ($input{'verbose'});

                        if ($p > $cl_perc){
                                if ($numCl > $nclusters){
                                        return 0;
                                }
                                else{
                                        return $numCl;
                                }
                        }
                }
        }
        close FILE;

	return 0;

}

