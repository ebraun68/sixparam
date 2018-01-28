#!/usr/bin/env perl

use warnings;
use strict;
use Data::Dumper;

############################################################################
# Set the global variables
############################################################################
my($progname) = $0;
my($version) = "0.9a";

my($hter);
my($iter);
my($jter);
my($kter);
my($lter);
my($mter);
my($nter);

my($iqexec) = "iqtree-omp-1.5.5"; # iqtree executable
my($threads) = " -nt 2";          # number of threads is multithread iqtree is used
# my($threads) = " ";             # remove comment to use with serial iqtree
my($verbose) = 0; # 0 for quiet, 1 for verbose

my($likelihood);

# Variables to store parameter values
my @paramcode; # v=volume; p=polarity; c=composition; a=aromaticity; t=transversion; g=gencode
my @optimize;  # same as @parameters (1=optimize; 0=do not optimize)
# Holds values from the previous round of optimization
my @oldparam; # 0=volume; 1=polarity; 2=composition; 3=aromaticity; 4=transversion; 5=gencode; 6=lnL

# Optimization increments (values for up to 8 rounds)
my @increment = (0.5,0.1,0.05,0.01,0.0005,0.0001,0.00001,0.00001);

# array to hold the parameter estimates
my @PV;	# 0=volume; 1=polarity; 2=composition; 3=aromaticity; 4=tv; 5=gencode
$PV[6] = "+F"; # amino acid frequency optimization
$PV[7] = ""; # rate heterogeneity

my @RV;	# used to round the values to avoid floating point problems
		# values rounded to precision of X.XXXXX
		
# array to hold identity of parameters to optimize
my @OptPV;	# 0=volume; 1=polarity; 2=composition; 3=aromaticity; 4=tv; 5=gencode
			# Note - 1=optimize; 0=do not optimize
			
			
############################################################################
# Read the run parameters (the control file and the output file)
############################################################################
if ( @ARGV != 2) {
	print "Usage:\n  \$ $progname <ctlfile> <outfile>\n";
	print "  ctlfile = control file\n";
	print "  outfile = tab-delimited output file\n";
	print "exiting...\n";
	exit;
}

my($ctlfile)=$ARGV[0];
my($outfile)=$ARGV[1];

############################################################################
# Read the control file
############################################################################
print "Reading the ctlfile file... ";
open (my $CTLF, $ctlfile) or die "Could not open file $ctlfile for input.\n";
my @ctllist = <$CTLF>; # Read the control file
close($CTLF) or die "Could not close file $ctlfile.\n";
my($modelnum) = $#ctllist;
print "Optimize $modelnum models\n";
$modelnum++; # incremented because the first line of the control file is the filenames and run parameters, not a model

# Now assign the filename variables for the analysis
my($datafname);
my($treefname);
my($basemodel) = "poiss"; # base model (typical default given)
my($genmat) = "universal_code"; # genetic code matrix (".dat" or "_tv.dat" will be appended)

# Variables for the control and estimation of rate heterogeneity and amino acid frequencies
my($aafmethod)  = "ml"; # "empirical" or "ml" - method to estimate amino acid frequencies
my($ratemethod) = "gamma";	# "all", "equal", "inv", "gamma", or "invgamma"
							# "all" tests among four models (E,I,G, and I+G), the other flags
							# limit analysis to a specific rate heterogeneity methods 
my($aafreqs) = "+F"; # holds the optimization output, starting value is just empirical frequencies
my($ratehet) = "";

($datafname,$treefname,$basemodel,$genmat,$aafmethod,$ratemethod) = split(/\s+/, $ctllist[0]);

############################################################################
# Read the dataset (relaxed phylip format) to obtain the number of taxa
# and sites (to calculate AICc). Calculate "base" number of free parameters
############################################################################
open (my $PHYF, $datafname) or die "Could not open file $datafname for input.\n";
my @phylist = <$PHYF>; # Read the phylip infile file
close($PHYF) or die "Could not close file $datafname.\n";

my($ntaxa);
my($nchar);
($ntaxa,$nchar) = split(/\s+/, $phylist[0]);
chomp($nchar);

my($nbranches) = ( 2 * $ntaxa ) - 3;
my($nbaseparam);
if ( $aafmethod ne "empirical" ) { $nbaseparam = $nbranches + 19; }
if ( $ratemethod eq "gamma" || $ratemethod eq "inv" ) { $nbaseparam++; }
if ( $ratemethod eq "invgamma" ) { $nbaseparam = $nbranches + 2; }

print "Analyze dataset $datafname using tree $treefname\n";
print "ntaxa = $ntaxa -- nsites = $nchar\n";
print "Free parameters:\n";
print "  branch lengths = $nbranches\n";
print "  total          = $nbaseparam\n";
print "(plus the rate matrix parameters)\n\n";
print " *** STARTING ANALYSES ***\n\n\n";

my($aic);	my($aicc);	my($nfreeparam);

############################################################################
# First, estimate the equilibrium amino acid frequencies and among-sites
# rate heterogeneity parameters (alpha for +G; proportion inv for +I)
# using a starting matrix (a "good enough" matrix like LG).
############################################################################
my($startingmat) = "PAMmatrices/LG.dat"; # starting matrix
my($alpha);	my($pinv);	my @aafreq;

print "Optimize amino acid frequency and rate heterogeneity parameters using $startingmat\n";
system("cp $startingmat TEMPMAT"); 	# copy the starting matrix to "TEMPMAT"
									# TEMPMAT is used by "est_rates_aaf()"
($ratehet,$aafreqs) = est_rates_aaf($aafmethod,$ratemethod);
print "  -- $ratehet\n";
print "  -- $aafreqs\n";

if ( $verbose == 1 ) { print "\nrate heterogeneity extension = $ratehet\n"; }

############################################################################
# Optimize the models
############################################################################
my($reoptprecision) = 0.0005;	
my($oldlike) = 0;	my($likediff) = 1;

open (my $OUTF, ">$outfile") or die "Could not open file $outfile for output.\n";

# Print the header for the outfile
print $OUTF "#\tModel\tV\tP\tC\tA\tT\tG\t -lnL\t";
print $OUTF "V\tP\tC\tA\tT\tG\tRate Het\tAA freqs\t";
print $OUTF "Dataset\tTree\tBase Model\tTaxa\tSites\t";
print $OUTF "K\tAIC\tAICc\n";

for ($lter=1; $lter<$modelnum; $lter++) {

	# Set the model before doing the optimization
	(@OptPV) = split(/\s+/, $ctllist[$lter]);
	#	$OptPV[0] = volume
	#	$OptPV[1] = polarity
	#	$OptPV[2] = composition
	#	$OptPV[3] = aromaticity
	#	$OptPV[4] = transversion parameter
	#	$OptPV[5] = gencode
	#	$OptPV[6] = model name
	chomp($OptPV[6]);
	print "\n\n*** OPTIMIZING MODEL #";
	print "$lter (";
	print "$OptPV[6]";
	print ") ***\n\n";
	
	# Zero out the starting value for any parameter that is not optimized
	for ($iter=0; $iter<6; $iter++) {
		if ( $OptPV[$iter] == 0 ) { $PV[$iter] = 0; }
	}
	
	$oldlike  = 0;
	$likediff = 1;
	$hter = 0;

	while ( $likediff > $reoptprecision ) {

		# note that $hter increases each round because we assume each round will yield
		# parameter estimates closer to the final estimates

		# the $iter and $jter loops optimize the parameters of interest (V, P, C, A, T, and G)
		for ($iter=$hter; $iter<8; $iter++) { # outer loop is for the increment precision
			for ($jter=0; $jter<6; $jter++) { # inner loop is for the parameter
				if ( $OptPV[$jter] > 0 ) {
					$PV[$jter] = optimize($jter,$increment[$iter],$PV[0],$PV[1],$PV[2],$PV[3],$PV[4],$PV[5],$aafreqs,$ratehet);
				}
			}
		}

		# Rounds parameter estimates to a precision of X.XXXXX
		for ($jter=0; $jter<6; $jter++) {
			if ( $PV[$jter] < 0 ) {
				$RV[$jter] = int( ($PV[$jter] * -100000.0) + 0.5 );
				$PV[$jter] = $RV[$jter] / -100000.0;
			}
			else {
				$RV[$jter] = int( ($PV[$jter] * 100000.0) + 0.5 );
				$PV[$jter] = $RV[$jter] / 100000.0;
			}
		}
		$PV[6] = "$aafreqs";
		$PV[7] = "$ratehet";

		print "\n";
		print "volume      = $PV[0]\n";
		print "polarity    = $PV[1]\n";
		print "composition = $PV[2]\n";
		print "aromaticity = $PV[3]\n";
		print "tv_penalty  = $PV[4]\n";
		print "gencode     = $PV[5]\n";
		print "rates       = $ratehet\n";
		print "aa_freqs    = $aafreqs\n";

		$likelihood = calc_lnL(@PV); # lnL may change slightly with rounding of parameters, so recalculate it
		print "lnL = $likelihood\n\n";
	
		# Check the difference between the old lnL and the new likelihood
		$likediff = abs ( $likelihood - $oldlike );
		$oldlike = $likelihood;
	
		if ( $likediff > $reoptprecision ) {
			# This reoptimizes the amino acid frequencies and rate heterogeneity parameters
			# using the matrix generated by generate-biophysical-matrix (the program called
			# by the calc_lnL subroutine
			print "Reoptimizing amino acid frequency and rate heterogeneity parameters\n";
			system("cp temp-biophys.dat TEMPMAT"); 	# copy the temp-biophys.dat to "TEMPMAT"
			($ratehet,$aafreqs) = est_rates_aaf($aafmethod,$ratemethod);
		}
	
		$hter++;
		if ( $hter == 8 ) {
			print "Optimization exceeded maximum number of rounds.\n";
			print "Terminated with final likelihood difference = $likediff\n";
			$likediff = 0; # zero out $likediff to terminate the while loop
		}

	} # while ( $likediff > $reoptprecision...

	# Output the data to the outfile in tab-delimited format:
	# first, print the model name and the appropriate flags for parameters
	# to optimize (1=optimize; 0=do not optimize)
	print $OUTF "$lter\t$OptPV[6]\t"; # model # and model name
	for ($iter=0; $iter<6; $iter++) { print $OUTF "$OptPV[$iter]\t"; }
	# second, print the log-likelihood 
	print $OUTF "$likelihood\t";
	# third, print the estimated model parameters
	for ($iter=0; $iter<6; $iter++) { print $OUTF "$PV[$iter]\t"; }
	# fourth, print the rate heterogeneity or the amino acid frequencies
	print $OUTF "$ratehet\t$aafreqs\t";
	# fifth, print the base model, data file, and tree file
	print $OUTF "$datafname\t$treefname\t$basemodel\t";
	
	# finally, calculate and output the AIC and AICc
	$nfreeparam = $nbaseparam + $OptPV[0] + $OptPV[1] + $OptPV[2] + $OptPV[3] + $OptPV[4] + $OptPV[5];
	$aic = (2.0 * $likelihood) + (2.0 * $nfreeparam);
	$aicc = $aic + ((2.0 * $nfreeparam * ($nfreeparam+1.0)) / ( $nchar - $nfreeparam - 1.0));
	print $OUTF "$ntaxa\t$nchar\t$nfreeparam\t$aic\t$aicc\n";

} # end for ($lter=1; $lter<$modelnum; $lter++ ...

close($OUTF) or die "Could not close file $outfile\n";


exit;

############################################################################
# subroutine to estimate rate heterogeneity (E,I,G,I+G) and amino acid  
# frequencies. This matrix assumes that a PAML format amino acid rate 
# matrix file is available and named "TEMPMAT" 
############################################################################
sub est_rates_aaf {

	# Set estimation methods using input variables - this allows user choice between the
	# use of empirical amino acid frequencies or ML estimates and (for rate heterogeneity)
	# between choosing the best model (among four possibilities) or limiting choice to a
	# single rate heterogeneity model.
	my($AAmethod,$RTmethod) = @_;
	print "  optimization methods: $AAmethod -- $RTmethod\n";
	my($AAfreqmethod) = "FO";
	if ( $AAmethod eq "empirical" ) { $AAfreqmethod = "F"; }
	my($RATEmethod) = "E,I,G,I+G";
	if ( $RTmethod eq "equal" ) { $RATEmethod = "E"; }
	elsif ( $RTmethod eq "inv" ) { $RATEmethod = "I"; }
	elsif ( $RTmethod eq "gamma" ) { $RATEmethod = "G"; }
	elsif ( $RTmethod eq "invgamma" ) { $RATEmethod = "I+G"; }	
	
	# variables for analysis
	my($lgmodel);
	my($lgalpha) = 0;
	my($lgpinv) = 0;
	my @lgaaf;
	
	# output variable
	my @localrateaaf;
	$localrateaaf[0] = "";
	$localrateaaf[1] = "+F{,";
	# Note - IQ-TREE 1.5.5 throws an error when a value like:
	#        +F{0.0977,0.0477,0.0393,0.0628,0.0149,0.0342,0.0728,0.0670,0.0156,0.0804,0.1000,0.0763,0.0256,0.0230,0.0269,0.0555,0.0555,0.0031,0.0153,0.0863}
	# is passed. So I have added a comma after the first { and one before the last }. This 
	# allowed me to fix the problem by appending the amino acid frequencies to the rate matrix
	
	# set up the command to estimate
	my($cmdr) = "$iqexec -s $datafname -pre est_rates_aaf_temp -st AA -m TESTONLY -mset TEMPMAT -mfreq $AAfreqmethod -mrate $RATEmethod";
	if ( $verbose == 1 ) { print "$cmdr\n"; }
	system("$cmdr -merit AICc -te $treefname -redo $threads > est_rates_aaf_temp.screen");
	
	system("grep \"^Model of substitution:\" est_rates_aaf_temp.iqtree > est_rates_aaf_temp.model");
	system("grep \"^  pi(\" est_rates_aaf_temp.iqtree >> est_rates_aaf_temp.model");
	system("grep \"^Proportion of invariable sites:\" est_rates_aaf_temp.iqtree >> est_rates_aaf_temp.model");
	system("grep \"^Gamma shape alpha:\" est_rates_aaf_temp.iqtree >> est_rates_aaf_temp.model");
	
	open (my $LGF, "est_rates_aaf_temp.model") or die "Could not open file oest_rates_aaf_temp.model for input\n";
	my @lgfline = <$LGF>; # Read the temporary file with the model information
	close($LGF) or die "Could not close file est_rates_aaf_temp.model\n";

	my @lgvals = split(/\s+/, $lgfline[0]);
	$lgmodel = $lgvals[3];
	for ($nter=0; $nter<20; $nter++) { # iterate through the amino acid frequencies
		@lgvals = split(/\s+/, $lgfline[$nter+1]);
		$lgaaf[$nter] = $lgvals[3];
		$localrateaaf[1] = "$localrateaaf[1]" . "$lgaaf[$nter]";
		if ( $nter < 19 ) { $localrateaaf[1] = "$localrateaaf[1]" . ","; }
		else { $localrateaaf[1] = "$localrateaaf[1]" . ",}"; } # see above for comment about comma before }
	}
	if ( $lgmodel eq "TEMPMAT+FO+I+G4" ) {
		@lgvals = split(/\s+/, $lgfline[21]);
		$lgpinv = $lgvals[4];
		@lgvals = split(/\s+/, $lgfline[22]);
		$lgalpha = $lgvals[3];
		$localrateaaf[0] = "+I{" . "$lgpinv" . "}+G4{" . "$lgalpha" . "}";
	}
	elsif ( $lgmodel eq "TEMPMAT+FO+G4" ) {
		@lgvals = split(/\s+/, $lgfline[21]);
		$lgalpha = $lgvals[3];
		$localrateaaf[0] = "+G4{" . "$lgalpha" . "}";
	}
	elsif ( $lgmodel eq "TEMPMAT+FO+I" ) {
		@lgvals = split(/\s+/, $lgfline[21]);
		$lgpinv = $lgvals[4];
		$localrateaaf[0] = "+I{" . "$lgpinv" . "}";
	}
	elsif ( $lgmodel eq "TEMPMAT+F+I+G4" ) { # These elsif statements cover the use of empirical frequencies
		@lgvals = split(/\s+/, $lgfline[21]);
		$lgpinv = $lgvals[4];
		@lgvals = split(/\s+/, $lgfline[22]);
		$lgalpha = $lgvals[3];
		$localrateaaf[0] = "+I{" . "$lgpinv" . "}+G4{" . "$lgalpha" . "}";
	}
	elsif ( $lgmodel eq "TEMPMAT+F+G4" ) {
		@lgvals = split(/\s+/, $lgfline[21]);
		$lgalpha = $lgvals[3];
		$localrateaaf[0] = "+G4{" . "$lgalpha" . "}";
	}
	elsif ( $lgmodel eq "TEMPMAT+F+I" ) {
		@lgvals = split(/\s+/, $lgfline[21]);
		$lgpinv = $lgvals[4];
		$localrateaaf[0] = "+I{" . "$lgpinv" . "}";
	}
	
	print "model = $lgmodel\n";
	print "  alpha = $lgalpha\n";
	print "  pinv  = $lgpinv\n";
	print "@lgaaf\n";
	
	return @localrateaaf;
}

############################################################################
# optimize parameter (V, P, C, A, T, or G) subroutine
############################################################################
sub optimize {

	my @LP;	my @NP;	my($oldparam); 
	my($newLNL);	my($bestLNL);
	my($parameter,$incr,$Vv,$Pv,$Cv,$Av,$Tv,$Gv,$AAFv,$RATEv) = @_;
	$LP[0] = $Vv;		$NP[0] = $Vv;
	$LP[1] = $Pv;		$NP[1] = $Pv;
	$LP[2] = $Cv;		$NP[2] = $Cv;
	$LP[3] = $Av;		$NP[3] = $Av;
	$LP[4] = $Tv;		$NP[4] = $Tv;
	$LP[5] = $Gv;		$NP[5] = $Gv;	
	$LP[6] = $AAFv;		$NP[6] = $AAFv;
	$LP[7] = $RATEv;	$NP[7] = $RATEv;
	# $parameter is the parameter to be optimized (all other parameters are fixed)
	# $incr is the increment used to optimize the focal parameter, $AAFv and $RATEv
	# are the amino acid frequency and rate heterogeneity settings
	
	# check whether it is better to increase or decrease the parameter
	my($direction) = $incr; # + for increase; - for decrease
	
	my($origLNL) = calc_lnL(@LP);
	print "orig\t$origLNL\t$parameter\t$LP[0]\t$LP[1]\t$LP[2]\t$LP[3]\t$LP[4]\t$LP[5]\n";
	$bestLNL = $origLNL;

	$NP[$parameter] = $LP[$parameter] + $incr;
	my($upLNL) = calc_lnL(@NP);
	print "up\t$upLNL\t$parameter\t$NP[0]\t$NP[1]\t$NP[2]\t$NP[3]\t$NP[4]\t$NP[5]\n";
	
	$NP[$parameter] = $LP[$parameter] - $incr;
	my($downLNL) = calc_lnL(@NP);
	print "down\t$downLNL\t$parameter\t$NP[0]\t$NP[1]\t$NP[2]\t$NP[3]\t$NP[4]\t$NP[5]\n";
	
	if ( $origLNL == $downLNL && $origLNL == $upLNL ) { return $LP[$parameter]; }
	if ( $origLNL < $downLNL && $origLNL < $upLNL ) { return $LP[$parameter]; }
	if ( $downLNL < $origLNL ) { 
		$newLNL = $downLNL;
		$direction = -1.0 * $incr; 
	}
	else { 
		$NP[$parameter] = $LP[$parameter] + $incr;
		$newLNL = $upLNL;
	}
	
	# Hill climbing algorithm:
	# increase (or decrease) the parameter value until the likelihood begins to decrease
	while ( $newLNL < $bestLNL ) {
		$bestLNL=$newLNL;
		$LP[0] = $NP[0];	$LP[1] = $NP[1];
		$LP[2] = $NP[2];	$LP[3] = $NP[3];
		$LP[4] = $NP[4];	$LP[5] = $NP[5];
		$NP[$parameter] = $LP[$parameter] + $direction;
		
		$newLNL = calc_lnL(@NP);
		print "$bestLNL\t$newLNL\t$parameter\t$NP[0]\t$NP[1]\t$NP[2]\t$NP[3]\t$NP[4]\t$NP[5]\n";
		if ( $newLNL > $bestLNL ) { 
			$NP[$parameter] = $LP[$parameter];
		}
	}
	
	return $NP[$parameter];
	
}

############################################################################
# calc_lnL subroutine (uses iqtree for likelihood calculations)
############################################################################
sub calc_lnL {
	my($vol,$pol,$comp,$arom,$tv,$gc,$aaf,$rth) = @_;
	
	# Generate the rate matrix
	if ( $verbose == 1 ) { print "./generate-biophysical-matrix $basemodel $genmat temp-biophys.dat $vol $pol $comp $arom $tv $gc 0"; }
	# Note: generate-biophysical-matrix is called with 0 as the last term to suppress printing amino acid frequecies
	system("./generate-biophysical-matrix $basemodel $genmat temp-biophys.dat $vol $pol $comp $arom $tv $gc 0");
	
	# append the amino acid frequencies
	# This code takes a string like +F{,0.0977,0.0477,0.0393,0.0628,0.0149,0.0342,0.0728,0.0670,0.0156,0.0804,0.1000,0.0763,0.0256,0.0230,0.0269,0.0555,0.0555,0.0031,0.0153,0.0863,}
	# and appends the amino acid frequencies included in the string (between commas) to a rate matrix files
	my @localaaf = split(/,/, $aaf);
	open (my $AAFF, ">>temp-biophys.dat") or die "Could not open file temp-biophys.dat for append\n";
	print $AAFF "\n";
	for ($mter=1; $mter<21; $mter++) { print $AAFF "$localaaf[$mter] "; }
	print $AAFF "\n";
	close($AAFF) or die "Could not close file temp-biophys.dat\n";

	# Run IQ-TREE with the test matrix
	my($cmd) = "$iqexec -s $datafname -pre opt_biophys_temp -st AA -m temp-biophys.dat" . "$rth";
	if ( $verbose == 1 ) { print "$cmd\n"; }
	system("$cmd -te $treefname -redo $threads > opt_biophys_temp.screen");
	
	# Extract the log likelihood 
	system("grep \"^Log-likelihood of the tree:\" opt_biophys_temp.iqtree  > opt_biophys_temp.lnL");
	open (my $OPF, "opt_biophys_temp.lnL") or die "Could not open file opt_biophys_temp.lnL for input\n";
	my @iqfline = <$OPF>; # Read the temporary file with the lnL
	my @iqvals = split(/\s+/, $iqfline[0]);
	my($like) = $iqvals[4];
	close($OPF) or die "Could not close file opt_biophys_temp.lnL\n";
	
	return -1.0*$like;
}

