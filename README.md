# sixparam

### Six-parameter model of protein evolution ###

This goal of this project is estimate parameters for the six-parameter model of protein
evolution proposed in:

Braun EL (submitted) An evolutionary model motivated by physicochemical properties of 
amino acids reveals variation among proteins.


################
### Programs ###

There are three programs necessary to estimate the model parameters:

	--	1. optimize-biophysical-models.pl --
Perl program that reads a control file and optimizes a set of submodels based on eq3 and
eq4 from the manuscript. 

This program expects a relaxed phylip format protein datafile and a newick format treefile. This program calls a number of programs and requires some data files, including those listed below. The path to those files must be saved in the program. For example:

(the path to IQ-TREE is on lines 21-23)
	
	my($iqexec) = "iqtree-omp-1.5.5"; # iqtree executable
	my($threads) = " -nt 2";          # number of threads is multithread iqtree is used
	# my($threads) = " ";             # remove comment to use with serial iqtree

(the path to the generate-biophysical-matrix program is saved in the calc_lnL subroutine on line 444)

	system("./generate-biophysical-matrix $basemodel $genmat temp-biophys.dat $vol $pol $comp $arom $tv $gc 0");

To use the optimize-biophysical-models.pl program run as follows:

$ ./optimize-biophysical-models.pl <ctlfile> <outfile>
	
	Where:
		--  ctlfile = control file
		--  outfile = tab-delimited output file

	--	2. generate-biophysical-matrix.cpp --
C++ program that takes the model parameters as input and generates a PAML format R matrix
using those parameters. 

To compile using the GNU compiler simply use:
 $ g++ -c generate-biophysical-matrix.cpp
 $ g++ -o generate-biophysical-matrix generate-biophysical-matrix.o

	--	3. IQ-TREE: Download this program from http://www.iqtree.org --
The sixparam code has been tested with versions 1.5.5 and 1.5.6


#######################################
### Brief description of the models ###

The eq3 models have six parameters that reflect weights on changes in four different
amino acid properties:
		V = molecular volume
			P = polarity
		C = composition (atomic weight ratio of hetero to carbon atoms)
		A = aromaticity (PC III from Sneath 1966)
The other two parameters capture mutational input:
		G = weights the minimum # of nucleotide changes necessary for an amino acid interchange
		T = weights the impact of transversions
	
eq3 models are based on this equation:

	Rij = exp(-VdeltaVij) exp(-PdeltaP) exp(-CdeltaC) exp(-AdeltaA) exp(-TdeltaT) 1/Nij^G

Rij is the relevant element of the R matrix. The "delta" values are normalized differences 
between amino acid i and j in each amino acid property. deltaT is 1 if interchanging amino
acid i and j requires at least one transversion; otherwise it is 0. Nij is the minimum 
number of nucleotide changes necessary for an interchange of amino acids i and j.

The eq4 models add information from an empirical model of amino acid substitution, as
shown in this equation:

	Rij = Kij exp(-VdeltaVij) exp(-PdeltaP) exp(-CdeltaC) exp(-AdeltaA) exp(-TdeltaT) 1/Nij^G

Kij is the relevant element from an empirical model, such as the LG (Le and Gascuel 2008)
or JTT (Jones et al. 1992) models.


#####################
### Required data ###

IQ-TREE tests 18 empirical models (listed below). PAML format data matrices with the data
for these models are saved in the "empirical_models" directory

	Blosum62	HIVw		mtREV
	cpREV		JTT		mtZOA
	Dayhoff		JTTDCMut	PMB
	DCMut		LG		rtREV
	FLU		mtART		VT
	HIVb		mtMAM		WAG

Matrices for the universal genetic code are also provided.


#################
### Test data ###

Test data are provided in the following directories:

-- test_data_Chen	(from Chen et al. 2015)
The data matrices are reduced to a subset of 30 taxa with dense sampling. Trees are
the ML trees for each protein alignment using the best-fitting empirical model.

-- test_data_Jarvis	(from Jarvis et al. 2014)
The data matrices are those that were judged to have no homology errors by Springer
and Gatesy (2017). In a few cases one taxon with a potential homology error has been
removed (based on Springer and Gatesy 2017). Trees are the ML trees for each protein
alignment using the best-fitting empirical model.

-- test_data_Rokas	(from Rokas et al. 2005)
The data matrices were exported from the nexus file available from Rokas. Trees are
the ML trees for each protein alignment using the best-fitting empirical model.

-- test_data_Wolf	(from Wolf et al. 2004)
Eight concatenated data matrices and two trees (Coelomata and Ecdysozoa topologies)

All test_data folders contain a subfolder with control files to run analyses of the test
datasets. The format for the control files is:

	infile.phy	treefile.tre	poiss	universal_code	ml	gamma
	0	0	0	0	0	0	poiss
	1	0	0	0	0	0	V
	0	1	0	0	0	0	P
	etc.
	1	1	1	1	1	1	VPCATG

-or-

	infile.phy	treefile.tre	empirical_models/MODEL.dat	universal_code	ml	gamma
	0	0	0	0	0	0	MODEL
	1	0	0	0	0	0	MODEL+V
	0	1	0	0	0	0	MODEL+P
	etc.
	1	1	1	1	1	1	MODEL+VPCATG

For the input files the path should be specified unless the files are in the same
directory as optimize-biophysical-models.pl

For the later rows, 0 indicates that the parameter is not optimized (it is fixed at a 
value of zero) whereas a 1 indicates the model parameter is optimized. The order of the 
parameters is V, P, C, A, T, and G. The seventh column is simply the model name. "poiss" 
is used for the R matrix populated with 1's (i.e., the F81-like model)


##################
### References ###

Chen,M.Y. et al. (2015). Selecting question-specific genes to reduce incongruence in 
phylogenomics: a case study of jawed vertebrate backbone phylog-eny. Syst. Biol., 64, 
1104-1120.

Jarvis,E.D. et al. (2014) Whole-genome analyses resolve early branches in the tree of 
life of modern birds. Science, 346, 1320-1331.

Jones,D.T., Taylor,W.R., & Thornton,J.M. (1992) The rapid generation of mutation data 
matrices from protein sequences. CABIOS, 8, 275-282.

Le,S.Q. and Gascuel,O. (2008) An improved general amino acid replacement matrix. Mol. 
Biol. Evol., 25, 1307-1320.

Rokas,A. & Carroll,S.B. (2005). More genes or more taxa? The relative contribution of
gene number and taxon number to phylogenetic accuracy. Mol. Biol. Evol., 22, 1337-1344.

Sneath,P.H.A. (1966) Relations between chemical structure and biological activity in 
peptides. J. Theor. Biol., 12, 157-195.

Wolf,Y.I. et al. (2004) Coelomata and not Ecdysozoa: evidence from genome-wide 
phylogenetic analysis. Genome Res., 14, 29-36.
