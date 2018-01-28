#include <iostream>
#include <fstream>
#include <cstring>
#include <iomanip>
#include <cstdlib>
#include <cmath>
using namespace std;

/*
 *  generate_matrix.cpp v 0.9b
 *  
 *
 *  Created by Edward Braun on Oct 2017.
 *  Copyright 2017. All rights reserved.
 *
 */

class matrixmod {
	double pam_matrix[20][20];
	double aa_freq[20];

	double vol_diff[20][20];	// volume difference
	double pol_diff[20][20];	// polarity difference
	double comp_diff[20][20];   // compositional difference
	double arom_diff[20][20];   // aromaticity difference
	double tv_matrix[20][20];     // transversion difference given genetic code
	double gencode_matrix[20][20];  // genetic code matrix
	
public:

/******************************************************/
/* The following routines set (read or generate) the  */
/* PAM and gencode matrices for various analyses.     */
/******************************************************/
	
	/* Read a matrix for modification */
	void read_matrix(bool readfile, char* infn) {
		unsigned i, j;
		if (readfile) {
			ifstream pamf(infn);
			for (i=1; i<20; i++) {
				for (j=0; j<i; j++) {
					pamf >> pam_matrix[i][j];
					pam_matrix[j][i]=pam_matrix[i][j];
				}
			}
			for (i=0; i<20; i++) {
				pamf >> aa_freq[i];
				aa_freq[i]=aa_freq[i];
			}
			pamf.close();
		} // end if (readfile...
		else {
			for (i=0; i<20; i++) {
				for (j=0; j<20; j++) {
					pam_matrix[i][j]=1;
				}
				aa_freq[i]=0.05;
			}
		}
	} // end ** read_matrix()
	
	/* Fill a gencode matrix with ones, use the universal code, or read a genetic code matrix from a file */
	void set_gencode(bool readfile, char* infn) {
		unsigned i, j;
		if (readfile) {
			ifstream gcf(infn);
			for (i=1; i<20; i++) {
				for (j=0; j<i; j++) {
					gcf >> gencode_matrix[i][j];
					gencode_matrix[j][i]=gencode_matrix[i][j];
				}
			}
			gcf.close();
		} // end if (readfile...
		else {
			for (i=0; i<20; i++) {
				for (j=0; j<20; j++) {
					gencode_matrix[i][j]=1;
				}
			}
		} // end else...
	} // end ** read_gencode()
	
	/* Fill a transversion matrix with zeros or read a genetic code matrix from a file */
	void set_transversion(bool readfile, char* infn) {
		unsigned i, j;
		if (readfile) {
			ifstream tvcf(infn);
			for (i=1; i<20; i++) {
				for (j=0; j<i; j++) {
					tvcf >> tv_matrix[i][j];
					tv_matrix[j][i]=tv_matrix[i][j];
				}
			}
			tvcf.close();
		} // end if (readfile...
		else {
			for (i=0; i<20; i++) {
				for (j=0; j<20; j++) {
					tv_matrix[i][j]=0;
				}
			}
		} // end else...
	} // end ** set_transversion()


/******************************************************/
/* The following routines set the difference matrices */
/* for various amino acid characteristics.            */
/******************************************************/

	/* Set up the side chain volume (size) difference matrix */
    void set_volume() {
		unsigned i, j;
		double vol[20]; // vol stores volume values from Grantham (1974) Science 185:862-864
		vol[0]=31;		vol[1]=124;		vol[2]=56;		vol[3]=54;		vol[4]=55;
		vol[5]=85;		vol[6]=83;		vol[7]=3;		vol[8]=96;		vol[9]=111;
		vol[10]=111;	vol[11]=119;	vol[12]=105;	vol[13]=132;	vol[14]=32.5;
		vol[15]=32;		vol[16]=61;		vol[17]=170;	vol[18]=136;	vol[19]=84;
		
		// now set vol_max and the vol_diff array
		double vol_max=0;
		for (i=0; i<20; i++) {
			for (j=i; j<20; j++) {
				vol_diff[i][j]=vol[i]-vol[j];
				if (vol_diff[i][j]<0) vol_diff[i][j]=-1*vol_diff[i][j];
				if (vol_diff[i][j]>vol_max) vol_max=vol_diff[i][j];
			}
		}
		
		// normalize the vol_diff array
		for (i=0; i<20; i++) {
			for (j=i; j<20; j++) {
				vol_diff[i][j]=vol_diff[i][j]/vol_max;
				vol_diff[j][i]=vol_diff[i][j];  // generate a symmetric and normalized matrix
			}
		}
	} // end ** set_volume()
	
	/* Set up the amino acid polarity difference matrix */
	void set_polarity() {
		unsigned i, j;
		double pol[20]; // pol stores polarity values from Grantham (1974) Science 185:862-864
		pol[0]=8.1;		pol[1]=10.5;	pol[2]=11.6;	pol[3]=13.0;	pol[4]=5.5;
		pol[5]=10.5;	pol[6]=12.3;	pol[7]=9.0;		pol[8]=10.4;	pol[9]=5.2;
		pol[10]=4.9;	pol[11]=11.3;   pol[12]=5.7;	pol[13]=5.2;	pol[14]=8.0;
		pol[15]=9.2;	pol[16]=8.6;	pol[17]=5.4;	pol[18]=6.2;	pol[19]=5.9;
		
		// now set pol_max and the pol_diff array
		double pol_max=0;
		for (i=0; i<20; i++) {
			for (j=i; j<20; j++) {
				pol_diff[i][j]=pol[i]-pol[j];
				if (pol_diff[i][j]<0) pol_diff[i][j]=-1*pol_diff[i][j];
				if (pol_diff[i][j]>pol_max) pol_max=pol_diff[i][j];
			}
		}
		
		// normalize the pol_diff array
		for (i=0; i<20; i++) {
			for (j=i; j<20; j++) {
				pol_diff[i][j]=pol_diff[i][j]/pol_max;
				pol_diff[j][i]=pol_diff[i][j];  // generate a symmetric and normalized matrix
			}
		}
	} // end ** set_polarity()
	
	/* Set up the amino acid composition difference matrix */
	void set_composition() {
		unsigned i, j;
		double comp[20]; // comp stores composition values from Grantham (1974) Science 185:862-864
		comp[0]=0.00;	comp[1]=0.65;	comp[2]=1.33;	comp[3]=1.38;	comp[4]=2.75;
		comp[5]=0.89;	comp[6]=0.92;	comp[7]=0.74;	comp[8]=0.58;	comp[9]=0.00;
		comp[10]=0.00;	comp[11]=0.33;	comp[12]=0.00;	comp[13]=0.00;	comp[14]=0.39;
		comp[15]=1.42;	comp[16]=0.71;	comp[17]=0.13;	comp[18]=0.20;	comp[19]=0.00;
		
		// now set comp_max and the comp_diff array
		double comp_max=0;
		for (i=0; i<20; i++) {
			for (j=i; j<20; j++) {
				comp_diff[i][j]=comp[i]-comp[j];
				if (comp_diff[i][j]<0) comp_diff[i][j]=-1*comp_diff[i][j];
				if (comp_diff[i][j]>comp_max) comp_max=comp_diff[i][j];
			}
		}
		
		// normalize the comp_diff array
		for (i=0; i<20; i++) {
			for (j=i; j<20; j++) {
				comp_diff[i][j]=comp_diff[i][j]/comp_max;
				comp_diff[j][i]=comp_diff[i][j];  // generate a symmetric and normalized matrix
			}
		}
	} // end ** set_composition()
	
	/* Set up the amino acid aromaticity difference matrix */
	void set_aromaticity() {
		unsigned i, j;
		double arom[20]; // arom stores aromaticity (PC III) values from Sneath (1966) J Theor Biol 12:157-195
		arom[0]=-0.110;		arom[1]=0.079;		arom[2]=-0.136;		arom[3]=-0.285;	arom[4]=-0.184;
		arom[5]=-0.067;		arom[6]=-0.246;		arom[7]=-0.073;		arom[8]=0.320;	arom[9]=0.001;
		arom[10]=-0.008;	arom[11]=0.049;		arom[12]=-0.041;	arom[13]=0.438;	arom[14]=-0.016;
		arom[15]=-0.153;	arom[16]=-0.208;	arom[17]=0.493;		arom[18]=0.381;	arom[19]=-0.155;
		
		// now set arom_max and the arom_diff array
		double arom_max=0;
		for (i=0; i<20; i++) {
			for (j=i; j<20; j++) {
				arom_diff[i][j]=arom[i]-arom[j];
				if (arom_diff[i][j]<0) arom_diff[i][j]=-1*arom_diff[i][j];
				if (arom_diff[i][j]>arom_max) arom_max=arom_diff[i][j];
			}
		}
		
		// normalize the arom_diff array
		for (i=0; i<20; i++) {
			for (j=i; j<20; j++) {
				arom_diff[i][j]=arom_diff[i][j]/arom_max;
				arom_diff[j][i]=arom_diff[i][j];  // generate a symmetric and normalized matrix
			}
		}
	} // end ** set_aromaticity()
		
	/* Set up the volume, polarity, composition and aromaticity difference matrices */
	void set_parameter_matrices() {
		set_volume();
		set_polarity();
		set_composition();
		set_aromaticity();
	} // end ** set_parameter_matrices()

/******************************************************/
/* The following routine generates the modified       */
/* matrix and writes it to outfn.                     */
/******************************************************/
	
	/* Generate and write a matrix to "outfn" using the volume, polarity, composition, */
	/* aromaticity and genetic code parameters that are passed to the procedure.       */
	void generate_matrix(char* outfn, double vol_param, double pol_param, double comp_param, 
						 double arom_param, double tv_param, double gencode, unsigned aafappend) {
		ofstream pamf(outfn);
		unsigned i, j;
		double VF, PF, CF, AF, TVF, GCF;
		for (i=1; i<20; i++) {
			for (j=0; j<i; j++) {
				VF=exp(-1.0*vol_param*vol_diff[i][j]);
				PF=exp(-1.0*pol_param*pol_diff[i][j]);
				CF=exp(-1.0*comp_param*comp_diff[i][j]);
				AF=exp(-1.0*arom_param*arom_diff[i][j]);
				TVF=exp(-1.0*tv_param*tv_matrix[i][j]);
				if (gencode==0) GCF=1;
				else GCF=1.0/pow(gencode_matrix[i][j],gencode);
				/* write the modified PAM matrix value */
				pamf << pam_matrix[i][j]*VF*PF*CF*AF*TVF*GCF*1000.0 << ' ';
			}
			pamf << endl;
		}
		pamf << endl;
		
		/* Write the amino acid frequencies if aafappend == 1 */
		if ( aafappend == 1 ) {
			for (i=0; i<20; i++) {
				pamf << aa_freq[i] << ' ';
			}
			pamf << endl;
		}
		pamf.close();
	} // end ** generate_matrix()
	
}; // end matrixmod class

int main(int argc, char* argv[]) { 

	if (argc!=11) {
		cerr << "usage: $ ./generate_matrix <pamfile> <gencodefile> <outfile> <vol> <pol> <comp> <arom> <tv> <gencode> <aaf>\n"; 
		cerr << "   aaf == 0 ... do not append amino acid frequencies\n";
		cerr << "   aaf == 1 ... do append amino acid frequencies\n";
		return 0;
	}
	
	char pamfn[256];
	bool usepam = true;
	strcpy(pamfn,argv[1]);
	if ( strcmp(pamfn,"poiss")==0 || strcmp(pamfn,"poisson")==0 ) usepam=false;
	
	char codefn[256];
	char tvfn[256];
	bool usecode = true;
	bool usetv = true;
	strcpy(codefn,argv[2]);
	if ( strcmp(codefn,"none")==0 ) {
		usecode=false;
		usetv=false;
	}
	strcat(codefn,".dat");
	strcpy(tvfn,argv[2]);
	strcat(tvfn,"_tv.dat");
	
	char outfname[256];
	strcpy(outfname,argv[3]);
	
	double vpar = atof(argv[4]);
	double ppar = atof(argv[5]);
	double cpar = atof(argv[6]);
	double apar = atof(argv[7]);
	double tvpar = atof(argv[8]);
	double gpar = atof(argv[9]);
	unsigned append = atoi(argv[10]);
	
	matrixmod MAT;
	
	MAT.read_matrix(usepam,pamfn);
	MAT.set_gencode(usecode,codefn);
	MAT.set_transversion(usecode,tvfn);
	MAT.set_parameter_matrices();
	MAT.generate_matrix(outfname,vpar,ppar,cpar,apar,tvpar,gpar,append);
	
	return 0;
}
