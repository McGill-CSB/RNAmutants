#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <getopt.h>
#include <ctype.h>
#include <time.h>
#include <float.h>
#include <string.h>
#include "reader_energy_params.h"
#include "energy_functions.h"
#include "RNAmutants.h"
#include "util.h"
#include "sampling.h"
#include "mfe_backtrack.h"
#include "energy_params_tables.h"
#include "constraints.h"
#include "hybrid.h"

#define NO_DEBUG_SAMPLING
#define LEFT_CHECK 0
#define RIGHT_CHECK 26

/* various variables */

const double __infinity__ = DBL_MAX;

/* declared in turner_functions.c */

extern double __temperature;
extern double __RT__;

/* variables used in verbose mode */

unsigned long long int total_memory = 0;
time_t starting_time = 0;
time_t starting_algo_time = 0;

/* internal loop lookup table (mutants motifs) */

unsigned int *****loop_tab;
unsigned int *****stack_tab;

/* sequence constraints (structure constraints are declared in a separate file) */

int *cst_tape=NULL;
int **cst_segment=NULL;

/* special modes */
static int warning_flag=0;
int number_of_sample=0;
int framesize=0;
double cutpoint=0;
int include_intermolecular_interactions=1;

/* dangle flag (used in zuker_functions.c) */
int dangle_used = 1;

/* Version */

void version() {
	printf("RNAmutants version 2.0.1 (01/03/2012)\n");
	printf("Copyrights (C) 2011 Jerome Waldispuhl\n\n");
	printf("Academic version - may not be used by or for commercial organizations!\n\n");
	printf("For commercial use contact:\n\n");
	printf("Jerome Waldispuhl\n");
	printf("School of Computer Science\n");
	printf("McGill University\n");
	printf("Montreal, QC, CANADA\n");
	printf("E-mail: jeromew@cs.mcgill.ca\n\n");
	printf("THIS SOFTWARE COMES WITH ABSOLUTELY NO WARRANTY! USE AT YOUR OWN RISK!\n");
	exit(EXIT_SUCCESS);
}

/* Usage */

void usage(const char *softname) {
	fprintf(stdout,"Usage : %s [COMMANDS]\n",softname);
	fprintf(stdout,"COMMANDS (base):\n");
	fprintf(stdout,"--file\t\t[-f] <file>\t\tSet input sequence from fasta file <file>.\n");
	fprintf(stdout,"--string\t[-s] <string>\t\tSet input sequence from string.\n");
	fprintf(stdout,"--library\t[-l] <directory>\tEnergy paramaters library path.\n");
	fprintf(stdout,"--constraints\t[-C]\t\t\tUse sequence and structure contraints stored in input file.\n");
	fprintf(stdout,"--hybrid\t[-H]\t\t\tActivate the RNA duplex version and define the cutpoint.\n");
	
	fprintf(stdout,"COMMANDS (model):\n");
	fprintf(stdout,"--mutations\t[-m] <integer>\t\tSet the maximal number of mutations.\n");
	fprintf(stdout,"--sample-number\t[-n] <integer>\t\tSet the number of samples.\n");
	fprintf(stdout,"--sample-file\t[-N] <file>\t\tSet the parametric sampling command.\n");
	fprintf(stdout,"--bias\t\t[-M] <filename>\t\tMutational bias (4x4 matrix).\n");
	
	fprintf(stdout,"COMMANDS (energy):\n");
	fprintf(stdout,"--turner99\t\t\t\tUse 1999 Turner's nearest neighbor energy model.\n");
	fprintf(stdout,"\t\t\t\t\tAvailable only for a temperature of 37C.\n");
	fprintf(stdout,"--energy-model\t[-e] <directory>\tUse energy parameters stored in <directory>.\n");
	fprintf(stdout,"--temperature\t[-t] <integer>\t\tSet the temperature of energy parameters.\n");
	fprintf(stdout,"\t\t\t\t\tMust range from 0 to 100. Default is 37.\n");
	fprintf(stdout,"--absolute-temperature\t[-T] <number>\tFormal temperature in Kelvin degrees.\n");
	
	fprintf(stdout,"COMMANDS (utils):\n");
	fprintf(stdout,"--job-id\t[-J] <string>\t\tReport job ID <name> in error messages.\n");
	fprintf(stdout,"--help\t\t[-h]\t\t\tDisplay this menu.\n");
	fprintf(stdout,"--verbose\t[-v]\t\t\tEnable verbose mode.\n");
	fprintf(stdout,"--version\t[-V]\t\t\tShow version.\n");
	fprintf(stdout,"--warning\t[-w]\t\t\tEnable all warning.\n");
	
	fprintf(stdout,"COMMANDS (output):\n");
	fprintf(stdout,"--Z-output\t[-0] <file>\t\tSet the output file for the partition function values.\n");
	fprintf(stdout,"--mfe-output\t[-1] <file>\t\tSet the output file for the superoptimal mutants.\n");
	fprintf(stdout,"--sample-output\t[-2] <file>\t\tSet the output file for the samples.\n");
	fprintf(stdout,"--sample-stats\t\t\t\tOutput statistics computed on sample set.\n");	
	exit(EXIT_SUCCESS);
}

/* require by lexer */

int yywrap () {
	return 1;
}

/**********************************************************************************/
/*** utilities                                                                  ***/
/**********************************************************************************/

/*** counter for tri and tetraloop ***/

int init_triloop_cmpt_table() {
	int ii,jj,cmpt = 0;
	double tab[10] = {0,0,0,0,0,0,0,0,0,0};
	double val;
	
	for (ii=0;ii<1024;ii++) {
		if ((val=refTableTriLoop[ii])) {
			jj=0;
			while (tab[jj]) {
				if (val==tab[jj]) break; /* value already found */
				jj++;
				if ((jj==10)||(jj>cmpt)) {
					fprintf(stderr,"%s:line %d: index out of bound.\n",__FILE__,__LINE__);
					exit(EXIT_FAILURE);
				}
			}
			if (!tab[jj]) { tab[jj]=val; cmpt++; }
		}
	}
	
	triloop_cmpt_table = (double *)xmalloc((2*cmpt)*sizeof(double));
	triloop_weight_table = (double *)xmalloc((2*cmpt)*sizeof(double));
	memset(triloop_cmpt_table,0,(2*cmpt)*sizeof(double));
	memset(triloop_weight_table,0,(2*cmpt)*sizeof(double));
	for (ii=0;ii<cmpt;ii++) {
		if (!tab[ii]) {
			fprintf(stderr,"%s:line %d: Unexpected zero value.\n",__FILE__,__LINE__);
			exit(EXIT_FAILURE);
		}
		triloop_cmpt_table[2*ii] = tab[ii];
		triloop_weight_table[2*ii] = tab[ii];
	}
	
	nb_value_triloop_table = cmpt;
	
	return cmpt;
}

int init_tetraloop_cmpt_table() {
	int ii,jj,cmpt = 0;
	double tab[10] = {0,0,0,0,0,0,0,0,0,0};
	double val;
	
	for (ii=0;ii<4096;ii++) {
		if ((val=refTableTetraLoop[ii])) {
			jj=0;
			while (tab[jj]) {
				if (val==tab[jj]) break; /* value already found */
				jj++;
				if ((jj==10)||(jj>cmpt)) {
					fprintf(stderr,"%s:line %d: index out of bound.\n",__FILE__,__LINE__);
					exit(EXIT_FAILURE);
				}
			}
			if (!tab[jj]) { tab[jj]=val; cmpt++; }
		}
	}
	
	tetraloop_cmpt_table = (double *)xmalloc((2*cmpt)*sizeof(double));
	tetraloop_weight_table = (double *)xmalloc((2*cmpt)*sizeof(double));
	memset(tetraloop_cmpt_table,0,(2*cmpt)*sizeof(double));
	memset(tetraloop_weight_table,0,(2*cmpt)*sizeof(double));
	for (ii=0;ii<cmpt;ii++) {
		if (!tab[ii]) {
			fprintf(stderr,"%s:line %d: Unexpected zero value.\n",__FILE__,__LINE__);
			exit(EXIT_FAILURE);
		}
		tetraloop_cmpt_table[2*ii] = tab[ii];
		tetraloop_weight_table[2*ii] = tab[ii];
	}
	
	nb_value_tetraloop_table = cmpt;
	
	return cmpt;
}

int init_input_tape_from_string(const char *s){
	
	/* read arn tape, store it in the dynamic table and returns its length */
	
	char ch, *buffer;
	int i,n;
	
	/*----- Parse RNA sequence and place in input_tape, with len in rna_len --*/
	
	n = strlen(s);
	input_tape = (char *) calloc(n+1,sizeof(char));
	strcpy(input_tape,s);
	n=0; 
	for (i=0;i<strlen(input_tape);i++){
		ch = input_tape[i];
		switch(ch){
			case 'A':
			case 'a':
				input_tape[n++]='A';
				continue;
			case 'C':	
			case 'c':
				input_tape[n++]='C';
				continue;
			case 'G':
			case 'g':
				input_tape[n++]='G';
				continue;
			case 'U':	
			case 'u':
			case 'T':
			case 't':
				input_tape[n++]='U';
				continue;
			case '&':
			case '+':
				cutpoint=n-0.5; /* +0.5 for the cutpoint, but -1 b/c index start at 0 */
				continue; 
			default: continue;
		}
	}
	input_tape[n]='\0'; /* move string terminator at end */
	rna_len = n;

	/* empty constraints */
	
	buffer = (char *)xmalloc((rna_len+1)*sizeof(char));
	
	for (i=0;i<rna_len;i++) { buffer[i]='?'; } /* empty constraints */
	buffer[rna_len]='\0';
	init_structure_cst_tables(rna_len);
	fill_structure_cst_tables(buffer);

	for (i=0;i<rna_len;i++) { buffer[i]='o'; } /* empty constraints */
	buffer[rna_len]='\0';
	cst_tape = init_mutation_cst_tables(buffer);
	cst_segment = init_mutation_cst_segment(cst_tape);

	free(buffer);
	
	return rna_len; /* return length of RNA string */
}

int init_input_tape_from_file(const char *filename,int cst_flag) { /* read arn tape, store it in the dynamic table and returns its length */
	FILE *fp;
	char carac='\0', *buffer;
	int ii,index, canbe_seq,canbe_ss,canbe_mut,is_comment, seq_setup, ss_setup, mut_setup;
	
	/* read and store the input sequence */
	
	if (!filename) return -1; /* return empty constraint */
	
	if (!(fp=fopen(filename,"r"))) {
		fprintf(stderr,"Cannot open file %s.%s\n",filename, getJobID());
		return -1;
	}
	
	/* allocate buffer memory */
	
	buffer = (char *)xmalloc(max_sequence_size*sizeof(char));
	memset(buffer,0,max_sequence_size*sizeof(char));
	
	index=0;
	canbe_seq = 1;
	canbe_ss = 1;	
	canbe_mut = 1;
	seq_setup = 0;
	ss_setup = 0;
	mut_setup = 0;
	is_comment = 0;
	while ((carac = fgetc(fp)) != EOF) {
		
		/* skip comments */
		if ((is_comment)&&(carac!='\n')) {
			continue;
		}
		
		/* read */
		switch(carac){
			case ' ':
			case '\t':
				continue;
			case '>':
				is_comment = 1;
				continue;
			case '\r':
			case '\n':
				buffer[index]='\0';
				if ((canbe_seq)&&(!canbe_ss)&&(!canbe_mut)) {
					if (strlen(buffer)>0) {
						init_input_tape_from_string(buffer);
						seq_setup=1;
					}
				}
				else if ((!canbe_seq)&&(canbe_ss)&&(!canbe_mut)) { /* structure contraints */
					if (rna_len==index) {
						if (!cst_flag) {
							fprintf(stderr,"WARNING: Structure constraints ignored. Please use --constraints.%s\n", getJobID());
							for (ii=0;ii<rna_len;ii++) { buffer[ii]='?'; } /* empty constraints */
							buffer[rna_len]='\0';
						}
						init_structure_cst_tables(rna_len);
						fill_structure_cst_tables(buffer);
						ss_setup=1;
					}
					else {
						fprintf(stderr,"WARNING: Length of structure constraints do not match sequence size. Ignored.%s\n",getJobID());
					}
				}
				else if ((!canbe_seq)&&(!canbe_ss)&&(canbe_mut)) { /* sequence constraints */
					if (rna_len==index) {
						if (!cst_flag) {
							fprintf(stderr,"WARNING: Sequence constraints ignored. Please use --constraints.%s\n", getJobID());
							for (ii=0;ii<rna_len;ii++) { buffer[ii]='o'; } /* empty constraints */
							buffer[rna_len]='\0';
						}
						cst_tape = init_mutation_cst_tables(buffer);
						cst_segment = init_mutation_cst_segment(cst_tape);
						mut_setup=1;
					}
					else {
						fprintf(stderr,"WARNING: Length of sequence constraints do not match sequence size. Ignored.%s\n", getJobID());
					}
				}
				else if ((!is_comment)&&(strlen(buffer)>0)) {
					fprintf(stderr,"Corrupted input file %s.%s\n",filename,getJobID());
					exit(EXIT_FAILURE);
				}
				canbe_seq = 1;
				canbe_ss = 1;	
				canbe_mut = 1;
				is_comment = 0;
				memset(buffer,0,max_sequence_size*sizeof(char));
				buffer[0]='\0';
				index=0;
				continue;
			case 'A':
			case 'a':
			case 'C':
			case 'c':
			case 'G':
			case 'g':
			case 'U':
			case 'u':
			case 'T':
			case 't':
				buffer[index]=toupper(carac);
				index++;
				canbe_mut=0;
				canbe_ss=0;
				continue;
			case '&':
			case '+':
				if ((canbe_seq)&&(!canbe_ss)&&(!canbe_mut)) {
					buffer[index]=carac;
					index++;
				}
				else if ((cst_flag)&&((int)(cutpoint)!=index-1)) {
					fprintf(stderr,"ERROR: Corrupted input. Cut point mislocated.%s\n",getJobID());
				}
				continue;
			case '(':
			case ')':
			case '.':
			case '?':
				buffer[index]=carac; 
				index++;
				canbe_seq=0;
				canbe_mut=0;
				continue;
			case 'o':
			case 'O':
			case 'x':
			case 'X':
				buffer[index]=carac; 
				index++;
				canbe_seq=0;
				canbe_ss=0;
				continue;
			default:
				break;
		}
		
		if (index>=max_sequence_size) {
			fprintf(stderr,"Error: input exceed storage limit.%s\n",getJobID());
		}
	}

	/* last check in case of file does not end with a breakline */
	
	if (index) {
		buffer[index]='\0';
		if ((canbe_seq)&&(!canbe_ss)&&(!canbe_mut)) {
			if (strlen(buffer)>0) {
				init_input_tape_from_string(buffer);
				seq_setup=1;
			}
		}
		else if ((!canbe_seq)&&(canbe_ss)&&(!canbe_mut)) { /* structure contraints */
			if (rna_len==index) {
				if (!cst_flag) {
					fprintf(stderr,"WARNING: Structure constraints ignored. Please use --constraints.%s\n", getJobID());
					for (ii=0;ii<rna_len;ii++) { buffer[ii]='?'; } /* empty constraints */
					buffer[rna_len]='\0';
				}
				init_structure_cst_tables(rna_len);
				fill_structure_cst_tables(buffer);
				ss_setup=1;
			}
			else {
				fprintf(stderr,"WARNING: Length of structure constraints do not match sequence size. Ignored.%s\n",getJobID());
			}
		}
		else if ((!canbe_seq)&&(!canbe_ss)&&(canbe_mut)) { /* sequence constraints */
			if (rna_len==index) {
				if (!cst_flag) {
					fprintf(stderr,"WARNING: Sequence constraints ignored. Please use --constraints.%s\n", getJobID());
					for (ii=0;ii<rna_len;ii++) { buffer[ii]='o'; } /* empty constraints */
					buffer[rna_len]='\0';
				}
				cst_tape = init_mutation_cst_tables(buffer);
				cst_segment = init_mutation_cst_segment(cst_tape);
				mut_setup=1;
			}
			else {
				fprintf(stderr,"WARNING: Length of sequence constraints do not match sequence size. Ignored.%s\n", getJobID());
			}
		}
		else if ((!is_comment)&&(strlen(buffer)>0)) {
			fprintf(stderr,"Corrupted input file %s.%s\n",filename,getJobID());
			exit(EXIT_FAILURE);
		}
		canbe_seq = 1;
		canbe_ss = 1;	
		canbe_mut = 1;
		is_comment = 0;
		memset(buffer,0,max_sequence_size*sizeof(char));
		buffer[0]='\0';
		index=0;
	}
	
	/* check empty sequence */
	if ((!rna_len)||(!seq_setup)) {
		fprintf(stderr,"\nCannot find input sequence: Abort.%s\n\n", getJobID());
		usage("RNAmutants");
	}
	
	/* check and complete all assignments are done */
	
	if (!ss_setup) {
		if (cst_flag) {
			fprintf(stderr,"WARNING: Cannot find structure constraints. Ignored. %s\n", getJobID());
		}
		for (ii=0;ii<rna_len;ii++) { buffer[ii]='?'; } /* empty constraints */
		buffer[rna_len]='\0';
		init_structure_cst_tables(rna_len);
		fill_structure_cst_tables(buffer);
		ss_setup=1;
	}

	if (!mut_setup) {
		if (cst_flag) {
			fprintf(stderr,"WARNING: Cannot find sequence constraints. Ignored. %s\n", getJobID());
		}
		for (ii=0;ii<rna_len;ii++) { buffer[ii]='o'; } /* empty constraints */
		buffer[rna_len]='\0';
		cst_tape = init_mutation_cst_tables(buffer);
		cst_segment = init_mutation_cst_segment(cst_tape);
		mut_setup=1;
	}
	
	input_tape[rna_len] = '\0';
	free(buffer);
	fclose(fp);
	
	return rna_len;
}

/** read mutation matrix **/

int init_mutation_matrix_from_file(const char *filename) { 
	FILE *fp;
	int irow=0,icol=0, index=0,rcode=0,ii,jj;
	int list_aa_index[5] = {INDEX_A,INDEX_C,INDEX_G,INDEX_U,INDEX_STOP};
	double bias=-1.0;
	char *buffer;
	
	/* read and store the input sequence */
	
	if (!filename) {
		return -1; /* return empty constraint */
	}
	
	if (!(fp=fopen(filename,"r"))) {
		fprintf(stderr,"Cannot open file %s.%s\n",filename, getJobID());
		return -1;
	}

	buffer=(char *)xmalloc(100*sizeof(char));
	
	while ((rcode = fscanf(fp,"%s",buffer))!=EOF) {
		
		/* Checks */
		if ((rcode!=1)||(index>15)) {
			fprintf(stderr,"WARNING: Corrupted file %s. Ignored.%s\n",filename,getJobID());
			/** reset table **/
			for (ii=0;ii<4;ii++) {
				for (jj=0;jj<4;jj++) {
					mutbias[ii][jj]=1.0;
				}
			}
			free(buffer);
			return 0;
		}
				
		/* proceed */

		bias = atof(buffer);
		if (bias>=0) {
			mutbias[list_aa_index[irow]][list_aa_index[icol]] = bias;
		}
		else {
			fprintf(stderr,"ERROR: Corrupted mutation weigth %s. Ignored.%s\n",buffer,getJobID());
			/** reset table **/
			for (ii=0;ii<4;ii++) {
				for (jj=0;jj<4;jj++) {
					mutbias[ii][jj]=1.0;
				}
			}
			free(buffer);
			return -1;
		}
		
#if 0
		printf("M[%d,%d] <- %e (%d).\n",list_aa_index[irow],list_aa_index[icol],(double)(bias),index);
#endif
		index++;
		irow = index / 4;
		icol = index % 4;
		bias=-1.0;
		
	}

	if (index!=16) {
		fprintf(stderr,"WARNING: Incomplete mutation weight matrix in file \"%s\". Ignored.%s\n",filename,getJobID());
		/** reset table **/
		for (ii=0;ii<4;ii++) {
			for (jj=0;jj<4;jj++) {
				mutbias[ii][jj]=1.0;
			}
		}
		free(buffer);
		return 0;
	}

	fclose(fp);
	free(buffer);
	
	return 1;

}

/** build the list of possible configurations **/

unsigned int **build_generic_loop_list(int ileft, int iright, int nb_mutations) {
	int mutation_index[4] = {0,0,0,0}; 
	int ii,ind_ext=0,ind_int=0;
	int *liste_unbp_nt[2];
	int nb_mut_ext,nb_mut_int;
	int left_index_ext=char2index(ileft);
	int right_index_ext=char2index(iright);
	int left_index_int=char2index(ileft+1);
	int right_index_int=char2index(iright-1);
	int nt_left_ext,nt_right_ext,nt_left_int,nt_right_int;
	unsigned int hcode_config=0;
	unsigned int **list_of_config;
	int ii_list_of_config[3]={0,0,0};
	int all_config[5][7][4] = {{{0},{-1}},{{1},{2},{3},{4},{-1}},{{1,2},{1,3},{1,4},{2,3},{2,4},{3,4},{-1,-1}},
		{{1,2,3},{1,2,4},{1,3,4},{2,3,4},{-1,-1,-1}},{{1,2,3,4},{-1,-1,-1,-1}}};
	int i_config;
	
	/** allocate memory **/
	
	list_of_config=(unsigned int **)xmalloc(3*sizeof(unsigned int *));
	for(ii=0;ii<3;ii++) {
		list_of_config[ii]=(unsigned int *)xmalloc((int)(genereMutant(4,nb_mutations)+1.0)*sizeof(unsigned int));
	}
	
	/* enum configurations */
	
	i_config=0;
	while (all_config[nb_mutations][i_config][0]>=0) { /* enumerate possible distribution of mutation sites */
		
		/* init index mutation arrays */
		
		mutation_index[0]=mutation_index[1]=mutation_index[2]=mutation_index[3]=0;
		for (ii=0;ii<nb_mutations;ii++) {
			/* mutation_index[i] = 0 if WT and 1 if mutations at position i */
			mutation_index[all_config[nb_mutations][i_config][ii]-1]=1;
		}
		
		/*** 0: stack, 1: (i,j) base pair, 2: (i+1,j-1) base pair ***/
		
		/**************** stack ****************/
		
		/* 0: mutations, 1: right mutation, 2: left mutation, 3: mutations on both sides */
		nb_mut_ext=(mutation_index[0]<<1)|(mutation_index[3]);
		nb_mut_int=(mutation_index[1]<<1)|(mutation_index[2]);
		
		/* init exterior base pair */
		nt_left_ext = bp_mutation[nb_mut_ext][left_index_ext][right_index_ext][0][0];
		nt_right_ext = bp_mutation[nb_mut_ext][left_index_ext][right_index_ext][0][1];
		ind_ext = 1;
		
		while (nt_left_ext!=INDEX_STOP) {
			
			/* init interior base pair */
			nt_left_int = bp_mutation[nb_mut_int][left_index_int][right_index_int][0][0];
			nt_right_int = bp_mutation[nb_mut_int][left_index_int][right_index_int][0][1];
			ind_int=1;
			
			while (nt_left_int!=INDEX_STOP) {
				
				/****** REMOVE CONFIG WITH CST NT ******/
				
				if (!(((cst_tape[ileft])&&(kronecker(ileft,nt_left_ext)))||
					  ((cst_tape[iright])&&(kronecker(iright,nt_right_ext)))||
					  ((cst_tape[ileft+1])&&(kronecker(ileft+1,nt_left_int)))||
					  ((cst_tape[iright-1])&&(kronecker(iright-1,nt_right_int))))) {
					
					/* store hcode */
					hcode_config = (nt_left_ext<<6)|(nt_left_int<<4)|(nt_right_int<<2)|(nt_right_ext);
					list_of_config[0][ii_list_of_config[0]] = hcode_config;
					ii_list_of_config[0]++;
					
				}
				
				/* update opening base pair */
				
				nt_left_int = bp_mutation[nb_mut_int][left_index_int][right_index_int][ind_int][0];
				nt_right_int = bp_mutation[nb_mut_int][left_index_int][right_index_int][ind_int][1];
				ind_int++;
				
			}
			
			/* update closing base pair */
			
			nt_left_ext = bp_mutation[nb_mut_ext][left_index_ext][right_index_ext][ind_ext][0];
			nt_right_ext = bp_mutation[nb_mut_ext][left_index_ext][right_index_ext][ind_ext][1];
			ind_ext++;
		}
		
		/**************** only base pair exterior ****************/
		
		liste_unbp_nt[0] = uncst_mutation[mutation_index[1]][char2index(ileft+1)];
		liste_unbp_nt[1] = uncst_mutation[mutation_index[2]][char2index(iright-1)];
		
		/* init exterior base pair */
		nt_left_ext = bp_mutation[nb_mut_ext][left_index_ext][right_index_ext][0][0];
		nt_right_ext = bp_mutation[nb_mut_ext][left_index_ext][right_index_ext][0][1];
		ind_ext = 1;
		
		while (nt_left_ext!=INDEX_STOP) {
			int i_0=0;
			while ((nt_left_int = liste_unbp_nt[0][i_0])!=INDEX_STOP) {
				int i_1=0;
				while ((nt_right_int = liste_unbp_nt[1][i_1])!=INDEX_STOP) {
					
					/****** REMOVE CONFIG WITH CST NT ******/
					
					if (!(((cst_tape[ileft])&&(kronecker(ileft,nt_left_ext)))||
						  ((cst_tape[iright])&&(kronecker(iright,nt_right_ext)))||
						  ((cst_tape[ileft+1])&&(kronecker(ileft+1,nt_left_int)))||
						  ((cst_tape[iright-1])&&(kronecker(iright-1,nt_right_int))))) {
						
						/* store hcode */
						hcode_config = (nt_left_ext<<6)|(nt_left_int<<4)|(nt_right_int<<2)|(nt_right_ext);
						list_of_config[1][ii_list_of_config[1]] = hcode_config;
						ii_list_of_config[1]++;
						
					}
					i_1++;
				}
				i_0++;
			}
			
			/* update closing base pair */
			
			nt_left_ext = bp_mutation[nb_mut_ext][left_index_ext][right_index_ext][ind_ext][0];
			nt_right_ext = bp_mutation[nb_mut_ext][left_index_ext][right_index_ext][ind_ext][1];
			ind_ext++;
		}
		
		/**************** only base pair interior ****************/
		
		liste_unbp_nt[0] = uncst_mutation[mutation_index[0]][char2index(ileft)];
		liste_unbp_nt[1] = uncst_mutation[mutation_index[3]][char2index(iright)];
		
		/* init interior base pair */
		nt_left_int = bp_mutation[nb_mut_int][left_index_int][right_index_int][0][0];
		nt_right_int = bp_mutation[nb_mut_int][left_index_int][right_index_int][0][1];
		ind_int=1;
		
		while (nt_left_int!=INDEX_STOP) {
			int i_0=0;
			while ((nt_left_ext = liste_unbp_nt[0][i_0])!=INDEX_STOP) {
				int i_1=0;
				while ((nt_right_ext = liste_unbp_nt[1][i_1])!=INDEX_STOP) {
					
					/****** REMOVE CONFIG WITH CST NT ******/
					
					if (!(((cst_tape[ileft])&&(kronecker(ileft,nt_left_ext)))||
						  ((cst_tape[iright])&&(kronecker(iright,nt_right_ext)))||
						  ((cst_tape[ileft+1])&&(kronecker(ileft+1,nt_left_int)))||
						  ((cst_tape[iright-1])&&(kronecker(iright-1,nt_right_int))))) {
						
						/* store hcode */
						hcode_config = (nt_left_ext<<6)|(nt_left_int<<4)|(nt_right_int<<2)|(nt_right_ext);
						list_of_config[2][ii_list_of_config[2]] = hcode_config;
						ii_list_of_config[2]++;
						
					}
					i_1++;
				}
				i_0++;
			}
			
			
			/* update opening base pair */
			
			nt_left_int = bp_mutation[nb_mut_int][left_index_int][right_index_int][ind_int][0];
			nt_right_int = bp_mutation[nb_mut_int][left_index_int][right_index_int][ind_int][1];
			ind_int++;
			
		}
		
		/** ends of config **/
		
		i_config++;
		
	}
	
	/** setup the flags to mark the end of list **/
	
	for (ii=0;ii<3;ii++) {
		list_of_config[ii][ii_list_of_config[ii]] = 0; /* cannot be reached since AAAA is not a valid configuration */
	}
	
	
#if 0
	{
		int kk;
		printf("%d mutations\n",nb_mutations);
		for (kk=0;kk<3;kk++) {
			ii=0;
			printf("%d elements stored: ",ii_list_of_config[kk]);
			while(list_of_config[kk][ii]) {
				print_hcode(4,list_of_config[kk][ii]);
				printf(" ");
				ii++;
			}
			printf("\n");
		}
	}
#endif
	
	return list_of_config;
	
}

unsigned int ****enum_all_stacks(int nb_mutations) {
	int ii,jj;
	unsigned int ****tab_motifs;
    
	tab_motifs=(unsigned int ****)xmalloc((rna_len+1)*sizeof(unsigned int ***));
	for (ii=0;ii<=rna_len;ii++) {
		tab_motifs[ii]=(unsigned int ***)xmalloc((rna_len+1)*sizeof(unsigned int **)); 
		for (jj=0;jj<=rna_len;jj++) {
			tab_motifs[ii][jj]=NULL;
		}
	}
	
	for (ii=0;ii<=rna_len-4;ii++) {
		for (jj=ii+3;jj<rna_len;jj++) {
			
			if (cutpoint) {	/* cutpoint is used */
				
				if ((ii<cutpoint)&&(jj>cutpoint)) {
					if (include_intermolecular_interactions) {
						tab_motifs[ii][jj] = build_generic_loop_list(ii,jj,nb_mutations);
					}
					
				}
				else if (jj-ii>=__hairpin_min_length__+3){
					tab_motifs[ii][jj] = build_generic_loop_list(ii,jj,nb_mutations);
				}
				
			}
			else if (jj-ii>=__hairpin_min_length__+3){
				tab_motifs[ii][jj] = build_generic_loop_list(ii,jj,nb_mutations);
			}
		}
	}
	
	return tab_motifs;
	
}


/******* 1x1, 1x2, 2x1 and 2x2 internal loop cases ******/

int build_special_loop_list_for_given_config(unsigned int *list_of_config, int tabindex, int ileft, int iright, int nt_loop_left, int nt_loop_right, int nb_mutations,int *config) {
	int size_motif = nt_loop_left + nt_loop_right + 4; 
	int mutation_index[8] = {0,0,0,0,0,0,0,0}; 
	int ii,ind_close=0,ind_open=0, okconfig=0;
	int *liste_left_unbp_nt[2];
	int *liste_right_unbp_nt[2];
	int nb_mut_close,nb_mut_open;
	int left_index_close=char2index(ileft-1);
	int right_index_close=char2index(iright-1);
	int left_index_open=char2index(ileft+nt_loop_left);
	int right_index_open=char2index(iright-nt_loop_right-2);
	int nt_left_close,nt_right_close,nt_left_open,nt_right_open;
	unsigned int hcode_config=0;
	int shift, nt_lloop[2],nt_rloop[2],cmpt_lloop[2],cmpt_rloop[2],i_lloop,i_rloop;
	int ii_list_of_config=tabindex;
	
	/* init_arrays */
	
	for (ii=0;ii<nb_mutations;ii++) {
		mutation_index[config[ii]-1]=1;
	}
	
	nb_mut_close=(mutation_index[0]<<1)|(mutation_index[size_motif - 1]);
	nb_mut_open=(mutation_index[nt_loop_left+1]<<1)|(mutation_index[size_motif-nt_loop_right-2]);
	
	/**********************************/
	
	for (ii=0;ii<nt_loop_left;ii++) {
		liste_left_unbp_nt[ii] = uncst_mutation[mutation_index[ii+1]][char2index(ileft+ii)];
	}
	
	for (ii=0;ii<nt_loop_right;ii++) {
		liste_right_unbp_nt[ii] = uncst_mutation[mutation_index[size_motif-2-ii]][char2index(iright-2-ii)];
	}
	
	/* lookup mutated indices */
	
	/* init closing base pair */
	nt_left_close = bp_mutation[nb_mut_close][left_index_close][right_index_close][0][0];
	nt_right_close = bp_mutation[nb_mut_close][left_index_close][right_index_close][0][1];
	ind_close = 1;
	
	while (nt_left_close!=INDEX_STOP) {
		
		/* init opening base pair */
		nt_left_open = bp_mutation[nb_mut_open][left_index_open][right_index_open][0][0];
		nt_right_open = bp_mutation[nb_mut_open][left_index_open][right_index_open][0][1];
		ind_open=1;
		
		while (nt_left_open!=INDEX_STOP) {
			
			/*** enumerate nucleotides in loops ***/
			
			/** in left loop **/
			
			i_lloop=0;
			cmpt_lloop[0]=cmpt_lloop[1]=0;
			nt_lloop[0]=nt_lloop[1]=INDEX_GHOST;
			while (nt_lloop[0]!=INDEX_STOP) {
				
				if (i_lloop==nt_loop_left) {
					i_lloop--;
				}
				else {
					if ((nt_lloop[i_lloop]=liste_left_unbp_nt[i_lloop][cmpt_lloop[i_lloop]])!=INDEX_STOP) {
						
						cmpt_lloop[i_lloop]++;
						i_lloop++;
						
						/**  left loop completed **/
						if (i_lloop==nt_loop_left) {
							
							/** in right loop **/
							
							i_rloop=0;
							cmpt_rloop[0]=cmpt_rloop[1]=0;
							nt_rloop[0]=nt_rloop[1]=INDEX_GHOST;
							while (nt_rloop[0]!=INDEX_STOP) {
								
								if (i_rloop==nt_loop_right) {
									i_rloop--;
								}
								else {
									if ((nt_rloop[i_rloop]=liste_right_unbp_nt[i_rloop][cmpt_rloop[i_rloop]])!=INDEX_STOP) {
										cmpt_rloop[i_rloop]++;
										i_rloop++;
										
										/**  right loop completed **/
										if (i_rloop==nt_loop_right) {
											
											/****** REMOVE CONFIG WITH CST NT ******/
											
											okconfig=1;
											
											if (((cst_tape[ileft-1])&&(kronecker(ileft-1,nt_left_close)))||
												((cst_tape[iright-1])&&(kronecker(iright-1,nt_right_close)))||
												((cst_tape[ileft+nt_loop_left])&&(kronecker(ileft+nt_loop_left,nt_left_open)))||
												((cst_tape[iright-nt_loop_right-2])&&(kronecker(iright-nt_loop_right-2,nt_right_open)))) {
												okconfig=0;
											}
											
											for(ii=0;ii<nt_loop_left;ii++) {
												if ((cst_tape[ileft+ii])&&(kronecker(ileft+ii,nt_lloop[ii]))) {
													okconfig=0;
													break;
												}
											}
											
											for(ii=0;ii<nt_loop_right;ii++) {
												if ((cst_tape[iright-ii-2])&&(kronecker(iright-ii-2,nt_rloop[ii]))) {
													okconfig=0;
													break;
												}
											}
											
											
											/***************************************/
											
											if (okconfig) {
												
												/* compute hcode */
												
												hcode_config=nt_right_close;
												shift=2;
												for(ii=0;ii<nt_loop_right;ii++) {
													hcode_config |= (nt_rloop[ii])<<shift;
													shift+=2;
												}
												hcode_config |= nt_right_open<<shift;
												shift+=2;
												hcode_config |= nt_left_open<<shift;
												shift+=2;
												for(ii=nt_loop_left-1;ii>=0;ii--) {
													hcode_config |= (nt_lloop[ii])<<shift;
													shift+=2;
												}
												hcode_config |= nt_left_close<<shift;
												
												/* print hcode */
#if 0
												printf("(%d(%d,%d(%d,%d)%d,%d)%d) ---> %u\n",
													   nt_left_close,nt_lloop[0],nt_lloop[1],nt_left_open,nt_right_open,
													   nt_rloop[1],nt_rloop[0],nt_right_close,hcode_config);
#endif
												
												list_of_config[ii_list_of_config] = hcode_config;
												ii_list_of_config++;
												
												/**************************************/
											}
										}
									}
									else {
										cmpt_rloop[i_rloop]=0;
										i_rloop--;
									}
								} 
							}
						}
					}
					else {
						cmpt_lloop[i_lloop]=0;
						i_lloop--;
					}
				} 
			}
			
			/* update opening base pair */
			
			nt_left_open = bp_mutation[nb_mut_open][left_index_open][right_index_open][ind_open][0];
			nt_right_open = bp_mutation[nb_mut_open][left_index_open][right_index_open][ind_open][1];
			ind_open++;
			
		}
		
		/* update closing base pair */
		
		nt_left_close = bp_mutation[nb_mut_close][left_index_close][right_index_close][ind_close][0];
		nt_right_close = bp_mutation[nb_mut_close][left_index_close][right_index_close][ind_close][1];
		ind_close++;
	}
	
	/* end of list flag */
	list_of_config[ii_list_of_config] = 0;
	
	return ii_list_of_config;
	
}

unsigned int ****enum_all_special_internal_loops(int nb_mutations) {
	unsigned int ****tab_config;
	int ii,jj,cc,min_gap;
	static int size_tab = 9*sizeof(int), first_call = 1;
	static int *tab;

	if (cutpoint) { min_gap = 0; }
	else { min_gap = __hairpin_min_length__; }
	
	if (first_call) {
		tab = (int *)xmalloc(size_tab);
		first_call = 0;
	}
	
#if 0
	printf("start %d\n",nb_mutations);
#endif
	
	/** 0:1x1, 1:1x2, 2:2x1, 3:2x2 **/
	tab_config=(unsigned int ****)xmalloc(4*sizeof(unsigned int ***));
	
	for (cc=0;cc<4;cc++) {
		int index = 0;
		int lleft,lright,size_motif,**index_free=NULL;
		int max_value_1st;
		
		memset(tab,0,size_tab);
		lleft=1+(cc>>1);
		lright=1+(cc&1);
		size_motif = 4 + lleft + lright;
		max_value_1st = size_motif - nb_mutations + 1;
		
		/* we dont really need to optimize ii and jj */
		
		tab_config[cc]=(unsigned int ***)xmalloc(rna_len*sizeof(unsigned int **));
		index_free=(int **)xmalloc(rna_len*sizeof(int *));
		
		for (ii=0;ii<rna_len;ii++) {
			
			tab_config[cc][ii]=(unsigned int **)xmalloc(rna_len*sizeof(unsigned int *));
			index_free[ii]=(int *)xmalloc(rna_len*sizeof(int));
			
			for (jj=0;jj<rna_len;jj++) {
				
				if (jj-ii+1>=min_gap+size_motif-2) {
					tab_config[cc][ii][jj]=(unsigned int *)xmalloc((int)(genereMutant(size_motif,nb_mutations)+1)*sizeof(unsigned int));
				}
				else {
					tab_config[cc][ii][jj]=(unsigned int *)xmalloc(sizeof(unsigned int));
				}
				
				tab_config[cc][ii][jj][0] = 0; /* STOP signal: this config cannot be reached */
				index_free[ii][jj]=0; /* counter for the first available index */
				
			}
		}
		
		/* do NOT fill the table if number of mutations is too high */
		
		if (size_motif>=nb_mutations) {
			
			if (nb_mutations) {
				while(tab[0]<=max_value_1st) {
					
					if (tab[index]>size_motif) {
						tab[index]=0;
						index--;
					}
					else {
						if (index>=nb_mutations) {
							for (ii=0;ii<=rna_len-size_motif-min_gap;ii++) {
								for (jj=ii+size_motif+min_gap-1;jj<rna_len;jj++) {
									
									/* tab_config: list of configurations. 0 is stored at the last position. */
									/* index_free: index of the first available position (e.g. first cell with 0) */
									
									index_free[ii][jj] = build_special_loop_list_for_given_config(tab_config[cc][ii][jj],index_free[ii][jj],ii+1,jj+1,lleft,lright,nb_mutations,tab);
									tab_config[cc][ii][jj][index_free[ii][jj]]=0; /* STOP signal: this config (i.e. "0") cannot be reached */
									
								}
							}
							index--;
						} 
						else {
							tab[index]++;
							tab[index+1]=tab[index];
							index++;
						}
					}
				}
			}
			else {
				for (ii=0;ii<=rna_len-size_motif-min_gap;ii++) {
					for (jj=ii+size_motif+min_gap-1;jj<rna_len;jj++) {
						
						index_free[ii][jj] = build_special_loop_list_for_given_config(tab_config[cc][ii][jj],index_free[ii][jj],ii+1,jj+1,lleft,lright,0,NULL);
						tab_config[cc][ii][jj][index_free[ii][jj]]=0; /* STOP signal: this config (i.e. "0") cannot be reached */
						
					}
				}
			}
		}
		
		/* free local tables */
		for (ii=0;ii<rna_len;ii++) {
			free(index_free[ii]);
		}
		free(index_free);
	}
	
#if 0
	printf("%p\n",tab_config);
#endif
	
	/* exit */
	
	return tab_config;
	
}

/*********************************************************************/
/*                                                                   */
/* auxiliar functions                                                */
/*                                                                   */
/*********************************************************************/

/* init function used for Mutant algorithm */

void init_Mutant_arrays() { /* data structure can be optimized */
	int nn,ii,jj,kk,ll,local_len;
	
	Zs=(struct Zcell *****)xmalloc(max_mutations*sizeof(struct Zcell ****));
	Zes=(struct Zcell *****)xmalloc(max_mutations*sizeof(struct Zcell *****));
	Zms=(struct Zcell *****)xmalloc(max_mutations*sizeof(struct Zcell ****));
	Zm=(struct Zcell *****)xmalloc(max_mutations*sizeof(struct Zcell ****));
	Ze=(struct Zcell *****)xmalloc(max_mutations*sizeof(struct Zcell ****));
#ifdef DANGLE_ARRAY
	/* new arrays: internal loop and dangles */
	ZdangN=(struct Zcell *****)xmalloc(max_mutations*sizeof(struct Zcell ****));
	ZdangB=(struct Zcell *****)xmalloc(max_mutations*sizeof(struct Zcell ****));
	ZdangL=(struct Zcell *****)xmalloc(max_mutations*sizeof(struct Zcell ****));
	ZdangR=(struct Zcell *****)xmalloc(max_mutations*sizeof(struct Zcell ****));
#endif
	
	
	for (nn=0;nn<max_mutations;nn++) {
		
		Zs[nn]=(struct Zcell ****)xmalloc(rna_len*sizeof(struct Zcell ***));
		Zes[nn]=(struct Zcell ****)xmalloc(rna_len*sizeof(struct Zcell ***));
		Zms[nn]=(struct Zcell ****)xmalloc(rna_len*sizeof(struct Zcell ***));
		Zm[nn]=(struct Zcell ****)xmalloc(rna_len*sizeof(struct Zcell ***));
		Ze[nn]=(struct Zcell ****)xmalloc(rna_len*sizeof(struct Zcell ***));
#ifdef DANGLE_ARRAY
		ZdangN[nn]=(struct Zcell ****)xmalloc(rna_len*sizeof(struct Zcell ***));
		ZdangB[nn]=(struct Zcell ****)xmalloc(rna_len*sizeof(struct Zcell ***));
		ZdangL[nn]=(struct Zcell ****)xmalloc(rna_len*sizeof(struct Zcell ***));
		ZdangR[nn]=(struct Zcell ****)xmalloc(rna_len*sizeof(struct Zcell ***));
#endif
		
		for (ii=0;ii<rna_len;ii++) {
			
			Zs[nn][ii]=(struct Zcell ***)xmalloc(rna_len*sizeof(struct Zcell **));
			Zes[nn][ii]=(struct Zcell ***)xmalloc(rna_len*sizeof(struct Zcell **));
			Zms[nn][ii]=(struct Zcell ***)xmalloc(rna_len*sizeof(struct Zcell **));
			Zm[nn][ii]=(struct Zcell ***)xmalloc(rna_len*sizeof(struct Zcell **));
			Ze[nn][ii]=(struct Zcell ***)xmalloc(rna_len*sizeof(struct Zcell **));
#ifdef DANGLE_ARRAY
			ZdangN[nn][ii]=(struct Zcell ***)xmalloc(rna_len*sizeof(struct Zcell **));
			ZdangB[nn][ii]=(struct Zcell ***)xmalloc(rna_len*sizeof(struct Zcell **));
			ZdangL[nn][ii]=(struct Zcell ***)xmalloc(rna_len*sizeof(struct Zcell **));
			ZdangR[nn][ii]=(struct Zcell ***)xmalloc(rna_len*sizeof(struct Zcell **));
#endif
			
			for (jj=0;jj<rna_len;jj++) {
				local_len = jj - ii + 1;
				
				if (local_len < __hairpin_min_length__) { /* do not allocate memory for unexpected subsequences */
					Zs[nn][ii][jj]=0;
					Zes[nn][ii][jj]=0;
					Zms[nn][ii][jj]=0;
					Zm[nn][ii][jj]=0;
					Ze[nn][ii][jj]=0;
#ifdef DANGLE_ARRAY
					ZdangN[nn][ii][jj]=0;					
					ZdangB[nn][ii][jj]=0;
					ZdangL[nn][ii][jj]=0;
					ZdangR[nn][ii][jj]=0;
#endif
					continue;
				}
				
				Zs[nn][ii][jj]=(struct Zcell **)xmalloc(4*sizeof(struct Zcell *));
				Zes[nn][ii][jj]=(struct Zcell **)xmalloc(4*sizeof(struct Zcell *));
				Zms[nn][ii][jj]=(struct Zcell **)xmalloc(4*sizeof(struct Zcell *));
				Zm[nn][ii][jj]=(struct Zcell **)xmalloc(4*sizeof(struct Zcell *));
				Ze[nn][ii][jj]=(struct Zcell **)xmalloc(4*sizeof(struct Zcell *));
#ifdef DANGLE_ARRAY
				ZdangN[nn][ii][jj]=(struct Zcell **)xmalloc(4*sizeof(struct Zcell *));
				ZdangB[nn][ii][jj]=(struct Zcell **)xmalloc(4*sizeof(struct Zcell *));
				ZdangL[nn][ii][jj]=(struct Zcell **)xmalloc(4*sizeof(struct Zcell *));
				ZdangR[nn][ii][jj]=(struct Zcell **)xmalloc(4*sizeof(struct Zcell *));
#endif
				
				/* allocate memory for the visiblity fields */
				
				for (kk=0;kk<4;kk++) {
					
					Zs[nn][ii][jj][kk]=(struct Zcell *)xmalloc(4*sizeof(struct Zcell));
					Zes[nn][ii][jj][kk]=(struct Zcell *)xmalloc(4*sizeof(struct Zcell));
					Zms[nn][ii][jj][kk]=(struct Zcell *)xmalloc(4*sizeof(struct Zcell));
					Zm[nn][ii][jj][kk]=(struct Zcell *)xmalloc(4*sizeof(struct Zcell));
					Ze[nn][ii][jj][kk]=(struct Zcell *)xmalloc(4*sizeof(struct Zcell));
#ifdef DANGLE_ARRAY
					ZdangN[nn][ii][jj][kk]=(struct Zcell *)xmalloc(4*sizeof(struct Zcell));
					ZdangB[nn][ii][jj][kk]=(struct Zcell *)xmalloc(4*sizeof(struct Zcell));
					ZdangL[nn][ii][jj][kk]=(struct Zcell *)xmalloc(4*sizeof(struct Zcell));
					ZdangR[nn][ii][jj][kk]=(struct Zcell *)xmalloc(4*sizeof(struct Zcell));
#endif
					
					for (ll=0;ll<4;ll++) {
						Zs[nn][ii][jj][kk][ll].pf=0.0;
						Zs[nn][ii][jj][kk][ll].mfe=__infinity__;
						Zes[nn][ii][jj][kk][ll].pf=0.0;
						Zes[nn][ii][jj][kk][ll].mfe=__infinity__;
						Zms[nn][ii][jj][kk][ll].pf=0.0;
						Zms[nn][ii][jj][kk][ll].mfe=__infinity__;
						Zm[nn][ii][jj][kk][ll].pf=0.0;
						Zm[nn][ii][jj][kk][ll].mfe=__infinity__;
						Ze[nn][ii][jj][kk][ll].pf=0.0;	    
						Ze[nn][ii][jj][kk][ll].mfe=__infinity__;	    
#ifdef DANGLE_ARRAY
						ZdangN[nn][ii][jj][kk][ll].pf=0.0;	    
						ZdangN[nn][ii][jj][kk][ll].mfe=__infinity__;	    
						ZdangB[nn][ii][jj][kk][ll].pf=0.0;	    
						ZdangB[nn][ii][jj][kk][ll].mfe=__infinity__;	    
						ZdangL[nn][ii][jj][kk][ll].pf=0.0;	    
						ZdangL[nn][ii][jj][kk][ll].mfe=__infinity__;	    
						ZdangR[nn][ii][jj][kk][ll].pf=0.0;	    
						ZdangR[nn][ii][jj][kk][ll].mfe=__infinity__;	    
#endif
					}
				}
			}
		}
	}
}

/*************** precompute hairpin contribution (can also be integrate in main algo) ********************/

void precomputeHairpin() {
	int kk,ii,jj,xx,yy,uu,vv,ll,maxii,newK,nt2mut,ind;
	double bonus,n_total,n_partial,n_remaining,bound_bias;
	double hairpin_base_energy,mfe_hairpin_base_energy;
	
	for (kk=0;kk<max_mutations;kk++) {
		for (jj=__hairpin_min_length__+1;jj<rna_len;jj++) {
			maxii=jj-__hairpin_min_length__;
			for (ii=0;ii<maxii;ii++) {
				
				if (is_hairpin(ii,jj)) { /* lookup structure constraints */
					
					ll = jj-ii-1;
					nt2mut = ll-2 - cst_segment[ii+2][ll-3];
					
					for (xx=0;xx<4;xx++) {
						for (yy=0;yy<4;yy++) {
							if (validBasePair(xx,yy)) { /* special case when xx and yy base pair */
								
								/**** enumerate all interior base pair cases ****/
								
								Zs[kk][ii][jj][xx][yy].pf = 0.0;
								//Zs[kk][ii][jj][xx][yy].pfGC = 0.0;
								Zs[kk][ii][jj][xx][yy].mfe = __infinity__;
								
								for (uu=0;uu<4;uu++) {
									for (vv=0;vv<4;vv++) {
										
										hairpin_base_energy = boltzmannHairpin(xx,uu,vv,yy,ll);
										mfe_hairpin_base_energy = EHairpin(xx,uu,vv,yy,ll);
										newK=kk-kronecker(ii,xx)-kronecker(ii+1,uu)-kronecker(jj-1,vv)-kronecker(jj,yy);
										bound_bias =  nucleotide_bias(ii,xx) * nucleotide_bias(ii+1,uu) * nucleotide_bias(jj-1,vv) * nucleotide_bias(jj,yy);
										
										if ((newK>=0)&&(newK<=nt2mut)) {
											if ((__hairpin_min_length__<=ll)||(ll<=__hairpin_max_length__)) {
												/* check cst */
												if (!(((cst_tape[ii])&&(kronecker(ii,xx)))||
													  ((cst_tape[jj])&&(kronecker(jj,yy)))||
													  ((cst_tape[ii+1])&&(kronecker(ii+1,uu)))||
													  ((cst_tape[jj-1])&&(kronecker(jj-1,vv))))) {
													
													n_total = n_partial = 0;
													
													if (ll==3) { /** triloop **/
														triloop_counter(newK,ii,jj,xx,uu,vv,yy);
														for (ind=0;ind<nb_value_triloop_table;ind++) {
															/* n_partial = triloop_cmpt_table[2*ind+1]; */
															/* bonus = triloop_cmpt_table[2*ind]; */
															n_partial = triloop_weight_table[2*ind+1] * bound_bias;
															bonus = triloop_weight_table[2*ind];
															n_total += n_partial;
															if (n_partial>0) {
																Zs[kk][ii][jj][xx][yy].pf +=  hairpin_base_energy * exp(-(bonus/__RT__)) * n_partial;
																//Zs[kk][ii][jj][xx][yy].pfGC +=  hairpin_base_energy * exp(-(bonus/__RT__)) * n_partial;
																Zs[kk][ii][jj][xx][yy].mfe =  minimum_double(Zs[kk][ii][jj][xx][yy].mfe,mfe_hairpin_base_energy + bonus);
															}
														}
													} 
													else if (ll==4) { /** tetraloop **/ 
														tetraloop_counter(newK,ii,jj,xx,uu,vv,yy);
														n_total = 0;
														for (ind=0;ind<nb_value_tetraloop_table;ind++) {
															/* n_partial = tetraloop_cmpt_table[2*ind+1]; */
															/* bonus = tetraloop_cmpt_table[2*ind]; */
															n_partial = tetraloop_weight_table[2*ind+1] * bound_bias;
															bonus = tetraloop_weight_table[2*ind];
															n_total += n_partial;
															if (n_partial>0) {
																Zs[kk][ii][jj][xx][yy].pf +=  hairpin_base_energy * exp(-(bonus/__RT__)) * n_partial;
																Zs[kk][ii][jj][xx][yy].mfe =  minimum_double(Zs[kk][ii][jj][xx][yy].mfe,mfe_hairpin_base_energy + bonus);
															}
														}
													}
													
													/*** other cases ***/
#ifdef OLD_GGG_HAIRPIN
													n_partial = 0;
													if((xx==INDEX_G)&&(uu==INDEX_G)&&(yy==INDEX_U)) { /*** GGG loop ***/
														if (kronecker(ii+2,INDEX_G)) { /*** 3rd nt needs to mutate in G ***/
															if ((newK)&&(!(cst_tape[ii+2]))&&(nt2mut>0)) {
																n_partial = sequence_bias(ii+2,jj-1,newK-1) * bound_bias; 
																n_total += n_partial;
																if (n_partial>0) {
																	Zs[kk][ii][jj][xx][yy].pf +=  hairpin_base_energy * exp(-(ggg_hairpin_bonus/__RT__)) * n_partial;
																	Zs[kk][ii][jj][xx][yy].mfe =  minimum_double(Zs[kk][ii][jj][xx][yy].mfe,mfe_hairpin_base_energy + ggg_hairpin_bonus);
																}
															}
														}
														else { /*** already a GGG loop ***/
															if (nt2mut>0) {
																n_partial = sequence_bias(ii+2,jj-1,newK) * bound_bias;
																n_total += n_partial;
																if (n_partial>0) {
																	Zs[kk][ii][jj][xx][yy].pf +=  hairpin_base_energy * exp(-(ggg_hairpin_bonus/__RT__)) * n_partial;
																	Zs[kk][ii][jj][xx][yy].mfe =  minimum_double(Zs[kk][ii][jj][xx][yy].mfe,mfe_hairpin_base_energy + ggg_hairpin_bonus);
																}
															}
														}
													}
#endif
													
													/*** generic cases ***/
													
													n_remaining = sequence_bias(ii+2,jj-2,newK) * bound_bias - n_total;
													Zs[kk][ii][jj][xx][yy].pf +=  hairpin_base_energy  * n_remaining;
													Zs[kk][ii][jj][xx][yy].mfe =  minimum_double(Zs[kk][ii][jj][xx][yy].mfe,mfe_hairpin_base_energy);
    
                                                    /* update associated tables */
                                                    Zes[kk][ii][jj][xx][yy].pf  = Zms[kk][ii][jj][xx][yy].pf  = Zs[kk][ii][jj][xx][yy].pf;
                                                    Zes[kk][ii][jj][xx][yy].mfe = Zms[kk][ii][jj][xx][yy].mfe = Zs[kk][ii][jj][xx][yy].mfe;
													
												}
											}
										}
									}
								}
							}
						}
					}
				}	
			}
		}
	}
}


/*******************************************************************************/
/*                                                                             */
/* Mutant algorithm for computing of the boltzmann distribution                */
/*                                                                             */
/*******************************************************************************/

void RNAmutants_algorithm(int verbose) {
	
	int len, ii, jj, newii, newjj, rr, kk, xx, yy, nn, uu, vv, aa, bb, cc, dd, ee, ll, newK,nb_free_nt,nt2mut;
	unsigned int *list_config, *list_config_close, *list_config_open;
	int cloop,nb_mut_motif,nb_nt,nb_max_mut,ii_list,ii_list_open,ii_list_close,config,config_close,config_open;
	int nb_max_mut_unk_loop,nb_mut_unk_loop,nb_max_mut_close_bp,nb_mut_close_bp,nb_max_mut_open_bp,nb_mut_open_bp;
	int nb_max_mut_in_bulge, nb_mut_in_bulge, size_bulge,nb_mut_inside_boundary,min_rr,max_rr;
	double cmpt_mutants;
	clock_t start=0, current=0;
	double cpu_time_used;
	
	long totCount = ((rna_len-(__hairpin_min_length__+3))*(rna_len-(__hairpin_min_length__+2))/2)*((max_mutations)*(max_mutations-1))/2;
	long curCount = 0;
	
#ifdef DEBUG_SAMPLING
	double debug,pdebug;
#endif	
	
	if (verbose) {
		start = clock();
		current = start;
	}
	
	for (kk=0;kk<max_mutations;kk++) { /* number of mutations */
		
#ifdef INFO
		printf(".");fflush(stdout);
#endif


		/* SPECIAL CASE: triloops with an unpaired nucleotide on right side */
		
		len=__hairpin_min_length__+3;
		for (ii=0;ii<=rna_len-len;ii++) { /* ii is the index of the first nucleic acid of the subsequence (indices start at 0) */
			jj=ii+len-1;                  /* jj is the index of the last nucleic acid of the subsequence */
			
			for (xx=0;xx<4;xx++) {
				
				if (kronecker(ii,xx) && cst_tape[ii]) continue;
				
				for (yy=0;yy<4;yy++) {
					
					if (kronecker(jj,yy) && cst_tape[jj]) continue;
					
					/* right */
					if (!((cst_tape[jj])&&(kronecker(jj,yy)))) {
						for (uu=0;uu<4;uu++) {
							newK = kk - kronecker(jj,yy);
							if (newK>=0) {
								Zes[kk][ii][jj][xx][yy].pf += Zs[newK][ii][jj-1][xx][uu].pf * nucleotide_bias(jj,yy);
								Zes[kk][ii][jj][xx][yy].mfe = minimum_double(Zes[kk][ii][jj][xx][yy].mfe,Zs[newK][ii][jj-1][xx][uu].mfe);
								Zms[kk][ii][jj][xx][yy].pf += boltzmannExtendMultiLoop(1) * Zs[newK][ii][jj-1][xx][uu].pf * nucleotide_bias(jj,yy);
								Zms[kk][ii][jj][xx][yy].mfe = minimum_double(Zms[kk][ii][jj][xx][yy].mfe,EExtendMultiLoop(1) + Zs[newK][ii][jj-1][xx][uu].mfe);
                                /* FIXME: final delete
                                Ze[kk][ii][jj][xx][yy].pf += Ze[newK][ii][jj-1][xx][uu].pf * nucleotide_bias(jj,yy);
								Ze[kk][ii][jj][xx][yy].mfe = minimum_double(Ze[kk][ii][jj][xx][yy].mfe,Ze[newK][ii][jj-1][xx][uu].mfe);
								Zm[kk][ii][jj][xx][yy].pf += boltzmannExtendMultiLoop(1) * Zm[newK][ii][jj-1][xx][uu].pf * nucleotide_bias(jj,yy);
								Zm[kk][ii][jj][xx][yy].mfe = minimum_double(Zm[kk][ii][jj][xx][yy].mfe,EExtendMultiLoop(1) + Zm[newK][ii][jj-1][xx][uu].mfe);
                                */
							}
						}
					}
					
					/* left */
					if (!((cst_tape[ii])&&(kronecker(ii,xx)))) {
						for (uu=0;uu<4;uu++) {										
							newK = kk - kronecker(ii,xx);
							if (newK>=0) {
                                /* FIXME: Ze & Zm as above */
								Zes[kk][ii][jj][xx][yy].pf += Zs[newK][ii+1][jj][uu][yy].pf * boltzmann_penalty_close_helix(uu,yy) * nucleotide_bias(ii,xx);
								Zes[kk][ii][jj][xx][yy].mfe = minimum_double(Zes[kk][ii][jj][xx][yy].mfe,Zs[newK][ii+1][jj][uu][yy].mfe + penalty_close_helix(uu,yy));
								Zms[kk][ii][jj][xx][yy].pf += Zs[newK][ii+1][jj][uu][yy].pf * boltzmann_penalty_close_helix(uu,yy) * nucleotide_bias(ii,xx) * boltzmannExtendMultiLoop(1);
								Zms[kk][ii][jj][xx][yy].mfe = minimum_double(Zms[kk][ii][jj][xx][yy].mfe,Zs[newK][ii+1][jj][uu][yy].mfe + penalty_close_helix(uu,yy) + EExtendMultiLoop(1));
							}
						}
					}
				}
			}
		}
		
		/* GENERAL CASE */
		
		for (len=__hairpin_min_length__+4;len<=rna_len;len++) { /* len is the number of nucleic acids in the subsequence */
			/* minimum length is the minimal size for stack since hairpins have already been computed. */
			
			for (ii=0;ii<=rna_len-len;ii++) { /* ii is the index of the first nucleic acid of the subsequence (indices start at 0) */
				jj=ii+len-1;                  /* jj is the index of the last nucleic acid of the subsequence */
				curCount += kk;
				
				if (verbose) fprintf(stderr,"\r%4.1f%%    ",(double)(100.0*curCount)/(double)totCount);fflush(stderr);
				
#ifdef DEBUG_SAMPLING
				debug = pdebug = 0;
				if ((ii==LEFT_CHECK)&&(jj==RIGHT_CHECK)) {
					for (xx=0;xx<4;xx++) {
						for (yy=0;yy<4;yy++) {
							if (validBasePair(xx,yy)) { /* special case when xx and yy base pair */
								debug += Zs[kk][ii][jj][xx][yy].pf;
							}
						}
					}
					printf("hairpin: %e\n",debug);
					pdebug = debug;
				}
#endif

				if (is_basepair(ii,jj)) { /* NEW: check base pair to improve design performances */

					/* lookup stacks */
					
					nb_max_mut=minimum(4-cst_tape[ii]-cst_tape[ii+1]-cst_tape[jj-1]-cst_tape[jj],kk);
					
					if ((len>=__hairpin_min_length__+4)&&(is_stack(ii,jj))){ /* lookup structure constraints */
						
						for (nb_mut_motif=0;nb_mut_motif<=nb_max_mut;nb_mut_motif++) {
							list_config = stack_tab[nb_mut_motif][ii][jj][0];
							ii_list = 0;
							while ((config=list_config[ii_list])) {
								xx=config>>6;
								uu=(config>>4)&3;
								vv=(config>>2)&3;
								yy=config&3;
								newK=kk-nb_mut_motif + kronecker(ii+1,uu) + kronecker(jj-1,vv);
								Zs[kk][ii][jj][xx][yy].pf += boltzmannStack(xx,uu,vv,yy) * Zs[newK][ii+1][jj-1][uu][vv].pf * nucleotide_bias(ii,xx) * nucleotide_bias(jj,yy);
								Zs[kk][ii][jj][xx][yy].mfe = minimum_double(Zs[kk][ii][jj][xx][yy].mfe, EStack(xx,uu,vv,yy) + Zs[newK][ii+1][jj-1][uu][vv].mfe);
#ifdef DEBUG_SAMPLING
								if ((ii==LEFT_CHECK)&&(jj==RIGHT_CHECK)) debug += boltzmannStack(xx,uu,vv,yy) * Zs[newK][ii+1][jj-1][uu][vv].pf;
#endif
								ii_list++;
							}
						}
					}
					
#ifdef DEBUG_SAMPLING
					if ((ii==LEFT_CHECK)&&(jj==RIGHT_CHECK)) printf("stack: %e\n",debug-pdebug);
					pdebug = debug;
#endif
					
					/* lookup special internal loops */
					for (cloop=0;cloop<4;cloop++) {
						
						nb_nt = 6 + (cloop>>1) + (cloop&1); /* number of nt in the internal loop */
						nb_max_mut=minimum(nb_nt,kk);
						
						for (nb_mut_motif=0;nb_mut_motif<=nb_max_mut;nb_mut_motif++) {
							list_config = loop_tab[nb_mut_motif][cloop][ii][jj];
							ii_list = 0;
							while ((config=list_config[ii_list])) {
								switch (cloop) {
									case 0:
										if (!is_internal_loop(ii,ii+2,jj-2,jj)) break; /* lookup structure constraints */
										xx=config>>10;
										aa=(config>>8)&3;
										uu=(config>>6)&3;
										vv=(config>>4)&3;
										bb=(config>>2)&3;
										yy=config&3;
										newK = kk - nb_mut_motif + kronecker(ii+2,uu) + kronecker(jj-2,vv);
										Zs[kk][ii][jj][xx][yy].pf += boltzmannInternal_1x1(xx,aa,uu,vv,bb,yy) * Zs[newK][ii+2][jj-2][uu][vv].pf *
										nucleotide_bias(ii,xx) * nucleotide_bias(jj,yy) * nucleotide_bias(ii+1,aa) * nucleotide_bias(jj-1,bb);
										Zs[kk][ii][jj][xx][yy].mfe = minimum_double(Zs[kk][ii][jj][xx][yy].mfe,EInternal_1x1(xx,aa,uu,vv,bb,yy) + Zs[newK][ii+2][jj-2][uu][vv].mfe);  
#ifdef DEBUG_SAMPLING
										if ((ii==LEFT_CHECK)&&(jj==RIGHT_CHECK)) debug += boltzmannInternal_1x1(xx,aa,uu,vv,bb,yy) * Zs[newK][ii+2][jj-2][uu][vv].pf;
#endif
										break;
									case 1:
										if (!is_internal_loop(ii,ii+2,jj-3,jj)) break; /* lookup structure constraints */
										xx=config>>12;
										aa=(config>>10)&3;
										uu=(config>>8)&3;
										vv=(config>>6)&3;
										cc=(config>>4)&3;
										bb=(config>>2)&3;
										yy=config&3;
										newK = kk - nb_mut_motif + kronecker(ii+2,uu) + kronecker(jj-3,vv);
										Zs[kk][ii][jj][xx][yy].pf += boltzmannInternal_1x2(xx,aa,uu,vv,cc,bb,yy) * Zs[newK][ii+2][jj-3][uu][vv].pf *
										nucleotide_bias(ii,xx) * nucleotide_bias(jj,yy) * nucleotide_bias(ii+1,aa) * nucleotide_bias(jj-2,cc) * nucleotide_bias(jj-1,bb);
										Zs[kk][ii][jj][xx][yy].mfe = minimum_double(Zs[kk][ii][jj][xx][yy].mfe,EInternal_1x2(xx,aa,uu,vv,cc,bb,yy) + Zs[newK][ii+2][jj-3][uu][vv].mfe);
#ifdef DEBUG_SAMPLING
										if ((ii==LEFT_CHECK)&&(jj==RIGHT_CHECK)) debug += boltzmannInternal_1x2(xx,aa,uu,vv,cc,bb,yy) * Zs[newK][ii+2][jj-3][uu][vv].pf;
#endif
										break;
									case 2:
										if (!is_internal_loop(ii,ii+3,jj-2,jj)) break; /* lookup structure constraints */
										xx=config>>12;
										aa=(config>>10)&3;
										cc=(config>>8)&3;
										uu=(config>>6)&3;
										vv=(config>>4)&3;
										bb=(config>>2)&3;
										yy=config&3;
										newK = kk - nb_mut_motif + kronecker(ii+3,uu) + kronecker(jj-2,vv);
										Zs[kk][ii][jj][xx][yy].pf += boltzmannInternal_2x1(xx,aa,cc,uu,vv,bb,yy) * Zs[newK][ii+3][jj-2][uu][vv].pf *
										nucleotide_bias(ii,xx) * nucleotide_bias(jj,yy) * nucleotide_bias(ii+1,aa) * nucleotide_bias(ii+2,cc) * nucleotide_bias(jj-1,bb);
										Zs[kk][ii][jj][xx][yy].mfe = minimum_double(Zs[kk][ii][jj][xx][yy].mfe,EInternal_2x1(xx,aa,cc,uu,vv,bb,yy) + Zs[newK][ii+3][jj-2][uu][vv].mfe);
#ifdef DEBUG_SAMPLING
										if ((ii==LEFT_CHECK)&&(jj==RIGHT_CHECK)) debug += boltzmannInternal_2x1(xx,aa,cc,uu,vv,bb,yy) * Zs[newK][ii+3][jj-2][uu][vv].pf;
#endif
										break;
									case 3:
										if (!is_internal_loop(ii,ii+3,jj-3,jj)) break; /* lookup structure constraints */
										xx=config>>14;
										aa=(config>>12)&3;
										cc=(config>>10)&3;
										uu=(config>>8)&3;
										vv=(config>>6)&3;
										dd=(config>>4)&3;
										bb=(config>>2)&3;
										yy=config&3;
										newK = kk - nb_mut_motif + kronecker(ii+3,uu) + kronecker(jj-3,vv);
										Zs[kk][ii][jj][xx][yy].pf += boltzmannInternal_2x2(xx,aa,cc,uu,vv,dd,bb,yy) * Zs[newK][ii+3][jj-3][uu][vv].pf *
										nucleotide_bias(ii,xx) * nucleotide_bias(jj,yy) * nucleotide_bias(ii+1,aa) * nucleotide_bias(ii+2,cc) * nucleotide_bias(jj-2,dd) * nucleotide_bias(jj-1,bb);
										Zs[kk][ii][jj][xx][yy].mfe = minimum_double(Zs[kk][ii][jj][xx][yy].mfe,EInternal_2x2(xx,aa,cc,uu,vv,dd,bb,yy) + Zs[newK][ii+3][jj-3][uu][vv].mfe);
#ifdef DEBUG_SAMPLING
										if ((ii==LEFT_CHECK)&&(jj==RIGHT_CHECK)) debug += boltzmannInternal_2x2(xx,aa,cc,uu,vv,dd,bb,yy) * Zs[newK][ii+3][jj-3][uu][vv].pf;
#endif
										break;
								}
								ii_list++;
							}
						}
					}
					
#ifdef DEBUG_SAMPLING
					if ((ii==LEFT_CHECK)&&(jj==RIGHT_CHECK)) printf("special cases: %e\n",debug-pdebug);
					pdebug = debug;
#endif
					
					/** 1xn and nx1 internal loops **/
					
					nb_free_nt = minimum(__max_size_bulge__,len - __hairpin_min_length__ - 7);
					
					/** enumerate **/
					
					for (nn=1;nn<=nb_free_nt;nn++) { /* number of undefined nucleotides in internal loop */
						
						nb_max_mut_close_bp=minimum(4,kk);
						
						for (nb_mut_close_bp=0;nb_mut_close_bp<=nb_max_mut_close_bp;nb_mut_close_bp++) {
							
							list_config_close = stack_tab[nb_mut_close_bp][ii][jj][1];
							ii_list_close = 0;
							
							while ((config_close=list_config_close[ii_list_close])) {
								
								xx=config_close>>6;
								aa=(config_close>>4)&3;
								bb=(config_close>>2)&3;
								yy=config_close&3;
								
								/** 1xn **/
								
								if (is_internal_loop(ii,ii+2,jj-3-nn,jj)) { /* lookup structure constraints */
									
									nt2mut=nn-cst_segment[jj-1-nn][nn-1];
									nb_max_mut_unk_loop = minimum(nt2mut,kk-nb_mut_close_bp); /* to be adjusted later */
									
									for (nb_mut_unk_loop=0;nb_mut_unk_loop<=nb_max_mut_unk_loop;nb_mut_unk_loop++) { /* number of mutations in unknown part of loops */
										
										/** dont forget to add kronecker because mutation on single loop is be counted twice **/
										nb_max_mut_open_bp=minimum(4,kk-nb_mut_unk_loop-nb_mut_close_bp+kronecker(ii+1,aa)); 
										
										for (nb_mut_open_bp=0;nb_mut_open_bp<=nb_max_mut_open_bp;nb_mut_open_bp++) {
											
											list_config_open = stack_tab[nb_mut_open_bp][ii+1][jj-2-nn][2];
											ii_list_open = 0;
											while ((config_open=list_config_open[ii_list_open])) {
												
												cc=(config_open>>6)&3;
												
#if DEEP_DEBUG
												if (kk==0) {
													printf("--------------------------------\n");
													printf("Close:   %d-%d\n",xx,yy);
													printf("Close:  %d   %d\n",aa,bb);
													printf("Open:   %d   %d\n",(config_open>>6)&3,(config_open>>0)&3);
													printf("Open:    %d-%d\n",(config_open>>4)&3,(config_open>>2)&3);
												}
#endif
												
												if (aa==cc) {
													
													uu=(config_open>>4)&3;
													vv=(config_open>>2)&3;
													dd=config_open&3;
													
													newK = kk - nb_mut_close_bp - nb_mut_open_bp - nb_mut_unk_loop + kronecker(ii+1,aa) + kronecker(ii+2,uu) + kronecker(jj-3-nn,vv);
													
													if (newK>=0) {
														
														/* cmpt_mutants = genereMutant(nt2mut,nb_mut_unk_loop); */
														cmpt_mutants = sequence_bias(jj-1-nn,jj-2,nb_mut_unk_loop);
														
														if (cmpt_mutants>0) {
															Zs[kk][ii][jj][xx][yy].pf += boltzmannInternal_generic(xx,aa,cc,uu,vv,dd,bb,yy,1,nn+2) *
															Zs[newK][ii+2][jj-3-nn][uu][vv].pf * cmpt_mutants * nucleotide_bias(ii,xx) * nucleotide_bias(ii+1,aa) *
															nucleotide_bias(jj-2-nn,dd) * nucleotide_bias(jj-1,bb) * nucleotide_bias(jj,yy);
															Zs[kk][ii][jj][xx][yy].mfe = minimum_double(Zs[kk][ii][jj][xx][yy].mfe,
																										EInternal_generic(xx,aa,cc,uu,vv,dd,bb,yy,1,nn+2) + Zs[newK][ii+2][jj-3-nn][uu][vv].mfe);
#ifdef DEBUG_SAMPLING
															if ((ii==LEFT_CHECK)&&(jj==RIGHT_CHECK)) debug += boltzmannInternal_generic(xx,aa,cc,uu,vv,dd,bb,yy,1,nn+2) *
																Zs[newK][ii+2][jj-3-nn][uu][vv].pf * cmpt_mutants;
#endif
														}
													}
												}
												ii_list_open++;
											}
										}
									}
								}
								
								/** nx1 **/
								
								if (is_internal_loop(ii,ii+3+nn,jj-2,jj)) { /* lookup structure constraints */
									
									nt2mut=nn-cst_segment[ii+2][nn-1];
									nb_max_mut_unk_loop = minimum(nt2mut,kk-nb_mut_close_bp); /* to be adjusted later */
									
									for (nb_mut_unk_loop=0;nb_mut_unk_loop<=nb_max_mut_unk_loop;nb_mut_unk_loop++) { /* number of mutations in unknown part of loops */
										
										/** dont forget to add kronecker because mutation on single loop is be counted twice **/
										nb_max_mut_open_bp=minimum(4,kk-nb_mut_unk_loop-nb_mut_close_bp+kronecker(jj-1,bb)); 
										
										for (nb_mut_open_bp=0;nb_mut_open_bp<=nb_max_mut_open_bp;nb_mut_open_bp++) {
											
											list_config_open = stack_tab[nb_mut_open_bp][ii+2+nn][jj-1][2];
											ii_list_open = 0;
											while ((config_open=list_config_open[ii_list_open])) {
												
												dd=config_open&3;
												
												if (bb==dd) {
													
													cc=(config_open>>6)&3;
													uu=(config_open>>4)&3;
													vv=(config_open>>2)&3;
													
													newK = kk - nb_mut_close_bp - nb_mut_open_bp - nb_mut_unk_loop + kronecker(jj-1,bb) + kronecker(ii+3+nn,uu) + kronecker(jj-2,vv);
													if (newK>=0) {
														
														/* cmpt_mutants = genereMutant(nt2mut,nb_mut_unk_loop); */
														cmpt_mutants = sequence_bias(ii+2,ii+nn+1,nb_mut_unk_loop);
														
														if (cmpt_mutants>0) {
															Zs[kk][ii][jj][xx][yy].pf += boltzmannInternal_generic(xx,aa,cc,uu,vv,dd,bb,yy,nn+2,1) *
															Zs[newK][ii+3+nn][jj-2][uu][vv].pf * cmpt_mutants * nucleotide_bias(ii,xx) * nucleotide_bias(ii+1,aa) *
															nucleotide_bias(ii+2+nn,cc) * nucleotide_bias(jj-1,bb) * nucleotide_bias(jj,yy);
															Zs[kk][ii][jj][xx][yy].mfe = minimum_double(Zs[kk][ii][jj][xx][yy].mfe,
																										EInternal_generic(xx,aa,cc,uu,vv,dd,bb,yy,nn+2,1) + Zs[newK][ii+3+nn][jj-2][uu][vv].mfe);
#ifdef DEBUG_SAMPLING
															if ((ii==LEFT_CHECK)&&(jj==RIGHT_CHECK)) debug += boltzmannInternal_generic(xx,aa,cc,uu,vv,dd,bb,yy,nn+2,1) *
																Zs[newK][ii+3+nn][jj-2][uu][vv].pf * cmpt_mutants * nucleotide_bias(ii,xx) * nucleotide_bias(ii+1,aa) *
																nucleotide_bias(ii+2+nn,cc) * nucleotide_bias(jj-1,bb) * nucleotide_bias(jj,yy);
#endif
														}
													}
												}   
												ii_list_open++;
											}
										}
									}
								}
								ii_list_close++;
							}
						}
					}
					
#ifdef DEBUG_SAMPLING
					if ((ii==LEFT_CHECK)&&(jj==RIGHT_CHECK)) printf("1xn & nx1: %e\n",debug-pdebug);
					pdebug = debug;
#endif
					/** mxn internal loops (m+n>4) **/
					
					nb_free_nt = minimum(2*__max_size_bulge__,len - __hairpin_min_length__ - 8);
					
					/** lookup closing pair first because it's more expensive **/
					nb_max_mut_close_bp=minimum(4-cst_tape[ii]-cst_tape[ii+1]-cst_tape[jj-1]-cst_tape[jj],kk);
					
					for (nb_mut_close_bp=0;nb_mut_close_bp<=nb_max_mut_close_bp;nb_mut_close_bp++) {
						
						list_config_close = stack_tab[nb_mut_close_bp][ii][jj][1];
						ii_list_close = 0;
						
						while ((config_close=list_config_close[ii_list_close])) {
							
							xx=config_close>>6;
							aa=(config_close>>4)&3;
							bb=(config_close>>2)&3;
							yy=config_close&3;
							
							for (nn=1;nn<=nb_free_nt;nn++) { /* number of undefined nucleotides in internal loop */
								
								min_rr = maximum(0,nn-__max_size_bulge__);
								max_rr = minimum(nn,__max_size_bulge__);
								
								for (rr=min_rr;rr<=max_rr;rr++) {
									ll=nn-rr;
									
									if (!is_internal_loop(ii,ii+3+ll,jj-3-rr,jj)) continue; /* lookup structure constraints */
									
									nt2mut=nn;
									if (ll>0) { nt2mut -= cst_segment[ii+2][ll-1]; }
									if (rr>0) { nt2mut -= cst_segment[jj-1-rr][rr-1]; }
									nb_max_mut_unk_loop = minimum(nt2mut,kk-nb_mut_close_bp);
									
									for (nb_mut_unk_loop=0;nb_mut_unk_loop<=nb_max_mut_unk_loop;nb_mut_unk_loop++) { /* number of mutations in unknown part of loops */
										
										nb_max_mut_open_bp=minimum(4-cst_tape[ii+2+ll]-cst_tape[ii+3+ll]-cst_tape[jj-2-rr]-cst_tape[jj-3-rr],kk-nb_mut_unk_loop-nb_mut_close_bp);
										
										for (nb_mut_open_bp=0;nb_mut_open_bp<=nb_max_mut_open_bp;nb_mut_open_bp++) {
											list_config_open = stack_tab[nb_mut_open_bp][ii+2+ll][jj-2-rr][2];
											ii_list_open = 0;
											while ((config_open=list_config_open[ii_list_open])) {
												
												cc=(config_open>>6)&3;
												uu=(config_open>>4)&3;
												vv=(config_open>>2)&3;
												dd=config_open&3;
												
												newK = kk - nb_mut_close_bp - nb_mut_open_bp - nb_mut_unk_loop + kronecker(ii+3+ll,uu) + kronecker(jj-3-rr,vv);
												if (newK>=0) {
													
#if 0
													if (kk==0) {
														printf("kk=%d,newK=%d,close=%d,open=%d,unk=%d,ll=%d,rr=%d at %d,%d to %d,%d\n",
															   kk,newK,nb_mut_close_bp,nb_mut_open_bp,nb_mut_unk_loop,ll,rr,ii,jj,ii+3+ll,jj-3-rr);
													}
#endif
													
													/* cmpt_mutants = genereMutant(nt2mut,nb_mut_unk_loop); */
													cmpt_mutants = sequence_bias(ii+2,ii+1+ll,nb_mut_unk_loop) + sequence_bias(jj-1-rr,jj-2,nb_mut_unk_loop)
													+ sequence_bias(ii+2,jj-2,nb_mut_unk_loop)
													- sequence_bias(ii+2,jj-2-rr,nb_mut_unk_loop) - sequence_bias(ii+2+ll,jj-2,nb_mut_unk_loop)
													+ sequence_bias(ii+2+ll,jj-2-rr,nb_mut_unk_loop);
													
													if (cmpt_mutants>0) {
														Zs[kk][ii][jj][xx][yy].pf += boltzmannInternal_generic(xx,aa,cc,uu,vv,dd,bb,yy,ll+2,rr+2) *
														Zs[newK][ii+3+ll][jj-3-rr][uu][vv].pf * cmpt_mutants * nucleotide_bias(ii,xx) * nucleotide_bias(ii+1,aa) *
														nucleotide_bias(ii+2+ll,cc) * nucleotide_bias(jj-2-rr,dd) * nucleotide_bias(jj-1,bb) * nucleotide_bias(jj,yy);
														Zs[kk][ii][jj][xx][yy].mfe = minimum_double(Zs[kk][ii][jj][xx][yy].mfe,
																									EInternal_generic(xx,aa,cc,uu,vv,dd,bb,yy,ll+2,rr+2) + Zs[newK][ii+3+ll][jj-3-rr][uu][vv].mfe);
#ifdef DEBUG_SAMPLING
														if ((ii==LEFT_CHECK)&&(jj==RIGHT_CHECK)) {
															debug += boltzmannInternal_generic(xx,aa,cc,uu,vv,dd,bb,yy,ll+2,rr+2)
															* Zs[newK][ii+3+ll][jj-3-rr][uu][vv].pf * cmpt_mutants
															* nucleotide_bias(ii,xx) * nucleotide_bias(ii+1,aa)
															* nucleotide_bias(ii+2+ll,cc) * nucleotide_bias(jj-2-rr,dd)
															* nucleotide_bias(jj-1,bb) * nucleotide_bias(jj,yy);
														}
#endif
													}
												}
												ii_list_open++;
											}
										}
									}
								}
							}
							ii_list_close++;
						}
					}
					
				}
					
#ifdef DEBUG_SAMPLING
				if ((ii==LEFT_CHECK)&&(jj==RIGHT_CHECK)) printf("generic loop: %e\n",debug-pdebug);
				pdebug = debug;
#endif
				
				for (xx=0;xx<4;xx++) {
					
					if (kronecker(ii,xx) && cst_tape[ii]) continue;
					
					for (yy=0;yy<4;yy++) {
						
						if (kronecker(jj,yy) && cst_tape[jj]) continue;
						
						/* lookup bulge & multi-loop WARNING: THIS SECTION HAS NOT BEEN OPTIMIZED */
						
						/* bulges */
						
						if (validBasePair(xx,yy) && is_basepair(ii,jj)) { /* special case when xx and yy base pair */
							
							nb_mut_inside_boundary = kk - kronecker(ii,xx) - kronecker(jj,yy);
							
							if (nb_mut_inside_boundary>=0) {
								
								if (len > __hairpin_min_length__ +4) { /* min length required for a stack */
									
									for (uu=0;uu<4;uu++) {
										
										for (vv=0;vv<4;vv++) {
											
											if (validBasePair(uu,vv)) {
												
												nb_free_nt = minimum(__max_size_bulge__,len - __hairpin_min_length__ - 4);
												
												/* bulge */
												
												for (size_bulge=1;size_bulge<=nb_free_nt;size_bulge++) { /* size of the bulge */
													
													/* we need to distinguish the bulge because the asymetry of constraints */
													
													/* bulge on the left */
													
													newii = ii+1+size_bulge;
													newjj = jj-1;
													
													if (is_internal_loop(ii,newii,newjj,jj)) { /* lookup structure constraints */
														
														if (!(kronecker(newii,uu) && cst_tape[newii])||
															!(kronecker(newjj,vv) && cst_tape[newjj])) {
															
															nt2mut=size_bulge-cst_segment[ii+1][size_bulge-1];
															nb_max_mut_in_bulge = minimum(nt2mut,nb_mut_inside_boundary);
															
															for (nb_mut_in_bulge=0;nb_mut_in_bulge<=nb_max_mut_in_bulge;nb_mut_in_bulge++) { /* number of mutations in bulge */
																newK = nb_mut_inside_boundary  - nb_mut_in_bulge;
																if (newK>=0) {
																	if (Zs[newK][newii][newjj][uu][vv].pf) {
																		
																		/* cmpt_mutants = genereMutant(nt2mut,nb_mut_in_bulge); */
																		cmpt_mutants = sequence_bias(ii+1,ii+size_bulge,nb_mut_in_bulge);
																		
																		if (cmpt_mutants>0) {
																			Zs[kk][ii][jj][xx][yy].pf += boltzmannBulge(xx,uu,vv,yy,size_bulge) *
																			Zs[newK][newii][newjj][uu][vv].pf * cmpt_mutants * nucleotide_bias(ii,xx) * nucleotide_bias(jj,yy);
																			Zs[kk][ii][jj][xx][yy].mfe = minimum_double(Zs[kk][ii][jj][xx][yy].mfe,
																														EBulge(xx,uu,vv,yy,size_bulge) + Zs[newK][newii][newjj][uu][vv].mfe);
#ifdef DEBUG_SAMPLING
																			if ((ii==LEFT_CHECK)&&(jj==RIGHT_CHECK)) debug += boltzmannBulge(xx,uu,vv,yy,size_bulge) *
																				Zs[newK][newii][newjj][uu][vv].pf * cmpt_mutants;
#endif
																		}
																	}
																}
															}
														}
													}
													
													/* bulge on the right */
													
													newii = ii+1;
													newjj = jj-1-size_bulge;
													
													if (is_internal_loop(ii,newii,newjj,jj)) { /* lookup structure constraints */
														
														if (!(kronecker(newii,uu) && cst_tape[newii])||
															!(kronecker(newjj,vv) && cst_tape[newjj])) {
															
															nt2mut=size_bulge-cst_segment[jj-size_bulge][size_bulge-1];
															nb_max_mut_in_bulge = minimum(nt2mut,nb_mut_inside_boundary);
															
															for (nb_mut_in_bulge=0;nb_mut_in_bulge<=nb_max_mut_in_bulge;nb_mut_in_bulge++) { /* number of mutations in bulge */
																newK = nb_mut_inside_boundary  - nb_mut_in_bulge;
																if (newK>=0) {
																	if (Zs[newK][newii][newjj][uu][vv].pf) {
																		
																		/* cmpt_mutants = genereMutant(nt2mut,nb_mut_in_bulge); */
																		cmpt_mutants = sequence_bias(jj-size_bulge,jj-1,nb_mut_in_bulge);
																		
																		if (cmpt_mutants>0) {
																			Zs[kk][ii][jj][xx][yy].pf += boltzmannBulge(xx,uu,vv,yy,size_bulge) *
																			Zs[newK][newii][newjj][uu][vv].pf * cmpt_mutants * nucleotide_bias(ii,xx) * nucleotide_bias(jj,yy);
																			Zs[kk][ii][jj][xx][yy].mfe = minimum_double(Zs[kk][ii][jj][xx][yy].mfe,
																														EBulge(xx,uu,vv,yy,size_bulge) + Zs[newK][newii][newjj][uu][vv].mfe);
#ifdef DEBUG_SAMPLING
																			if ((ii==LEFT_CHECK)&&(jj==RIGHT_CHECK)) debug += boltzmannBulge(xx,uu,vv,yy,size_bulge) *
																				Zs[newK][newii][newjj][uu][vv].pf * cmpt_mutants;
#endif
																		}
																	}
																}
															}
														}	      
													}
												}
											}
											
											/* close multi-loop*/
											
											if ((ii<jj)&&(is_basepair(ii,jj))) { /* lookup structure constraints */
												
												Zs[kk][ii][jj][xx][yy].pf += boltzmannEndMultiLoop() * Zm[nb_mut_inside_boundary][ii+1][jj-1][uu][vv].pf * nucleotide_bias(ii,xx) * nucleotide_bias(jj,yy);
												Zs[kk][ii][jj][xx][yy].mfe = minimum_double(Zs[kk][ii][jj][xx][yy].mfe,EEndMultiLoop() + Zm[nb_mut_inside_boundary][ii+1][jj-1][uu][vv].mfe);
												
#if 0
												if (Zm[nb_mut_inside_boundary][ii+1][jj-1][uu][vv].pf)
													printf("%d,%d: %e (%e * %e)\n",ii,jj,Zs[kk][ii][jj][xx][yy].pf,
														   boltzmannEndMultiLoop(),Zm[nb_mut_inside_boundary][ii+1][jj-1][uu][vv].pf);
#endif
#ifdef DEBUG_SAMPLING
												if ((ii==LEFT_CHECK)&&(jj==RIGHT_CHECK)) debug += boltzmannEndMultiLoop() * Zm[nb_mut_inside_boundary][ii+1][jj-1][uu][vv].pf;
#endif
											}
										}
									}
								}
							}
							
							/* initialize related arrays */
							
							Zms[kk][ii][jj][xx][yy].pf  = Zes[kk][ii][jj][xx][yy].pf  = Zs[kk][ii][jj][xx][yy].pf * boltzmann_penalty_close_helix(xx,yy);
							Zms[kk][ii][jj][xx][yy].mfe = Zes[kk][ii][jj][xx][yy].mfe = Zs[kk][ii][jj][xx][yy].mfe + penalty_close_helix(xx,yy);
							
						}
						
						/* try intermediate pairings in order to build multi-loop */
						
						nb_mut_inside_boundary = kk - kronecker(ii,xx);
						
						for (rr=ii+1;rr<=jj-__hairpin_min_length__-1;rr++) {
							
							for (uu=0;uu<4;uu++) {
								
								/* initiate array and extend on left-hand side */
								
								if (single_strand_length(rr-1)>=(rr-ii-1)) { /* lookup structure constraints */
									
									if ((!((cst_tape[ii])&&(kronecker(ii,xx))))&&(validBasePair(uu,yy))) {
										
										if (rr-ii>1) {
											
											nt2mut = rr-ii-1 - cst_segment[ii+1][rr-ii-2];
											newK = minimum(nt2mut,nb_mut_inside_boundary);
											
											for (ee=0;ee<=newK;ee++) {
												
												/* cmpt_mutants = genereMutant(nt2mut,ee); */
												cmpt_mutants = sequence_bias(ii+1,rr-1,ee) * nucleotide_bias(ii,xx);
												
												if (cmpt_mutants > 0) {
													Zes[kk][ii][jj][xx][yy].pf += Zs[nb_mut_inside_boundary-ee][rr][jj][uu][yy].pf * boltzmann_penalty_close_helix(uu,yy) * cmpt_mutants;
													Zes[kk][ii][jj][xx][yy].mfe = minimum_double(Zes[kk][ii][jj][xx][yy].mfe,
                                                                                                 Zs[nb_mut_inside_boundary-ee][rr][jj][uu][yy].mfe + penalty_close_helix(uu,yy)); 
													Zms[kk][ii][jj][xx][yy].pf += Zs[nb_mut_inside_boundary-ee][rr][jj][uu][yy].pf * boltzmann_penalty_close_helix(uu,yy) * cmpt_mutants * boltzmannExtendMultiLoop(rr-ii);
													Zms[kk][ii][jj][xx][yy].mfe = minimum_double(Zms[kk][ii][jj][xx][yy].mfe,
																								 Zs[nb_mut_inside_boundary-ee][rr][jj][uu][yy].mfe + penalty_close_helix(uu,yy) + EExtendMultiLoop(rr-ii)); 
												}
											}
										}
										else {
											newK = kk - kronecker(ii,xx);
											if (newK>=0) {
												Zes[kk][ii][jj][xx][yy].pf += Zs[newK][rr][jj][uu][yy].pf * boltzmann_penalty_close_helix(uu,yy) * nucleotide_bias(ii,xx);
												Zes[kk][ii][jj][xx][yy].mfe = minimum_double(Zes[kk][ii][jj][xx][yy].mfe,Zs[newK][rr][jj][uu][yy].mfe + penalty_close_helix(uu,yy));
												Zms[kk][ii][jj][xx][yy].pf += Zs[newK][rr][jj][uu][yy].pf * boltzmann_penalty_close_helix(uu,yy) * nucleotide_bias(ii,xx) * boltzmannExtendMultiLoop(1);
												Zms[kk][ii][jj][xx][yy].mfe = minimum_double(Zms[kk][ii][jj][xx][yy].mfe,Zs[newK][rr][jj][uu][yy].mfe + penalty_close_helix(uu,yy) + EExtendMultiLoop(1));
											}
										}
									}
								}
								
								if (rr-ii>__hairpin_min_length__+1) {
									for (vv=0;vv<4;vv++) {
										if (validBasePair(vv,yy)) {
											for (ee=0;ee<=kk;ee++) {
												
												/* more than 2 helices */
												Zm[kk][ii][jj][xx][yy].pf += Zm[ee][ii][rr-1][xx][uu].pf * Zs[kk-ee][rr][jj][vv][yy].pf * boltzmann_penalty_close_helix(vv,yy) * boltzmannAddStemInMultiLoop(1);
												Zm[kk][ii][jj][xx][yy].mfe = minimum_double(Zm[kk][ii][jj][xx][yy].mfe,
																							Zm[ee][ii][rr-1][xx][uu].mfe + Zs[kk-ee][rr][jj][vv][yy].mfe + penalty_close_helix(vv,yy) + EAddStemInMultiLoop(1));
												Ze[kk][ii][jj][xx][yy].pf += Ze[ee][ii][rr-1][xx][uu].pf * Zs[kk-ee][rr][jj][vv][yy].pf * boltzmann_penalty_close_helix(vv,yy);
												Ze[kk][ii][jj][xx][yy].mfe = minimum_double(Ze[kk][ii][jj][xx][yy].mfe,
																							Ze[ee][ii][rr-1][xx][uu].mfe + Zs[kk-ee][rr][jj][vv][yy].mfe + penalty_close_helix(vv,yy));
												/* exactly 2 helices */
												Zm[kk][ii][jj][xx][yy].pf += Zms[ee][ii][rr-1][xx][uu].pf * Zs[kk-ee][rr][jj][vv][yy].pf * boltzmann_penalty_close_helix(vv,yy) * boltzmannAddStemInMultiLoop(2);
												Zm[kk][ii][jj][xx][yy].mfe = minimum_double(Zm[kk][ii][jj][xx][yy].mfe,
																							Zms[ee][ii][rr-1][xx][uu].mfe + Zs[kk-ee][rr][jj][vv][yy].mfe + penalty_close_helix(vv,yy) + EAddStemInMultiLoop(2));
												Ze[kk][ii][jj][xx][yy].pf += Zes[ee][ii][rr-1][xx][uu].pf * Zs[kk-ee][rr][jj][vv][yy].pf * boltzmann_penalty_close_helix(vv,yy);
												Ze[kk][ii][jj][xx][yy].mfe = minimum_double(Ze[kk][ii][jj][xx][yy].mfe,
																							Zes[ee][ii][rr-1][xx][uu].mfe + Zs[kk-ee][rr][jj][vv][yy].mfe + penalty_close_helix(vv,yy));
											}
										}
									}
								}
							}
						}
						
						/* extend all tables with an unpaired nucleotide on right side */
						
						if (single_strand_length(jj)>0) { /* lookup structure constraints */
							
							if (len>__hairpin_min_length__+2) { /* minimum length required */
								if (!((cst_tape[jj])&&(kronecker(jj,yy)))) {
									for (uu=0;uu<4;uu++) {
										/* right */
										newK = kk - kronecker(jj,yy);
										if (newK>=0) {
											Zes[kk][ii][jj][xx][yy].pf += Zes[newK][ii][jj-1][xx][uu].pf * nucleotide_bias(jj,yy);
											Zes[kk][ii][jj][xx][yy].mfe = minimum_double(Zes[kk][ii][jj][xx][yy].mfe,Zes[newK][ii][jj-1][xx][uu].mfe);
											Ze[kk][ii][jj][xx][yy].pf += Ze[newK][ii][jj-1][xx][uu].pf * nucleotide_bias(jj,yy);
											Ze[kk][ii][jj][xx][yy].mfe = minimum_double(Ze[kk][ii][jj][xx][yy].mfe,Ze[newK][ii][jj-1][xx][uu].mfe);
											Zms[kk][ii][jj][xx][yy].pf += boltzmannExtendMultiLoop(1) * Zms[newK][ii][jj-1][xx][uu].pf * nucleotide_bias(jj,yy);
											Zms[kk][ii][jj][xx][yy].mfe = minimum_double(Zms[kk][ii][jj][xx][yy].mfe,EExtendMultiLoop(1) + Zms[newK][ii][jj-1][xx][uu].mfe);
											Zm[kk][ii][jj][xx][yy].pf += boltzmannExtendMultiLoop(1) * Zm[newK][ii][jj-1][xx][uu].pf * nucleotide_bias(jj,yy);
											Zm[kk][ii][jj][xx][yy].mfe = minimum_double(Zm[kk][ii][jj][xx][yy].mfe,EExtendMultiLoop(1) + Zm[newK][ii][jj-1][xx][uu].mfe);
										}
									}
								}
							}  
						}
					}
				}
				
#ifdef DEBUG_SAMPLING
				if ((ii==LEFT_CHECK)&&(jj==RIGHT_CHECK)) printf("bulge & multi: %e\n",debug-pdebug);
				pdebug = debug;
#endif
			}
		}
		
		
		if (verbose) {
			current = clock();
			cpu_time_used = ((double) (current - start)) / CLOCKS_PER_SEC;	
			printf(">%d\t%g\n",kk,cpu_time_used);
		}
		
	}
	
	fprintf(stderr,"\r                   \r");fflush(stderr);
#ifdef INFO
	printf("\n");
#endif
	
}


/*********************************************************************/
/*                                                                   */
/* Main                                                              */
/*                                                                   */
/*********************************************************************/

int main(int argc, char **argv) {
	
	int param_temperature=37; /* define as integer as in mfold */
	char *samplefile=NULL,*param_temperature_as_string = strdup("37");
	double formal_temperature=0.0;  /* used in Boltzmann computation */
	const char *input_filename=NULL, *mutation_matrix_filename=NULL;
	const char *input_string=NULL;
	static int verbose_flag=0, mfold3_flag=1,input_string_flag=0,input_file_flag=0,sample_output_stats_flag=0,hybrid_flag=0,cst_flag=0;
	char *input_param_directory=NULL, *library_path=strdup("lib");
	char *partition_function_output_filename=NULL, *mfe_output_filename=NULL, *sample_output_filename=NULL;
    FILE *mfe_file_flux=NULL, *pf_file_flux=NULL, *sample_file_flux=NULL;
	int *compatible_neighbors,compatible_neighbors_flag=0;
	
	static struct option long_options[] =
	{
		{"file", required_argument,	0, 'f'},
		{"string", required_argument,	0, 's'},
		{"constraints",no_argument,	0, 'C'},
		{"turner99", no_argument,	&mfold3_flag, 1},
		{"mutations",required_argument,	0, 'm'},
		{"sample-number",required_argument,	0, 'n'},
		{"sample-file",required_argument,	0, 'N'},
		{"sample-stats",no_argument,	&sample_output_stats_flag, 1},
		{"Z-output",required_argument,	0, '0'},
		{"mfe-output",required_argument,	0, '1'},
		{"sample-output",required_argument,	0, '2'},
		{"temperature",required_argument,	0, 't'},
		{"absolute-temperature",required_argument,	0, 'T'},
		{"help",   no_argument, 0, 'h'},
		{"version", no_argument, 0, 'V'},
		{"warning", no_argument, 0, 'w'},
		{"verbose", no_argument, 0, 'v'},
		{"energy-model",required_argument,	0, 'e'},
		{"library",required_argument,	0, 'l'},
		{"hybrid",no_argument, 0, 'H'},
		{"job-id", required_argument,	0, 'J'},
		{"bias", required_argument,	0, 'M'},
		{0, 0, 0, 0}
	};
		
	/* parse line command */
	
	while (1) {
		
		/* getopt_long stores the option index here. */
		int c, option_index = 0;
		
		c = getopt_long (argc, argv, "f:s:l:CHm:n:N:M:e:t:T:hJ:vVw0:1:2:", long_options, &option_index);
		
		/* Detect the end of the options. */
		if (c == -1)
			break;
		
		switch (c) {
			case 'f':
				input_filename=strdup(optarg);
				input_file_flag = 1;
				break;
			case 's':
				input_string=strdup(optarg);
				input_string_flag = 1;
				break;
			case 'C':
				cst_flag=1;
				break;
			case 'm':
				max_mutations = atoi(optarg)+1;
				if (max_mutations<1) {
					fprintf(stderr,"WARNING: Corrupted number of mutations. Ignored. %s\n", getJobID());
					max_mutations=0;
				}
				break;
			case 'M':
				mutation_matrix_filename = strdup(optarg);
				break;
			case 'n':
				number_of_sample = atoi(optarg);
				break;
			case 'N':
				samplefile = strdup(optarg);
				break;
			case '0':
				partition_function_output_filename = strdup(optarg);
				break;
			case '1':
				mfe_output_filename = strdup(optarg);
				break;
			case '2':
				sample_output_filename = strdup(optarg);
				break;
			case 'h':
				usage(argv[0]);
			case 'v':
				verbose_flag = 1;
				break;
			case 'V':
				version();
			case 't':
				param_temperature = atoi(optarg);
				param_temperature_as_string = strdup(optarg);
				mfold3_flag=0; /* remove mfold3 flag */
				break;
			case 'T':
				formal_temperature = atof(optarg);
				break;
			case 'e':
				input_param_directory = strdup(optarg);
				break;
			case 'l':
				library_path = strdup(optarg);
				break;
			case 'x':
				cutpoint = atoi(optarg);
				break;
			case 'H':
				hybrid_flag = 1;
				break;
			case 'J':
				setJobID(optarg);
				break;
		}
	}
	
	/* check if argv[optind] is valid */
	
	if ((warning_flag)&&(optind != (argc))) {
		int iargv;
		fprintf(stderr,"WARNING: invalid syntax. Following argument(s) are ignored: ");
		for (iargv=optind;iargv<argc;iargv++) {
			fprintf(stderr,"%s ",argv[iargv]);
		}
		fprintf(stderr,"\n");
	}
	
	/* check commands */
	
	if (mfold3_flag) {
		if (param_temperature != 37) {
			if (warning_flag) {
				fprintf(stderr,"WARNING: Temperature %dC not supported by mfold-3 energy model. Value set to 37C.%s\n",param_temperature,getJobID());
			}
		}
		param_temperature = 37;
		param_temperature_as_string = strdup("dat");
	}				
	
	if ((param_temperature < 0)||(param_temperature > 100)) {
		if (param_temperature < 0) {
			param_temperature = 0;
		}
		if (param_temperature > 100) {
			param_temperature = 100;
		}
		if (warning_flag) {
			fprintf(stderr,"WARNING: Temperature out of range. Value must range between 0 and 100. Updated to %dC.%s\n",param_temperature,getJobID());
		}
	}
	
	if (formal_temperature < 0) {
		formal_temperature = 0.0;
		if (warning_flag) {
			fprintf(stderr,"WARNING: Negative formal temperature is not supported. Value updated to %fK.%s\n",formal_temperature,getJobID());
		}
	}
	
	/* print info on system precision */
	
	if (verbose_flag) {
		printf("> system features : epsilon=%e, max. value=%e, precision=%d digit\n",
			   DBL_EPSILON, DBL_MAX, DBL_DIG);
	}
	
	/* read the data and store then in tables */
	
	{
		char basename[200], filename[200];
		
		if (verbose_flag) {
			printf("> read parameter files defined for %s Celsius degrees\n",param_temperature_as_string);
		}
		
		/* build path to the directory where files are stored */
		
		basename[0] = '\0';
		if (input_param_directory) {
			if ((library_path)&&(warning_flag)) {
				fprintf(stderr,"WARNING: Library path directory not used.%s\n",getJobID());
			}
			strcat(&basename[0],input_param_directory);
			param_temperature_as_string = strdup("dat");
		}
		else {
			if (library_path) {
				strcat(&basename[0],library_path);
				strcat(&basename[0],"/");
				strcat(&basename[0],param_temperature_as_string);
			}
			else {
				strcat(&basename[0],"./");
				strcat(&basename[0],param_temperature_as_string);
			}
		}
		
		strcat(&basename[0],"/");
		
		/* read each parameter files */
		
		/* asint1x2 */
		memcpy(&filename[0],&basename[0],200*sizeof(char));
		strcat(&filename[0],"asint1x2.lop.");
		strcat(&filename[0],param_temperature_as_string);
		readAsint1x2(filename);
		
		/* dangle */
		memcpy(&filename[0],&basename[0],200*sizeof(char));
		strcat(&filename[0],"dangle.lop.");
		strcat(&filename[0],param_temperature_as_string);
		readDangle(filename);
		
		/* loop */
		memcpy(&filename[0],&basename[0],200*sizeof(char));
		strcat(&filename[0],"loop.lop.");
		strcat(&filename[0],param_temperature_as_string);
		readLoop(filename);
		
		/* miscloop */
		memcpy(&filename[0],&basename[0],200*sizeof(char));
		strcat(&filename[0],"miscloop.lop.");
		strcat(&filename[0],param_temperature_as_string);
		readMiscLoop(filename);
		
		/* stack */
		memcpy(&filename[0],&basename[0],200*sizeof(char));
		strcat(&filename[0],"stack.lop.");
		strcat(&filename[0],param_temperature_as_string);
		readStack(filename);
		
		/* sint2 */
		memcpy(&filename[0],&basename[0],200*sizeof(char));
		strcat(&filename[0],"sint2.lop.");
		strcat(&filename[0],param_temperature_as_string);
		readSint2(filename);
		
		/* sint4 */
		memcpy(&filename[0],&basename[0],200*sizeof(char));
		strcat(&filename[0],"sint4.lop.");
		strcat(&filename[0],param_temperature_as_string);
		readSint4(filename);
		
		/* tloop */
		memcpy(&filename[0],&basename[0],200*sizeof(char));
		strcat(&filename[0],"tloop.lop.");
		strcat(&filename[0],param_temperature_as_string);
		readTetraLoop(filename);
		
		/* triloop */
		memcpy(&filename[0],&basename[0],200*sizeof(char));
		strcat(&filename[0],"triloop.lop.");
		strcat(&filename[0],param_temperature_as_string);
		readTriLoop(filename);
		
		/* tsatcki */
		memcpy(&filename[0],&basename[0],200*sizeof(char));
		strcat(&filename[0],"tstacki.lop.");
		strcat(&filename[0],param_temperature_as_string);
		readTstacki(filename);
		
		/* tstackh */
		memcpy(&filename[0],&basename[0],200*sizeof(char));
		strcat(&filename[0],"tstackh.lop.");
		strcat(&filename[0],param_temperature_as_string);
		readTstackh(filename);
		
	}
	
	/* update formal temperature */
	
	{ 
		if (formal_temperature) {
			if (__temperature != formal_temperature) {
				__temperature = formal_temperature;
			}
		}
		else {
			__temperature = (double)(param_temperature) + 273.0;
		}
		
		
		if (verbose_flag) {
			printf("> formal temperature set to %.2f Kelvin degrees\n",__temperature);
		}
		
		/* update Boltzmann parameters */
		
		__RT__ = (double)(__temperature * GAS_CONSTANT)/1000.0; /* RT in kCal */
		
	}
	
    /* open output file descriptors */
    
    if (mfe_output_filename) {
        mfe_file_flux = fopen(mfe_output_filename,"w");
    } else {
        mfe_file_flux = stdout;
    }
    
    if (partition_function_output_filename) {
        pf_file_flux = fopen(partition_function_output_filename,"w");
    } else {
        pf_file_flux = stdout;
    }
    
    if (sample_output_filename) {
        sample_file_flux = fopen(sample_output_filename,"w");
    } else {
        sample_file_flux = stdout;
    }
    
	/* init dynamic table and check it */
	
	if (verbose_flag) {
		starting_time = time(0);
		printf("> allocating memory for dynamic tables.\n");
	}
	
	/* init input tape */
	
	if (!(input_file_flag || input_string_flag)) {
		fprintf(stderr,"Error: Input not provided.%s\n",getJobID());
		exit(EXIT_FAILURE);
	}
	else if (input_file_flag && input_string_flag) {
		fprintf(stderr,"Error: Dual input are not allowed.%s\n",getJobID());
		exit(EXIT_FAILURE);
	}
	
	if (input_file_flag) {
		init_input_tape_from_file(input_filename,cst_flag);
	}
	else {
		init_input_tape_from_string(input_string);
	}
	
	if (rna_len<1) {
		fprintf(stderr,"Error: Invalid input sequence.%s\n",getJobID());
		exit(EXIT_FAILURE);
	}
	
	/* init mutation weights */
	
	if (mutation_matrix_filename) {
		init_mutation_matrix_from_file(mutation_matrix_filename);
	}

	/* set maximal number of mutations to maximum value if not set */
	if (!max_mutations) {
		if (cst_flag) { max_mutations = cmpt_open_mutations() + 1; }
		else { max_mutations = rna_len + 1; }
	}
	else {
		if (cst_flag) {
			int max_mutations_open = cmpt_open_mutations() + 1;
			if (max_mutations>max_mutations_open) {
				fprintf(stderr,"WARNING: Maximal number of mutations cannot exceed the number of mutation sites. Value set to %d.%s\n",max_mutations_open-1,getJobID());
				max_mutations = max_mutations_open;
			}
		}
		else if (max_mutations>rna_len+1) {
			fprintf(stderr,"WARNING: Maximal number of mutations cannot exceed sequence length. Ignored.%s\n",getJobID());
			max_mutations = rna_len + 1;
		}
	}
	
	
	/* check cutpoint */
	
	if (hybrid_flag) {
		if (!cutpoint) {
			if (warning_flag) {
				fprintf(stderr,"WARNING: Missing cutpoint betwween sequences. Hybrid mode ignored.%s\n",getJobID());
			}
			hybrid_flag=0;
		}
	}
	else {
		if (cutpoint) {
			if (verbose_flag) {
				/* printf("> cutpoint used without hybridization mode."); */
				fprintf(stderr,"WARNING: cutpoint used out of hybrid mode. Ignored.%s\n",getJobID());
			}
			cutpoint=0;
			include_intermolecular_interactions=0;
		}
	}
	
	/* init arrays */
	if (hybrid_flag) {
		init_Mutant_arrays_hybrid();		
	}
	else {
		init_Mutant_arrays();
	}
	
	/* init table for tri and tetraloop counting */
	init_triloop_cmpt_table();
	init_tetraloop_cmpt_table();
	
	if (verbose_flag) {
		printf("> allocated memory in dynamic tables : %lld bytes\n",total_memory);
		printf("> time used for memory allocation : %ld secondes\n",time(0)-starting_time);
		starting_time = time(0);
	}

	/* precompute internal loop tables */
	
	{
		int nb_mut,max_nb_mut;
		
		max_nb_mut=minimum(9,max_mutations);
		loop_tab = (unsigned int *****)xmalloc(max_nb_mut*sizeof(unsigned int ****));
		for (nb_mut=0;nb_mut<max_nb_mut;nb_mut++) {
			loop_tab[nb_mut]=enum_all_special_internal_loops(nb_mut);
		}
		
		max_nb_mut=minimum(5,max_mutations);
		stack_tab = (unsigned int *****)xmalloc(max_nb_mut*sizeof(unsigned int ****));
		for (nb_mut=0;nb_mut<max_nb_mut;nb_mut++) {
			stack_tab[nb_mut]=enum_all_stacks(nb_mut);
		}
	}
	
	if (verbose_flag) {
		printf("> time used for internal loop tables : %ld secondes\n",time(0)-starting_time);
		printf("> sequence length : %d nt.\n",rna_len);
		starting_algo_time = time(0);
	}
	
	/* compute Boltzmann partition function */
	
	/* init dynamic tables with hairpins */
	
	if (hybrid_flag) {
		precomputeHairpin_hybrid(cutpoint,include_intermolecular_interactions);
	}
	else {
		precomputeHairpin();
	}
	
	/* Main algorithm */
	
	if (hybrid_flag) {
		RNAmutants_algorithm_hybrid(verbose_flag);
	}
	else {
		RNAmutants_algorithm(verbose_flag);
	}
	
	/* completed display results */
	
	if (verbose_flag) {
		printf("> time used for computation : %ld secondes\n",time(0)-starting_algo_time);
		printf("> total time used : %ld secondes\n",time(0)-starting_time);
	}
	
	
	{
		int uu,vv,kk,iLast=rna_len-1;
		double boltzmann_distribution_simple[max_mutations], boltzmann_distribution_multi[max_mutations];
		double mfe_v[max_mutations];
		
		compatible_neighbors = (int *)xmalloc((rna_len+1)*sizeof(int));
		memset(compatible_neighbors,0,(rna_len+1)*sizeof(int));
		
		if (!ss_constraint_is_non_empty()) { /* add empty structure iff constraints allow it */
			for(kk=0;kk<max_mutations;kk++) {
				/* boltzmann_distribution_simple[kk] = combinat(rna_len,kk) * pow(3,kk); */
				boltzmann_distribution_simple[kk] = sequence_bias(0,iLast,kk);
				boltzmann_distribution_multi[kk] = 0; 
				mfe_v[kk] = 0.0;
			}
		}
		else {	
			for(kk=0;kk<max_mutations;kk++) {
				boltzmann_distribution_simple[kk] = 0; 
				boltzmann_distribution_multi[kk] = 0; 
				mfe_v[kk] = __infinity__;
			}
		}
		
		for(kk=0;kk<max_mutations;kk++) {
			for (uu=0;uu<4;uu++) {
				for (vv=0;vv<4;vv++) {
					if (Zes[kk][0][iLast][uu][vv].pf>DBL_EPSILON) {
						boltzmann_distribution_simple[kk] += Zes[kk][0][iLast][uu][vv].pf;
					}
					mfe_v[kk] = minimum_double(mfe_v[kk],Zes[kk][0][iLast][uu][vv].mfe);
					if (Ze[kk][0][iLast][uu][vv].pf>DBL_EPSILON) {
						boltzmann_distribution_multi[kk] += Ze[kk][0][iLast][uu][vv].pf;
					}
					mfe_v[kk] = minimum_double(mfe_v[kk],Ze[kk][0][iLast][uu][vv].mfe);
				} 
			}
		}
		
		/* check precision before printing */
		
		if (DBL_DIG != 15) {
			fprintf(stderr,"WARNING: printf format might not be suited to double precision.%s\n",getJobID());
		}
		
		/* print results */
		
        /* Partition function */
        
#ifdef DETAILED
		for(kk=0;kk<max_mutations;kk++) {
			fprintf(pf_file_flux,"%d %.14e %.14e %.14e\n",kk,boltzmann_distribution_simple[kk]+boltzmann_distribution_multi[kk],
				   boltzmann_distribution_simple[kk],boltzmann_distribution_multi[kk]);
		}
#else
		for(kk=0;kk<max_mutations;kk++) {
			fprintf(pf_file_flux,"%d %.14e\n",kk,boltzmann_distribution_simple[kk]+boltzmann_distribution_multi[kk]);
		}
#endif
		
		/* MFE */
		
		fprintf(mfe_file_flux,">> Superoptimal structure(s):\n");
		for(kk=0;kk<max_mutations;kk++) {
			if ((ss_constraint_is_non_empty())&&(boltzmann_distribution_simple[kk]+boltzmann_distribution_multi[kk]==0)) {
				fprintf(mfe_file_flux,"> %d-superoptimal structure\n",kk);
				fprintf(mfe_file_flux,"None: no secondary structure satisfy the constraints.\n");
			}
			else {
				compatible_neighbors_flag = 1;
				compatible_neighbors[kk]=1;
				startBacktrackKSuperOptimal(mfe_file_flux,kk,mfe_v[kk]);
			}
		}
		
	}

	/* Sampling */
	
	if (number_of_sample) {
		fprintf(sample_file_flux,">> Sampling %d mutants:\n",number_of_sample);
		basicSamplingEngine(sample_file_flux,number_of_sample,sample_output_stats_flag,warning_flag,compatible_neighbors_flag);
	}
	
	if (samplefile) {
		fprintf(sample_file_flux,">> Sampling mutants from file %s:\n",samplefile);
		sampleFromFile(sample_file_flux,samplefile,sample_output_stats_flag,warning_flag,compatible_neighbors);
	}
	
    /** end of RNAmutants -- Cleaning now **/
    
    /* close output file descriptors */
    
    if (mfe_output_filename) {
        fclose(mfe_file_flux);
    }
    
    if (partition_function_output_filename) {
        fclose(pf_file_flux);
    }

    if (sample_output_filename) {
        fclose(sample_file_flux);
    }

#ifdef CLEAN_DYNAMIC_TABLE
	
	{
		int ileft,iright,uu,vv,kk;
		for(kk=0;kk<max_mutations;kk++) {
			for (ileft=0;ileft<rna_len;ileft++) {
				for (iright=ileft+__hairpin_min_length__;iright<rna_len;iright++) {
					for (uu=0;uu<4;uu++) {
						for (vv=0;vv<4;vv++) {
							if (!(Zs[kk][ileft][iright][uu][vv].pf>DBL_EPSILON)) {
								Zs[kk][ileft][iright][uu][vv].pf=0.0;
							}
							if (!(Zms[kk][ileft][iright][uu][vv].pf>DBL_EPSILON)) {
								Zms[kk][ileft][iright][uu][vv].pf=0.0;
							}
							if (!(Zm[kk][ileft][iright][uu][vv].pf>DBL_EPSILON)) {
								Zm[kk][ileft][iright][uu][vv].pf=0.0;
							}
							if (!(Zes[kk][ileft][iright][uu][vv].pf>DBL_EPSILON)) {
								Zes[kk][ileft][iright][uu][vv].pf=0.0;
							}
							if (!(Ze[kk][ileft][iright][uu][vv].pf>DBL_EPSILON)) {
								Ze[kk][ileft][iright][uu][vv].pf=0.0;
							}
#ifdef DANGLE_ARRAY
							if (!(ZdangN[kk][ileft][iright][uu][vv].pf>DBL_EPSILON)) {
								ZdangN[kk][ileft][iright][uu][vv].pf=0.0;
							}
							if (!(ZdangB[kk][ileft][iright][uu][vv].pf>DBL_EPSILON)) {
								ZdangB[kk][ileft][iright][uu][vv].pf=0.0;
							}
							if (!(ZdangL[kk][ileft][iright][uu][vv].pf>DBL_EPSILON)) {
								ZdangL[kk][ileft][iright][uu][vv].pf=0.0;
							}
							if (!(ZdangR[kk][ileft][iright][uu][vv].pf>DBL_EPSILON)) {
								ZdangR[kk][ileft][iright][uu][vv].pf=0.0;
							}
#endif
						} 
					}
				}
			}
		}
	}

#endif
	
	
	return EXIT_SUCCESS;
	
}

