/*
 *  constraints.c
 *  
 *
 *  Created by Jerome Waldispuhl on 27/07/2010.
 *  Copyright 2010 McGill. All rights reserved.
 *
 */

#include<stdio.h>
#include<stdlib.h>
#include <string.h>
#include "util.h"
#include "constraints.h"
#include "energy_params_tables.h"

/* from RNAmutants.h */

extern const int __hairpin_min_length__;
extern const int __hairpin_max_length__;
extern int rna_len;
extern double cutpoint;
#define INDEX_A 0
#define INDEX_C 1
#define INDEX_G 2
#define INDEX_U 3
#define INDEX_STOP -1

/* structure constraints tables */

int *ss_cst;
int *single_strand_cst; /* Unpaired nt before or at this index */
int is_non_empty=0; /* boolean used to determine if the structure can be empty */

/* sequence constraints are declared in RNAmutants.c */
extern int *cst_tape;
extern int **cst_segment;

/* predicats */

int ss_constraint_is_non_empty() {
	return is_non_empty;
}

int is_hairpin(int i, int j) {
	if ((i<j)&&
		((ss_cst[i]==j)||((ss_cst[i]<0)&&(ss_cst[j]<0)))&&
		(single_strand_cst[j-1]>=j-i-1)) {
		return 1;
	}
	else {
		return 0;
	}
}

int is_basepair(int i, int j) {
	if ((ss_cst[i]==j)||((ss_cst[i]<0)&&(ss_cst[j]<0))) {
		return 1;
	}
	return 0;
}

int is_stack(int i, int j) {
	if ((i<j)&&
		((ss_cst[i]==j)||((ss_cst[i]<0)&&(ss_cst[j]<0)))&&
		((ss_cst[i+1]==j-1)||((ss_cst[i+1]<0)&&(ss_cst[j-1]<0)))) {
		return 1;
	}
	else {
		return 0;
	}
}

int is_internal_loop(int i, int m, int n, int j) { /* used for bulge too */
	if (((i<m)&&(m<n)&&(n<j)&&((m-i>1)||(j-n>1)))&&
		((ss_cst[i]==j)||((ss_cst[i]<0)&&(ss_cst[j]<0)))&&
		((ss_cst[m]==n)||((ss_cst[m]<0)&&(ss_cst[n]<0)))&&
		(single_strand_cst[m-1]>=(m-i-1))&&(single_strand_cst[j-1]>=(j-n-1))) {	
		return 1;
	}
	else {
		return 0;
	}
}

int single_strand_length(int ii) {
	return single_strand_cst[ii];
}


int is_triloop(int ii, int jj, char **ss_sample) {

	if (jj-ii==4) {
		int hcode =
			(convert2index(ss_sample[0][ii])<<8) |
			(convert2index(ss_sample[0][ii+1])<<6) |
			(convert2index(ss_sample[0][ii+2])<<4) |
			(convert2index(ss_sample[0][ii+3])<<2) |
			convert2index(ss_sample[0][ii+4]);
		if (refTableTriLoop[hcode]) { return 1; }
	}
	
	return 0;
}

int is_tetraloop(int ii, int jj, char **ss_sample) {

	if (jj-ii==5) {
		int hcode =
		(convert2index(ss_sample[0][ii])<<10) |
		(convert2index(ss_sample[0][ii+1])<<8) |
		(convert2index(ss_sample[0][ii+2])<<6) |
		(convert2index(ss_sample[0][ii+3])<<4) |
		(convert2index(ss_sample[0][ii+4])<<2) |
		convert2index(ss_sample[0][ii+5]);
		if (refTableTetraLoop[hcode]) { return 1; }
	}
	
	
	return 0;
}

int is_ggg_hairpin(int ii, int jj, char **ss_sample) {

	if ((convert2index(ss_sample[0][ii])==INDEX_G)&&
		(convert2index(ss_sample[0][ii+1])==INDEX_G)&&
		(convert2index(ss_sample[0][ii+2])==INDEX_G)&&
		(convert2index(ss_sample[0][jj])==INDEX_U)) {
		return 1;
	}
	return 0;
}


/* fill table functions */

void init_structure_cst_tables(int length_seq) {
	int ii;
	
	ss_cst = xmalloc(length_seq * sizeof(int));
	single_strand_cst = xmalloc(length_seq * sizeof(int));
	
	for (ii=0;ii<length_seq;ii++) {
		ss_cst[ii]=-1;
		single_strand_cst[ii]=ii+1;
	}
	
}

void clean_structure_cst_tables(int length_seq) {
	int ii;
	
	for (ii=0;ii<length_seq;ii++) {
		ss_cst[ii]=-1;
		single_strand_cst[ii]=ii+1;
	}
	
}

void fill_structure_cst_tables(char *buffer) {
	
	int head=0, stack[1000], ii=0, openbp, closebp, prevbp='.', carac;
	
	while (buffer[ii]!='\0') {
		
		if (head<0) {
			fprintf(stderr,"WARNING: Corrupted structure constraints. Unbalanced structure. Constraints ignored.%s\n",getJobID());
			clean_structure_cst_tables(strlen(buffer));
			return;
		}
		
		carac = buffer[ii];
		switch(carac) {
			case '.':
				ss_cst[ii] = ii;
			case '?':
				if (ii>0) {
					single_strand_cst[ii] = single_strand_cst[ii-1] + 1;
				}
				else {
					single_strand_cst[0] = 1;
				}
				break;
			case '(':
				stack[head] = ii;
				head++;
				single_strand_cst[ii] = 0;
                prevbp='(';
				break;
			case ')':
				openbp = stack[head-1];
				closebp = ii;
				head--;
				if (((prevbp==')')&&(cutpoint==0))||
                    ((prevbp=='(')&&((closebp-openbp-1)>=__hairpin_min_length__)&&((closebp-openbp-1)<=__hairpin_max_length__))||
					((cutpoint>0)&&(((openbp<cutpoint)&&(closebp>cutpoint))))) {
					ss_cst[openbp] = closebp;
					ss_cst[closebp] = openbp;
					single_strand_cst[ii] = 0;
					is_non_empty = 1;
				}
				else {
					fprintf(stderr,"WARNING: Corrupted structure constraints. Minimum or maximum hairpin length requirement is not satisfied. Please, update structure constraints. %s\n",getJobID());
					exit(EXIT_SUCCESS);
				}
                prevbp=')';
				break;
			default:
				fprintf(stderr,"WARNING: Corrupted structure constraints. Constraints ignored. %s\n",getJobID());
				clean_structure_cst_tables(strlen(buffer));
				return;
		}
		
		ii++;
	}
	
#if 0
	for (ii=0;ii<length_seq;ii++) printf("%d,",ss_cst[ii]);
		printf("\n");
	for (ii=0;ii<length_seq;ii++) printf("%d,",single_strand_cst[ii]);
		printf("\n");
#endif


}

/* mutation contraint tables */

int *init_mutation_cst_tables(char *buffer) {
	int *tab,ii=0;
	
	tab=(int *)xmalloc(strlen(buffer)*sizeof(int));
	memset(tab,0,strlen(buffer)*sizeof(int));
	while(buffer[ii]!='\0') {
		
		switch(buffer[ii]) {
			case 'x':
			case 'X':
				tab[ii]=1;
				break;
			case 'o':
			case 'O':
				break;
			default:
				fprintf(stderr,"WARNING: Corrupted mutation contraints. Unrecognized character \"%c\". Ignored.\n",buffer[ii]);
				break;
		}
		ii++;
	}

	return tab;
}

int **init_mutation_cst_segment(int *lintab) {
	int **tab, ii, jj;

	tab=(int **)xmalloc(rna_len*sizeof(int *));
	for (ii=0;ii<rna_len;ii++) {
		tab[ii]=(int *)xmalloc((rna_len-ii)*sizeof(int));
		memset(tab[ii],0,(rna_len-ii)*sizeof(int));
		tab[ii][0] = lintab[ii];
#ifdef SHOWTABLES
		printf("index %d: 0:%d",ii,tab[ii][0]);
#endif
		for (jj=1;jj<rna_len-ii;jj++) {
			tab[ii][jj] = tab[ii][jj-1] + lintab[ii+jj];
#ifdef SHOWTABLES
			printf(", %d:%d",jj,tab[ii][jj]);
#endif
		}
#ifdef SHOWTABLES
		printf("\n");
#endif
	}
	
	return tab;
}

int cmpt_open_mutations() {
	int ii, cmpt=0;
	
	for (ii=0;ii<rna_len;ii++) { if (cst_tape[ii]==0) { cmpt++; } }
	
	return cmpt;
}


