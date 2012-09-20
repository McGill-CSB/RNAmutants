#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<time.h>
#include<float.h>
#include<math.h>
#include<assert.h>
#include"util.h"

#define INDEX_A 0
#define INDEX_C 1
#define INDEX_G 2
#define INDEX_U 3
#define INDEX_GHOST 4
#define INDEX_STOP -1

/* constant defined in RNAloss header */

extern int rna_len;     /* length of the arn */
extern char *input_tape; /* ARN sequence      */
extern const int __hairpin_min_length__; /* the minimal length of a hairpin loop */
extern const int __hairpin_max_length__; /* the maximal length of a hairpin loop */
extern const int __vis_range__;          /* range of visibility */
extern int max_mutations;
extern unsigned long long int total_memory;
extern double mutbias[4][4];
extern int *cst_tape;
extern int **cst_segment;
extern int *ss_cst;
extern int *single_strand_cst;

const int size_pascal_table = 200;

/* generic memory allocation*/

void *xmalloc(int size) {
	void *output = NULL;
	output = malloc(size);
	if (output == NULL) {
		fprintf(stderr,"xmalloc: Cannot proceed to allocate %d byte%s\n",size,getJobID());
		exit(EXIT_FAILURE);
	}
	total_memory += size;
	return output;
}

/* mathematical functions */

inline int minimum(int a, int b) {
	if (b<a) return b;
	else return a;
}

inline double minimum_double(double a, double b) {
	if (b<a) return b;
	else return a;
}

inline int maximum(int a, int b) {
	if (b>a) return b;
	else return a;
}

inline double random_va() {
	static int first_call_random = 1;
	double va; 
	if (first_call_random) {
		srandom(time(0));
		first_call_random = 0;
	}

	/* return va in ]0,1] to avoid border cases */
	va=0;
	while (va==0) { va = (double)(random())/RAND_MAX; }
	
	return va;
}

/* RNA functions */

char *emptyRNAss() {
	int ii;
	char *output = (char *)xmalloc((rna_len+1)*sizeof(char));
	
	for(ii=0;ii<rna_len;ii++) {
		output[ii]='?';
	}
	output[rna_len]='\0';
	
	return output;
}

void cleanRNAss(char *RNAss) {
	int ii;
	
	for(ii=0;ii<rna_len;ii++) {
		RNAss[ii]='.';
	}
	RNAss[rna_len]='\0';
	
}

int canBasePair(int ind1, int ind2) { /* return 1 if a can basepair with b and 0 otherwise */
	
	char base1 = input_tape[ind1];
	char base2 = input_tape[ind2];
	
	if ((base1=='A' && base2=='U')||(base1=='U' && base2=='A')
		||(base1=='G' && base2=='C')||(base1=='C' && base2=='G')
		||(base1=='G' && base2=='U')||(base1=='U' && base2=='G')) {  
		return 1;
	}
	else {
		return 0;
	}
}

int validBasePair(int base1, int base2) { /* return 1 if a can basepair with b and 0 otherwise */
	
	if ((base1==INDEX_A && base2==INDEX_U)||(base1==INDEX_U && base2==INDEX_A)
		||(base1==INDEX_G && base2==INDEX_C)||(base1==INDEX_C && base2==INDEX_G)
		||(base1==INDEX_G && base2==INDEX_U)||(base1==INDEX_U && base2==INDEX_G)) {  
		return 1;
	}
	else {
		return 0;
	}
}

/* check that double are equal under precision range */
/* return 0 if match, non-null otherwise             */

int checkDouble(double a, double b) {
	static int warning_flag = 1;
	long long int v1,v2=0;
	
	if (warning_flag) {
		if (sizeof(double)!=sizeof(long long int)) {
			fprintf(stderr,"WARNING: sizes of checking memory do NOT match%s\n",getJobID());
		}
		warning_flag = 0;
	}
	
	memcpy(&v1,&a,sizeof(long long int));
	memcpy(&v2,&b,sizeof(long long int));
	
	v1 = v1 >> 10; /* extract mantix */
	v2 = v2 >> 10;
	
	return memcmp(&v1,&v2,sizeof(long long int));
	
}

/**** math functions ****/

int kronecker(int index, int nt) { /* return in fact the negation of the kronecker */
	switch (nt) {
		case INDEX_A: if (input_tape[index] == 'A') return 0; else return 1;
		case INDEX_C: if (input_tape[index] == 'C') return 0; else return 1;
		case INDEX_G: if (input_tape[index] == 'G') return 0; else return 1;
		case INDEX_U: if (input_tape[index] == 'U') return 0; else return 1;
		case INDEX_GHOST: return 0;
		default: fprintf(stderr,"kronecker: unexpected index %d%s\n",nt,getJobID()); exit(1);
	}
}

double combinat(int M, int N) {
	static int first_call_combinat = 1;
	static double **pascal;
	
	if ((N>M)||(M<0)||(N<0)) {
#ifdef DEBUG
		fprintf(stderr,"combinat: illegal case (M=%d,N=%d).%s\n",M,N,getJobID());
#endif 
		return 0.0;
	}
	
	if (M > size_pascal_table) {
		fprintf(stderr,"combinat: pascal table size exceed (M=%d,N=%d).%s\n",M,N,getJobID());
		return 0.0;
	}
	
	if (first_call_combinat) {
		int ii,jj;
		
		pascal=(double **)malloc((size_pascal_table+1) * sizeof(double *));
		for (ii=0;ii<=size_pascal_table;ii++) {
			pascal[ii]=(double *)malloc((size_pascal_table+1) * sizeof(double));
			for (jj=0;jj<=size_pascal_table;jj++) {
				pascal[ii][jj] = 0.0;
			}
		}
		
		for (ii=0;ii<=size_pascal_table;ii++) pascal[ii][0] = pascal[ii][ii] = 1;
		
		for (ii=1;ii<=size_pascal_table;ii++)
			for (jj=1;jj<ii;jj++)
				pascal[ii][jj] = pascal[ii-1][jj-1]+pascal[ii-1][jj];
		
		first_call_combinat = 0;
	}
	
	
	return pascal[M][N];
	
}


double genereMutant(int len, int mut) {
	static int first_call_gm = 1;
	static double **gmTable;
	double value;
	
#ifdef CHECK
	if ((len<0)||(mut<0)) fprintf(stderr,"genereMutant: invalid arguments.%s\n",getJobID());
	return 0.0;
#endif
	
	if (first_call_gm) {
		
		int ll,mm,depth = minimum(rna_len,size_pascal_table);
		
		/* init table */
		
		gmTable = (double **)malloc((rna_len+1) * sizeof(double *));
		for (ll=0;ll<=rna_len;ll++) {
			gmTable[ll] = (double *)malloc(max_mutations * sizeof(double));
			for (mm=0;mm<max_mutations;mm++) {
				gmTable[ll][mm] = 0.0;
			}
		}
		
		/* fill table */
		
		for (ll=0;ll<=depth;ll++) {
			for (mm=0;mm<max_mutations;mm++) {
				value = combinat(ll,mm) * pow(3,mm);
				if (value>DBL_EPSILON) {
					gmTable[ll][mm] = value;
				}
			}
		}
		
		first_call_gm = 0;
		
	}
	
	if (len==-1) {
		if (!mut) {
			return 1.0;
		}
		else {
			fprintf(stderr,"%s:%d: Incoherency in the code.%s\n",__FILE__,__LINE__,getJobID());
		}
	}
	
	return gmTable[len][mut];
	
}

/** print configurations from 2-bit representation **/

void print_hcode(int size, unsigned int hc) {
	int i,s=3;
	
	printf("%d",hc&s);
	s=s<<2;
	for(i=1;i<size;i++) {
		printf("%d",(hc&s)>>(2*i));
		s=s<<2;
	}
	fflush(stdout);    
}

int convert2index(char carac) {
	switch (carac) {
		case 'a':
		case 'A': return INDEX_A;
		case 'c':
		case 'C': return INDEX_C;
		case 'g':
		case 'G': return INDEX_G;
		case 'u':
		case 'U': return INDEX_U;
		default :
			fprintf(stderr,"%s:%d: Unknow character (%c)%s\n",__FILE__,__LINE__,carac,getJobID());
			exit(EXIT_FAILURE);
	}
}

inline int char2index(int pos) {
	switch (input_tape[pos]) {
		case 'A': return INDEX_A;
		case 'C': return INDEX_C;
		case 'G': return INDEX_G;
		case 'U': return INDEX_U;
		default :
			fprintf(stderr,"%s:%d: Unknow character (%c) at position %d%s\n",__FILE__,__LINE__,input_tape[pos],pos,getJobID());
			exit(EXIT_FAILURE);
	}
}

char index2char(int pos, int index) {
	char orig = input_tape[pos];
	
	switch (index) {
		case INDEX_A:
			if (orig=='A') return 'A';
			else return 'a';
		case INDEX_C:
			if (orig=='C') return 'C';
			else return 'c';
		case INDEX_G:
			if (orig=='G') return 'G';
			else return 'g';
		case INDEX_U:
			if (orig=='U') return 'U';
			else return 'u';
	}
	return '?';
}

char writeNt(int pos, int lcase) {
	int swValue;
	
	if (lcase) {
		return input_tape[pos];
	}
	else {
		swValue = ((int)(3*random_va())+char2index(pos)+1)%4;
		switch (swValue) {
			case INDEX_A:
				return 'a';
			case INDEX_C:
				return 'c';
			case INDEX_G:
				return 'g';
			case INDEX_U:
				return 'u';
		}
	}
	
	return '?';
}

char writeNt_with_bias(int pos, int lcase) {
	
	double weight[4],weightsum=0,partialsum=0,rValue=0;
	
	if (lcase) {
		return input_tape[pos];
	}
	else {

		weight[INDEX_A] = (char2index(pos)!=INDEX_A) ? mutbias[char2index(pos)][INDEX_A] : 0;
		weight[INDEX_C] = (char2index(pos)!=INDEX_C) ? mutbias[char2index(pos)][INDEX_C] : 0;
		weight[INDEX_G] = (char2index(pos)!=INDEX_G) ? mutbias[char2index(pos)][INDEX_G] : 0;
		weight[INDEX_U] = (char2index(pos)!=INDEX_U) ? mutbias[char2index(pos)][INDEX_U] : 0;
		weightsum = weight[INDEX_A] + weight[INDEX_C] + weight[INDEX_G] + weight[INDEX_U];
		
		rValue = random_va()*weightsum;
		partialsum += weight[INDEX_A];
		if (rValue < partialsum) { return 'a'; }
		partialsum += weight[INDEX_C];
		if (rValue < partialsum) { return 'c'; }
		partialsum += weight[INDEX_G];
		if (rValue < partialsum) { return 'g'; }
		partialsum += weight[INDEX_U];
		if (rValue <= partialsum) { return 'u'; }
		
	}
	
	fprintf(stderr,"WARNING:%s:%d: Single nucleotide sampling failed.%s\n",__FILE__,__LINE__,getJobID());
	
	return '?';
}

static char *job_id = NULL; /* job ID provided for error reporting */

const char *getJobID() {
	return job_id ? job_id : " missing job ID";
}

void setJobID(const char *new_job_id) {
	const char *pfx = " (Job ID: ";
	const char *sfx = ")";
	if (job_id) {
		free(job_id);
		job_id = NULL;
	}
	int len = strlen(pfx) + strlen(new_job_id) + strlen(sfx) + 1;
	job_id = calloc(len, 1);
	assert(job_id);
	snprintf(job_id, len, "%s%s%s", pfx, new_job_id, sfx);
}


/****************************************************************/

/*** mutational bias ***/

inline double nucleotide_bias(int position, int newaa) {
	return mutbias[char2index(position)][newaa];
}

/* enumerate partitions */

double sequence_bias(int lefti, int righti, int nmut) {
	int ii,jj,mm,local_max_mut;
	static double weight[4];
	static double ***biasTable;
	static int first_call_bias = 1;
	
	/* sequence index and mutation numbers start at 0 */
	
	
	if (first_call_bias) {
		
		/* init tables */
		
		weight[INDEX_A] = mutbias[INDEX_A][INDEX_C] + mutbias[INDEX_A][INDEX_G] + mutbias[INDEX_A][INDEX_U];
		weight[INDEX_C] = mutbias[INDEX_C][INDEX_A] + mutbias[INDEX_C][INDEX_G] + mutbias[INDEX_C][INDEX_U];
		weight[INDEX_G] = mutbias[INDEX_G][INDEX_A] + mutbias[INDEX_G][INDEX_C] + mutbias[INDEX_G][INDEX_U];
		weight[INDEX_U] = mutbias[INDEX_U][INDEX_A] + mutbias[INDEX_U][INDEX_C] + mutbias[INDEX_U][INDEX_G];
		
		biasTable = (double ***)malloc((rna_len) * sizeof(double **));
		for (ii=0;ii<rna_len;ii++) {
			biasTable[ii] = (double **)malloc(rna_len * sizeof(double *));
			for (jj=0;jj<rna_len;jj++) {
				biasTable[ii][jj] = (double *)malloc((max_mutations+1) * sizeof(double));
				for (mm=0;mm<=max_mutations;mm++) {
					biasTable[ii][jj][mm] = 0.0;
				}
			}
		}
		
		/* fill table */
		for (ii=0;ii<rna_len;ii++) {
			biasTable[ii][ii][0] = mutbias[char2index(ii)][char2index(ii)];
			if (!cst_tape[ii]) { biasTable[ii][ii][1] = weight[char2index(ii)]; }
			for (jj=ii+1;jj<rna_len;jj++) {
				if (jj-ii+1 >= max_mutations) { local_max_mut=max_mutations; }
				else { local_max_mut=jj-ii+1; }
				for (mm=0;mm<=local_max_mut;mm++) {
					biasTable[ii][jj][mm] += biasTable[ii][jj-1][mm] * mutbias[char2index(jj)][char2index(jj)];
					if ((mm>0)&&(!cst_tape[jj])) {
						biasTable[ii][jj][mm] += biasTable[ii][jj-1][mm-1] * weight[char2index(jj)];
					}
				}
			}
		}
		
		/* setup flag */
		
		first_call_bias = 0;
		
	}
	
	/* return value */
	
	return biasTable[lefti][righti][nmut];
	
}

/****************************************************************/

/*** return mutation partition between two subsequences (NB: return the number of mutations in left segment) ***/

int partition_mutations(int i_ext_left, int i_int_left, int i_int_right, int i_ext_right, int nb_mut_unk_loop) {

	double weight_left=0, weight_right=0,total_weight=0,va=0,anchor;
	int mut_left, mut_right, max_mut_left, min_mut_left, max_mut_right;
	int lenleft = i_int_left-i_ext_left+1;
	int lenright = i_ext_right-i_int_right+1;

	/* border cases */
	
	if ((lenleft<=0)&&(lenright<=0)) { fprintf(stderr,"%s:%d: ERROR. Partition cannot be computed on two empty strings.%s\n",__FILE__,__LINE__,getJobID()); }
	if (lenleft<=0) { return 0; }
	if (lenright<=0) { return nb_mut_unk_loop; }
	
	/* compute boundaries */
	
	if (lenleft>0) { max_mut_left = minimum(nb_mut_unk_loop,lenleft-cst_segment[i_ext_left][lenleft-1]); }
	else { max_mut_left = 0; }
	
	if (lenright>0) { max_mut_right = minimum(nb_mut_unk_loop,lenright-cst_segment[i_int_right][lenright-1]); }
	else { max_mut_right = 0; }

	min_mut_left = maximum(0,nb_mut_unk_loop-max_mut_right);
	
	
	/* compute total weight */
	
	for (mut_left=min_mut_left; mut_left<=max_mut_left; mut_left++) {
		mut_right = nb_mut_unk_loop - mut_left;
		weight_left = sequence_bias(i_ext_left,i_int_left,mut_left);
		weight_right = sequence_bias(i_ext_right,i_int_left,mut_right);
		total_weight += weight_left * weight_right;
	}			
	
	/* sample partition */

	va = random_va() * total_weight;
	anchor = 0;
	for (mut_left=min_mut_left; mut_left<=max_mut_left; mut_left++) {
		mut_right = nb_mut_unk_loop - mut_left;
		weight_left = sequence_bias(i_ext_left,i_int_left,mut_left);
		weight_right = sequence_bias(i_ext_right,i_int_left,mut_right);
		anchor += weight_left * weight_right;
		if (va<=anchor) {
			return mut_left;
		}
	}
	
	/* Should not be reached */
	
	fprintf(stderr,"%s: Partition of mutations failed... (left=[%d,%d], right=[%d,%d], number of mutations=%d, va=%f, partition number=%f)%s\n",
			__FILE__,i_ext_left,i_int_left,i_int_right,i_ext_right,nb_mut_unk_loop,va,anchor,getJobID());
	exit(1);


}

/*** return mutation partition between two subsequences (NB: return the number of mutations in left segment) ***/

int partition_mutations_uniform(int i_ext_left, int i_int_left, int i_int_right, int i_ext_right, int nb_mut_unk_loop) {
	
	double weight_left=0, weight_right=0,total_weight=0,va=0,anchor;
	int mut_left, mut_right, max_mut_left, min_mut_left, max_mut_right;
	int lenleft = i_int_left-i_ext_left+1;
	int lenright = i_ext_right-i_int_right+1;
	
	/* border cases */
	
	if ((lenleft<=0)&&(lenright<=0)) { fprintf(stderr,"%s:%d: ERROR. Partition cannot be computed on two empty strings.%s\n",__FILE__,__LINE__,getJobID()); }
	if (lenleft<=0) { return 0; }
	if (lenright<=0) { return nb_mut_unk_loop; }
	
	/* compute boundaries */
	
	if (lenleft>0) { max_mut_left = minimum(nb_mut_unk_loop,lenleft-cst_segment[i_ext_left][lenleft-1]); }
	else { max_mut_left = 0; }
	
	if (lenright>0) { max_mut_right = minimum(nb_mut_unk_loop,lenright-cst_segment[i_int_right][lenright-1]); }
	else { max_mut_right = 0; }
	
	min_mut_left = maximum(0,nb_mut_unk_loop-max_mut_right);
	
	
	/* compute total weight */
	
	for (mut_left=min_mut_left; mut_left<=max_mut_left; mut_left++) {
		mut_right = nb_mut_unk_loop - mut_left;
		weight_left = genereMutant(lenleft,mut_left);
		weight_right = genereMutant(lenright,mut_right);
		total_weight += weight_left * weight_right;
	}			
	
	/* sample partition */
	
	va = random_va() * total_weight;
	anchor = 0;
	for (mut_left=min_mut_left; mut_left<=max_mut_left; mut_left++) {
		mut_right = nb_mut_unk_loop - mut_left;
		weight_left = genereMutant(lenleft,mut_left);
		weight_right = genereMutant(lenright,mut_right);
		anchor += weight_left * weight_right;
		if (va<=anchor) {
			return mut_left;
		}
	}
	
	/* Should not be reached */
	
	fprintf(stderr,"%s: Partition of mutations failed... (left=[%d,%d], right=[%d,%d], number of mutations=%d, va=%f, partition number=%f)%s\n",
			__FILE__,i_ext_left,i_int_left,i_int_right,i_ext_right,nb_mut_unk_loop,va,anchor,getJobID());
	exit(1);
	
	
}
























