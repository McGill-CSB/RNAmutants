#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <string.h>
#include "reader_energy_params.h"
#include "energy_functions.h"
#include "util.h"
#include "sampling.h"
#include "energy_params_tables.h"
#include "constraints.h"

/* declared in turner_functions.c */

extern double __temperature;
extern double __RT__;

/* in RNAmutants.h */

#define INDEX_A 0
#define INDEX_C 1
#define INDEX_G 2
#define INDEX_U 3
#define INDEX_STOP -1

/* global variables */

struct Zcell {
	double pf;
	double mfe;
};

extern struct Zcell *****Zms;
extern struct Zcell *****Zes;
extern struct Zcell *****Zs;
extern struct Zcell *****Zm;
extern struct Zcell *****Ze;
extern const int __hairpin_min_length__; /* the minimal length of a hairpin loop */
extern const int __hairpin_max_length__; /* the minimal length of a hairpin loop */
extern const int __vis_range__;          /* range of visibility */
extern const int __max_size_bulge__;     /* maximal size of bulge in loop */
extern int rna_len; /* length of the arn*/
extern char *input_tape;
extern int max_mutations;
extern int print_warning;
extern int dangle_used;
extern unsigned int *****loop_tab;
extern unsigned int *****stack_tab;
extern int *cst_tape;
extern int **cst_segment;
extern int *ss_cst;
extern int *single_strand_cst;
extern int nb_value_triloop_table;
extern double *triloop_cmpt_table;
extern double *triloop_weight_table;
extern int nb_value_tetraloop_table;
extern double *tetraloop_cmpt_table;
extern double *tetraloop_weight_table;
extern int uncst_mutation[2][4][4];
extern int *cst_tape;

long int **bp_stat_table=NULL;
long int ***mutation_stat_table=NULL;

/**** variables for hybrid ****/

extern double mfe_hybrid_initialization_energy;
extern double pf_hybrid_initialization_energy;
extern int include_intermolecular_interactions;
extern double cutpoint;

#define NO_TRACEBACK
#define NO_DEBUG_SAMPLING
#define LEFT_CHECK 0
#define RIGHT_CHECK 26

/***********************************************************************************************/
/* fill statistical tables                                                                     */
/***********************************************************************************************/

int fill_stat_tables(char *seq,char *structure) {
	int ii,openbp,closebp,ssa,nt;
	int *tmp_bp_table=NULL, *heap=NULL, itop=0;
	
	/* build base pair table from structure */
	
	tmp_bp_table=(int *)xmalloc(rna_len*sizeof(int));
	for (ii=0;ii<rna_len;ii++) {
		tmp_bp_table[ii]=ii;
	}
	heap=(int *)xmalloc(rna_len*sizeof(int));
	
	for (ii=0;ii<rna_len;ii++) {
		ssa = structure[ii];
		if ((ssa=='(')||(ssa=='[')) {
			heap[itop]=ii+1;
			itop++;
		}
		else if ((ssa==')')||(ssa==']')) {
			openbp=heap[itop];
			closebp=ii;
			itop--;
			tmp_bp_table[openbp]=closebp;
			tmp_bp_table[closebp]=openbp;
			/* fill stat tab */
			bp_stat_table[openbp][closebp]+=1;
			bp_stat_table[closebp][openbp]+=1;		
		}
		else if (ssa=='.') {
			/* fill stat tab */
			bp_stat_table[ii][ii]+=1;
		}
	}
	
#if DEBUG
	if (itop) {
		fprintf(stderr,"WARNING:%s:%s: Corrupted sample structure. Heap not empty.%s\n",__FILE__,__LINE__,getJobID());
	}
#endif
	
	for (ii=0;ii<rna_len;ii++) {
		nt=seq[ii];
		if (nt=='a') {
			mutation_stat_table[ii][INDEX_A][tmp_bp_table[ii]] += 1;
		}
		else if (nt=='c') {
			mutation_stat_table[ii][INDEX_C][tmp_bp_table[ii]] += 1;
		}
		else if (nt=='g') {
			mutation_stat_table[ii][INDEX_G][tmp_bp_table[ii]] += 1;
		}
		else if (nt=='u') {
			mutation_stat_table[ii][INDEX_U][tmp_bp_table[ii]] += 1;
		}
	}
	
	return 1;
	
}

/***********************************************************************************************/
/* fill statistical tables                                                                     */
/***********************************************************************************************/

void print_stat_tables(int nb_mutations, int nb_samples) {
	int ii,jj,nn;
	char nt;
	
	/* print base pair occurrences */
	
	printf("> base pair statistics computed from %d samples with %d mutations\n",nb_samples,nb_mutations);
	for (ii=0;ii<rna_len;ii++) {
		for (jj=0;jj<rna_len;jj++) {
			if ((ii<=jj)&&(bp_stat_table[ii][jj])) {
				printf("%d\t%d\t%ld\n",ii+1,jj+1,bp_stat_table[ii][jj]);
			}
		}
	}
	
	/* print mutation occurrences */
	
	printf("> mutation statistics computed from %d samples with %d mutations\n",nb_samples,nb_mutations);
	for (ii=0;ii<rna_len;ii++) {
		for (nn=0;nn<4;nn++) {
			switch (nn) {
				case 0:
					nt='A';
					break;
				case 1:
					nt='C';
					break;
				case 2:
					nt='G';
					break;
				case 3:
					nt='U';
					break;
				default:
					fprintf(stderr,"%s:%d:corrupted loop.%s\n",__FILE__,__LINE__,getJobID());
					exit(EXIT_FAILURE);
			}
			for (jj=0;jj<rna_len;jj++) {
				if (mutation_stat_table[ii][nn][jj]) {
					printf("%d\t%c\t%d\t%ld\n",ii+1,nt,jj+1,mutation_stat_table[ii][nn][jj]);
				}
			}
		}
	}
}

/***********************************************************************************************/
/* sampling random mutation in a segment                                                       */
/***********************************************************************************************/

int check_sample(char *sample, int expected_nb_mutations) {
	int index=0,effective_nb_mutations=0;
	char carac;
	
	while ((carac=sample[index])!='\0') {
		switch(carac) {
			case 'a':
			case 'c':
			case 'g':
			case 'u':
				effective_nb_mutations++;
		}
		index++;
	}
	
	if (effective_nb_mutations==expected_nb_mutations) {
		return 1;
	}
	else {
		return 0;
	}
}

/***********************************************************************************************/
/* sampling random mutation in a segment                                                       */
/***********************************************************************************************/

void fill_weighted_random_mutations(int nb_mutations, int i, int j, char **ss_sample) {
	int nt, ilist, *nt_list;
	double va,anchor;
	
#ifdef TRACEBACK
	fprintf(stdout,"fill weighted random: mut=%d, ii=%d, jj=%d: %s\n",nb_mutations,i,j,getJobID());
	fflush(stdout);
#endif
	
	/* check */
	
	if (j-i+1<nb_mutations) {
		fprintf(stderr,"WARNING:%s:line %d: Bactrack failed (i=%d,j=%d,mut=%d).%s\n",__FILE__,__LINE__,i,j,nb_mutations,getJobID());
		return;		
	}
	
	/* special case of empty string*/
	if (j<i) { return; }
	
	/* backtrack */
	anchor = 0.0;
	va = random_va() * sequence_bias(i,j,nb_mutations);
	
	if (j-i+1>nb_mutations) { /* there is still some non-mutated positions */
		if (j-i>0) { anchor += sequence_bias(i+1,j,nb_mutations) * nucleotide_bias(i,char2index(i)); } /* end is not reached */
		else { anchor += nucleotide_bias(i,char2index(i)); } /* last nucleotide */
		if (va<=anchor) {
			ss_sample[0][i]=writeNt_with_bias(i,1);
			ss_sample[1][i]='.';
			fill_weighted_random_mutations(nb_mutations, i+1, j, ss_sample);
			return;
		}
	}
	
	if ((nb_mutations>0)&&(!cst_tape[i])) { /* mutations at position i */
		ilist=0;
		nt_list = uncst_mutation[1][char2index(i)];
		while ((nt=nt_list[ilist])!=INDEX_STOP) {
			if (j-i>0) { anchor += nucleotide_bias(i,nt) * sequence_bias(i+1,j,nb_mutations-1); }
			else { anchor += nucleotide_bias(i,nt); }
			if (va<=anchor) {
				ss_sample[0][i]=writeNt_with_bias(i,0);
				ss_sample[1][i]='.';
				fill_weighted_random_mutations(nb_mutations-1, i+1, j, ss_sample);
				return;			
			}
			ilist++;
		}
	}
	
	fprintf(stderr,"%s:line %d: sampling failed (i=%d,j=%d,k=%d,va=%e,anchor=%e,max=%f).%s\n",__FILE__,__LINE__,j,j,nb_mutations,va,anchor,sequence_bias(i,j,nb_mutations),getJobID());
	
#ifdef TRACEBACK
	fprintf(stdout ,"done%s\n",getJobID());
#endif
	
}


/***********************************************************************************************/
/* sampling when (i,j) base pair                                                               */
/***********************************************************************************************/

double samplingHelix(int nb_mut_remaining,int ii,int jj, int lnt, int rnt, char **ss_sample) {
	
	int cloop,nb_mut_motif,nb_nt,nb_max_mut,ii_list,ii_list_open,ii_list_close,config,config_close,config_open;
	int nb_max_mut_unk_loop,nb_mut_unk_loop,nb_min_mut_close_bp,nb_max_mut_close_bp,nb_mut_close_bp,nb_max_mut_open_bp,nb_mut_open_bp;
	int xx,aa,cc,uu,vv,dd,bb,yy,newK,nb_free_nt,len,nn,ll,rr,nb_mut_ext_aa,nb_max_mut_in_bulge;
	int size_bulge,nb_mut_in_bulge,nb_mut_inside_boundary,nt2mut,min_rr,max_rr;
	unsigned int *list_config, *list_config_close, *list_config_open;
	double base_hairpin_energy,E_base_hairpin_energy,bonus,n_total,n_partial,n_remaining,weight,bound_bias;
	int newii,newjj, ilist1=0, *list_nt1=NULL, ilist2=0, *list_nt2=NULL, maxk, nt1, nt2, kk, hcode, valid_random_hairpin;
	double pf_ref_left, pf_ref_right, Esample, Esample1, Esample2;
	int tab_mut_remaining[3][2][2] = {{{0,0}},{{0,1},{1,0}},{{1,1}}}; /* used for tetraloops */
	int min_length;
	
#ifdef DEBUG_SAMPLING
	double debug=0,pdebug=0;
#endif
	
#ifdef TRACEBACK
	fprintf(stdout ,"Helix sampling: mut=%d, ii=%d, jj=%d, lnt=%d, rnt=%d%s\n",nb_mut_remaining,ii,jj,lnt,rnt,getJobID());
#endif
	
	/* In that case, s and b are necessarily equals to 0 */
	
	double va = random_va() * Zs[nb_mut_remaining][ii][jj][lnt][rnt].pf;
	
#ifdef DEBUG_SAMPLING
	if ((ii==LEFT_CHECK)&&(jj==RIGHT_CHECK)) va = Zs[nb_mut_remaining][ii][jj][lnt][rnt].pf;
#endif
	
	/* initialization */
	
	double anchor = 0.0;
	
	/* check */
	
	if (!((ss_cst[ii]==jj)||((ss_cst[ii]<0)&&(ss_cst[jj]<0)))) {
		fprintf(stderr,"%s:line %d:helix sampling for [%d,%d] with %d base pair failed. Structure constraints does not allow base pairing.%s\n",__FILE__,__LINE__,ii,jj,nb_mut_remaining,getJobID());
	}
	
	if (!Zs[nb_mut_remaining][ii][jj][lnt][rnt].pf) {
		fprintf(stderr,"%s:line %d:helix sampling for [%d,%d] with %d base pair failed.%s\n",__FILE__,__LINE__,ii,jj,nb_mut_remaining,getJobID());
	}
	
	if (((jj<cutpoint)||(ii>=cutpoint))&&(jj-ii-1 < __hairpin_min_length__)) {
		fprintf(stderr,"%s:line %d:helix sampling failed... Hairpin threshold is not respected at position (%d,%d)%s\n",__FILE__,__LINE__,ii,jj,getJobID());
	}
	
	if ((cst_tape[ii])&&(kronecker(ii,lnt))) {
		fprintf(stderr,"%s:line %d:%s: %d: helix sampling failed... index %d cannot mutate%s\n",__FILE__,__LINE__,__FILE__,__LINE__,ii,getJobID());
		return 0;
	}
	
	if ((cst_tape[jj])&&(kronecker(jj,rnt))) {
		fprintf(stderr,"%s:line %d: helix sampling failed... index %d cannot mutate%s\n",__FILE__,__LINE__,jj,getJobID());
		return 0;
	}
	
	/* include base pairs */
	
	
	ss_sample[0][ii]=index2char(ii,lnt);
	ss_sample[0][jj]=index2char(jj,rnt);
	
	if ((ii<cutpoint)&&(jj>cutpoint)) { /* inter molecular hairpin */
		if (include_intermolecular_interactions) {
			ss_sample[1][ii]='[';
			ss_sample[1][jj]=']';
		}
		else {
			ss_sample[1][ii]='{';
			ss_sample[1][jj]='}';
		}
	}
	else {		
		ss_sample[1][ii]='(';
		ss_sample[1][jj]=')';
	}
	
	/* length */
	
	len = jj - ii + 1;
	
	/* number of mutation in external nucleotides */
	
	nb_mut_ext_aa = kronecker(ii,lnt) + kronecker(jj,rnt);
	nb_mut_inside_boundary = nb_mut_remaining - nb_mut_ext_aa;		  
	
	/*****************************************************************************************************/
	/** hairpin                                                                                         **/
	/*****************************************************************************************************/
	
	if (is_hairpin(ii,jj)) {
		if ((ii<cutpoint)&&(jj>cutpoint)) { /* inter molecular hairpin */
			if (include_intermolecular_interactions) {
				newK = nb_mut_remaining-kronecker(ii,lnt)-kronecker(jj,rnt);
				if (len>2) {
					nt2mut = len-2 - cst_segment[ii+1][jj-ii-2];
				}
				else {
					nt2mut=0;
				}
				
				if ((newK>=0)||(newK<=nt2mut)) { /* check configurations */
					
					int nmut_left,nmut_right;
					double pf_left,pf_right;
					int rightmostseq1=(int)(cutpoint-0.5);
					int leftmostseq2=(int)(cutpoint+0.5);
					int valid_index_left,valid_index_right,valid_index;
					
					/* single stranded intermolecular hairpin */
					
					/* n_remaining = genereMutant(nt2mut,newK); */
					if (jj-ii>1) {
						n_remaining = sequence_bias(ii+1,jj-1,newK) * nucleotide_bias(ii,lnt) * nucleotide_bias(jj,rnt);
					}
					else {
						n_remaining = nucleotide_bias(ii,lnt) * nucleotide_bias(jj,rnt);
					}
					anchor += n_remaining * pf_hybrid_initialization_energy;
					
					if ((n_remaining)&&(va<=anchor)) {
						
#ifdef TRACEBACK
						fprintf(stdout ,"line %d: sample hairpin: k=%d, i=%d, j=%d, va=%e, anchor=%e%s\n",__LINE__,newK,ii,jj,va,anchor,getJobID());
#endif
						if (len>2) {
							fill_weighted_random_mutations(newK,ii+1, jj-1, ss_sample);
						}
						return mfe_hybrid_initialization_energy;
					}
					
					/* intermolecular hairpin with secondary structure */
					
					valid_index_left=0;
					valid_index_right=0;
					valid_index=0;
					
					if (rightmostseq1-ii>=__hairpin_min_length__+2) {
						valid_index_left=1;
						valid_index=1;
					}
					if (jj-leftmostseq2>=__hairpin_min_length__+2) {
						valid_index_right=1;
						valid_index=1;
					}
					
					if (valid_index) {
						
						for (nmut_left=0;nmut_left<=newK;nmut_left++) {
							int ss,uu,vv,ww;
							nmut_right=newK-nmut_left;
							
							for(uu=0;uu<4;uu++) {
								for(ss=0;ss<4;ss++) {
									/* if: rightmostseq1 - ii < 3, then sequence_bias returns 0 */
									if (rightmostseq1 - ii >= 3) { 
										pf_ref_left = pf_left = sequence_bias(ii+2,rightmostseq1-1,nmut_left-kronecker(ii+1,uu)-kronecker(rightmostseq1,ss))
										* nucleotide_bias(ii+1,uu) * nucleotide_bias(rightmostseq1,ss);
									}
									else {
										pf_ref_left = pf_left = nucleotide_bias(ii+1,uu) * nucleotide_bias(rightmostseq1,ss);
									}
									if (valid_index_left) {
										pf_left+=Zes[nmut_left][ii+1][rightmostseq1][uu][ss].pf+Ze[nmut_left][ii+1][rightmostseq1][uu][ss].pf;
									}
									for(ww=0;ww<4;ww++) {
										for(vv=0;vv<4;vv++) {
											if (jj - leftmostseq2 >= 3) {
												pf_ref_right = pf_right = sequence_bias(leftmostseq2+1,jj-2,nmut_right-kronecker(leftmostseq2,ww)-kronecker(jj-1,vv))
												* nucleotide_bias(leftmostseq2,ww) * nucleotide_bias(jj-1,vv);
											}
											else {
												pf_ref_right = pf_right = nucleotide_bias(leftmostseq2,ww) * nucleotide_bias(jj-1,vv);
											}
											if (valid_index_right) {
												pf_right+=Zes[nmut_right][leftmostseq2][jj-1][ww][vv].pf+Ze[nmut_right][leftmostseq2][jj-1][ww][vv].pf;
											}
											
											if ((pf_left>pf_ref_left)||(pf_right>pf_ref_right)) {
												
												anchor += pf_left * pf_right * pf_hybrid_initialization_energy * nucleotide_bias(ii,lnt) * nucleotide_bias(jj,rnt);
												
												if (va<=anchor) {
													
#ifdef TRACEBACK
													fprintf(stdout ,"line %d: sample intermolecular hairpin: k=%d, i=%d, j=%d, u=%d, s=%d, w=%d, v=%d, va=%e, anchor=%e%s\n", __LINE__, newK, ii, jj, uu, ss, ww, vv, va, anchor,getJobID());
#endif
													
													Esample1 = 0.0;
													Esample2 = 0.0;
													if (pf_left>pf_ref_left) { /* non empty secndary structure in left region */
														Esample1 = samplingExteriorLoop(nmut_left,ii+1,rightmostseq1,uu,ss,ss_sample);
													}
													else { /* empty secndary structure in left region */													
														if (ii<rightmostseq1) {
															fill_weighted_random_mutations(nmut_left,ii+1,rightmostseq1, ss_sample);
														}
													}
													if (pf_right>pf_ref_right) { /* non empty secndary structure in left region */
														Esample2 = samplingExteriorLoop(nmut_right,leftmostseq2,jj-1,ww,vv,ss_sample);
													}
													else{ /* empty secndary structure in left region */
														if (leftmostseq2>jj) {
															fill_weighted_random_mutations(nmut_right,leftmostseq2,jj-1, ss_sample);
														}
													}
													return mfe_hybrid_initialization_energy + Esample1 + Esample2;
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
		else { /* single strand hairpin & check contraints */
			
			for (uu=0;uu<4;uu++) {
				if (!((cst_tape[ii+1])&&(kronecker(ii+1,uu)))) {
					for (vv=0;vv<4;vv++) {
						if (!((cst_tape[jj-1])&&(kronecker(jj-1,vv)))) {
							
							base_hairpin_energy = boltzmannHairpin(lnt,uu,vv,rnt,len-2);
							E_base_hairpin_energy = EHairpin(lnt,uu,vv,rnt,len-2);
							newK = nb_mut_remaining-kronecker(ii,lnt)-kronecker(ii+1,uu)-kronecker(jj-1,vv)-kronecker(jj,rnt);
							nt2mut = len-4 - cst_segment[ii+2][jj-ii-4];
							bound_bias =  nucleotide_bias(ii,lnt) * nucleotide_bias(ii+1,uu) * nucleotide_bias(jj-1,vv) * nucleotide_bias(jj,rnt);
							
							if ((newK<0)||(newK>nt2mut)) continue; /* check configurations */
							
							n_total = n_partial = 0;
							
							if (len==5) { /*** triloop ***/
								
								/* assign new value */
								if ((newK==0)||(!cst_tape[ii+2])) {
									ilist1=0;
									list_nt1=uncst_mutation[newK][char2index(ii+2)];
									while ((nt1=list_nt1[ilist1])!=INDEX_STOP) {
										hcode = (lnt<<8) | (uu<<6) | (nt1<<4) | (vv<<2) | rnt;
										if ((bonus=refTableTriLoop[hcode])) { /* triloop exists if bonus exists */										
											n_partial = bound_bias * nucleotide_bias(ii+2,nt1);
											n_total += n_partial;
											anchor +=  base_hairpin_energy * exp(-(bonus/__RT__)) * n_partial;
											if (va<=anchor) {
#ifdef TRACEBACK
												fprintf(stdout ,"line %d: sample triloop: k=%d, i=%d, j=%d, u=%d, nt1=%d, v=%d, va=%e, anchor=%e%s\n",__LINE__,newK,ii,jj,uu,nt1,vv,va,anchor,getJobID());
#endif
												ss_sample[0][ii+1]=index2char(ii+1,uu);
												ss_sample[1][ii+1]='.';
												ss_sample[0][ii+2]=index2char(ii+2,nt1);
												ss_sample[1][ii+2]='.';
												ss_sample[0][jj-1]=index2char(jj-1,vv);
												ss_sample[1][jj-1]='.';
												return E_base_hairpin_energy + bonus;
											}
										}
										ilist1++;
									}
								}
							}
							else if (len==6) { /*** tetraloop ***/
								
								if (newK==1) { maxk=2; } /* used to decide the number of configurations to try */
								else { maxk=1; }
								
								for (kk=0;kk<maxk;kk++) {
									if ((!tab_mut_remaining[newK][kk][0])||(!cst_tape[ii+2])) { /* continue if mutation allowed or unconstrained position */
										ilist1=0;
										list_nt1=uncst_mutation[tab_mut_remaining[newK][kk][0]][char2index(ii+2)];
										while ((nt1=list_nt1[ilist1])!=INDEX_STOP) {
											if ((!tab_mut_remaining[newK][kk][1])||(!cst_tape[ii+3])) { /* continue if mutation allowed or unconstrained position */
												ilist2=0;
												list_nt2=uncst_mutation[tab_mut_remaining[newK][kk][1]][char2index(ii+3)];
												while ((nt2=list_nt2[ilist2])!=INDEX_STOP) {
													hcode = (lnt<<10) | (uu<<8) | (nt1<<6) | (nt2<<4) | (vv<<2) | rnt;
													if ((bonus=refTableTetraLoop[hcode])) {
														n_partial = bound_bias * nucleotide_bias(ii+2,nt1) * nucleotide_bias(ii+3,nt2);
														n_total += n_partial;
														anchor +=  base_hairpin_energy * exp(-(bonus/__RT__)) * n_partial;
														if (va<=anchor) {
#ifdef TRACEBACK
															fprintf(stdout ,"line %d: sample tetraloop: k=%d, i=%d, j=%d, u=%d, nt1=%d, nt2=%d, v=%d, va=%e, anchor=%e%s\n",__LINE__,newK,ii,jj,uu,nt1,nt2,vv,va,anchor,getJobID());
#endif
															ss_sample[0][ii+1]=index2char(ii+1,uu);
															ss_sample[1][ii+1]='.';
															ss_sample[0][ii+2]=index2char(ii+2,nt1);
															ss_sample[1][ii+2]='.';
															ss_sample[0][ii+3]=index2char(ii+3,nt2);
															ss_sample[1][ii+3]='.';
															ss_sample[0][jj-1]=index2char(jj-1,vv);
															ss_sample[1][jj-1]='.';
															return  E_base_hairpin_energy + bonus;
														}
													}
													ilist2++;
												}
											}
											ilist1++;
										}
									}
								}
							}
							
							/*** special case ***/
							
#ifdef WITH_GGG_HAIRPIN
							if((lnt==INDEX_G)&&(uu==INDEX_G)&&(rnt==INDEX_U)) { /*** GGG loop ***/
								double ggg_newK = newK;
								n_partial = 0;
								if (kronecker(ii+2,INDEX_G)) { /*** 3rd nt needs to mutate ***/
									if ((newK)&&(!(cst_tape[ii+2]))) {
										ggg_newK--; /** remove the mutation for GGG loop **/
										/* n_partial = genereMutant(nt2mut-1,ggg_newK); */
										n_partial = sequence_bias(ii+2,jj-1,ggg_newK) * bound_bias;
										n_total += n_partial;
										anchor +=  base_hairpin_energy * exp(-(ggg_hairpin_bonus/__RT__)) * n_partial;
									}
								}
								else { /*** already a GGG loop ***/
									/* n_partial = genereMutant(nt2mut-1,ggg_newK); */
									n_partial = sequence_bias(ii+2,jj-1,ggg_newK) * bound_bias;
									n_total += n_partial;
									anchor +=  base_hairpin_energy * exp(-(ggg_hairpin_bonus/__RT__)) * n_partial;
								}
								if (va<=anchor) {
#ifdef TRACEBACK
									fprintf(stdout ,"line %d: sample GGG loop: k=%d, i=%d, j=%d, u=%d, v=%d, va=%e, anchor=%e%s\n",__LINE__,newK,ii,jj,uu,vv,va,anchor,getJobID());
#endif
									if ((len-2<__hairpin_min_length__)||(len-2>__hairpin_max_length__)) {
										fprintf(stderr,"Incoherency: Hairpin does not have the correct length!! (line %d)%s\n",__LINE__,getJobID());
										return 0;
									}
									ss_sample[0][ii+1]=index2char(ii+1,uu);
									ss_sample[1][ii+1]='.';
									ss_sample[0][ii+2]=index2char(ii+2,INDEX_G);
									ss_sample[1][ii+2]='.';
									ss_sample[0][jj-1]=index2char(jj-1,vv);
									ss_sample[1][jj-1]='.';
									if (len>5) { /*** fill region if the size of the hairpin is more than 3 ***/
										fill_weighted_random_mutations(ggg_newK,ii+3, jj-2, ss_sample);
									}
									return E_base_hairpin_energy + ggg_hairpin_bonus;
								}	  
							}
#endif
							/*** generic case ***/
							/* n_remaining = genereMutant(nt2mut,newK) - n_total; */
							n_remaining = sequence_bias(ii+2,jj-2,newK) * bound_bias - n_total;
							anchor +=  base_hairpin_energy * n_remaining;
							if (va<=anchor) {
#ifdef TRACEBACK
								fprintf(stdout ,"line %d: sample hairpin: k=%d, i=%d, j=%d, u=%d, v=%d, va=%e, anchor=%e%s\n",__LINE__,newK,ii,jj,uu,vv,va,anchor,getJobID());
#endif
								if ((len-2<__hairpin_min_length__)||(len-2>__hairpin_max_length__)) {
									fprintf(stderr,"Incoherency: Hairpin does not have the correct length!! (line %d)%s\n",__LINE__,getJobID());
									return 0;
								}
								
								/* filter tri- and tetra-loop (not rigorous) */
								valid_random_hairpin=0;
								while (!valid_random_hairpin) {
									
									ss_sample[0][ii+1]=index2char(ii+1,uu);
									ss_sample[1][ii+1]='.';
									ss_sample[0][jj-1]=index2char(jj-1,vv);
									ss_sample[1][jj-1]='.';
									fill_weighted_random_mutations(newK,ii+2, jj-2, ss_sample);
									
#ifdef WITH_GGG_HAIRPIN
									if ((!is_triloop(ii,jj,ss_sample))&&(!is_tetraloop(ii,jj,ss_sample))&&(!is_ggg_hairpin(ii,jj,ss_sample))) {
										valid_random_hairpin=1;
									}
#else
									if ((!is_triloop(ii,jj,ss_sample))&&(!is_tetraloop(ii,jj,ss_sample))) {
										valid_random_hairpin=1;
									}
#endif								
								}
								return E_base_hairpin_energy;
								
							}	  
						}
					}
				}
			}
		}
	}
	
	/* a simple check to ensure we can insert a regular loop (assert minimal length requirements) */
	
	if ((((ii<cutpoint)&&(cutpoint<jj))&&(len<4))||(((ii>=cutpoint)||(cutpoint>jj))&&(len<__hairpin_min_length__+2))) {
		fprintf(stderr,"%s:%d: ERROR: Cannot sample a regular loop between index %d and %d.%s\n",__FILE__,__LINE__,ii,jj,getJobID());
		return 0;
	}
	
	
	/*****************************************************************************************************/
	/** Stacking pairs                                                                                  **/
	/*****************************************************************************************************/
	
	if (is_stack(ii,jj)) {
		
		nb_max_mut=minimum(4-cst_tape[ii]-cst_tape[ii+1]-cst_tape[jj-1]-cst_tape[jj],nb_mut_remaining);
		
		for (nb_mut_motif=nb_mut_ext_aa;nb_mut_motif<=nb_max_mut;nb_mut_motif++) {
			list_config = stack_tab[nb_mut_motif][ii][jj][0];
			ii_list = 0;
			while ((config=list_config[ii_list])) {
				xx=config>>6;
				yy=config&3;
				
				if ((lnt==xx)&&(rnt==yy)) {
					uu=(config>>4)&3;
					vv=(config>>2)&3;
					/* we must add the cost of external nucleotides */
					newK=nb_mut_remaining-nb_mut_motif + kronecker(ii+1,uu) + kronecker(jj-1,vv);
					if ((newK<0)||(newK>nb_mut_remaining)) {
						fprintf(stderr,"Incoherency: newK CANNOT be negative!! (line %d)%s\n",__LINE__,getJobID());
						return 0;
					}
					anchor += boltzmannStack(lnt,uu,vv,rnt) * Zs[newK][ii+1][jj-1][uu][vv].pf * nucleotide_bias(ii,lnt) * nucleotide_bias(jj,rnt);
#ifdef DEBUG_SAMPLING
					if ((ii==LEFT_CHECK)&&(jj==RIGHT_CHECK)) debug += boltzmannStack(lnt,uu,vv,rnt) * Zs[newK][ii+1][jj-1][uu][vv].pf;
#endif
					if (va<=anchor) {
#ifdef TRACEBACK
						fprintf(stdout ,"line %d: sample stack: k=%d, i=%d, j=%d, u=%d, v=%d, va=%e, anchor=%e%s\n",__LINE__,newK,ii+1,jj-1,uu,vv,va,anchor,getJobID());
#endif
						Esample = samplingHelix(newK,ii+1,jj-1,uu,vv,ss_sample);
						return EStack(lnt,uu,vv,rnt) + Esample;
					}
				}
				ii_list++;
			}
		}
		
#ifdef DEBUG_SAMPLING
		if ((ii==LEFT_CHECK)&&(jj==RIGHT_CHECK)) printf("sample stack: %e\n",debug-pdebug);
		pdebug = debug;
#endif
		
	}
	
	/*****************************************************************************************************/
	/** bulge AND multi-loop                                                                            **/
	/*****************************************************************************************************/
	
	nb_mut_inside_boundary = nb_mut_remaining - kronecker(ii,lnt) - kronecker(jj,rnt);
	
	if (nb_mut_inside_boundary>=0) {
		
		if ((ii<cutpoint)&&(cutpoint<jj)) {
			nb_free_nt = minimum(__max_size_bulge__,len - 4);
			min_length = 4;
		} else {
			nb_free_nt = minimum(__max_size_bulge__,len - __hairpin_min_length__ - 4);
			min_length = __hairpin_min_length__ +4;
		}
		
		if (len > min_length) { /* min length required for a stack */
			
			for (uu=0;uu<4;uu++) {
				
				for (vv=0;vv<4;vv++) {
					
					if (validBasePair(uu,vv)) {
						
						/* bulge */
						
						for (size_bulge=1;size_bulge<=nb_free_nt;size_bulge++) { /* size of the bulge */
							
							/* we do need to distinguish the bulge because the asymetry of constraints */
							
							/* bulge on the left */
							
							newii = ii+1+size_bulge;
							newjj = jj-1;
							
							if (((newii<cutpoint)&&(cutpoint<newjj))||((len>=__hairpin_min_length__+4+size_bulge)&&((newjj<cutpoint)||(newii>=cutpoint)))) { /* check that loop respects cutpoint */
								
								if (is_internal_loop(ii,newii,newjj,jj)) { /* lookup structure constraints */						
									
									if (!(kronecker(newii,uu) && cst_tape[newii])||
										!(kronecker(newjj,vv) && cst_tape[newjj])) {
										
										nt2mut=size_bulge-cst_segment[ii+1][size_bulge-1];
										nb_max_mut_in_bulge = minimum(nt2mut,nb_mut_inside_boundary);
										
										for (nb_mut_in_bulge=0;nb_mut_in_bulge<=nb_max_mut_in_bulge;nb_mut_in_bulge++) { /* number of mutations in bulge */
											newK = nb_mut_inside_boundary  - nb_mut_in_bulge;
											if (newK>=0) {
												if (Zs[newK][newii][newjj][uu][vv].pf) {
													
													/* weight = genereMutant(nt2mut,nb_mut_in_bulge) */
													weight = sequence_bias(ii+1,ii+size_bulge,nb_mut_in_bulge);
													
													anchor += boltzmannBulge(lnt,uu,vv,rnt,size_bulge) * Zs[newK][newii][newjj][uu][vv].pf * weight * nucleotide_bias(ii,lnt) * nucleotide_bias(jj,rnt);
													
#ifdef DEBUG_SAMPLING
													if ((ii==LEFT_CHECK)&&(jj==RIGHT_CHECK)) debug += boltzmannBulge(lnt,uu,vv,rnt,size_bulge) *
														Zs[newK][newii][newjj][uu][vv].pf * genereMutant(nt2mut,nb_mut_in_bulge);
#endif
													if (va<=anchor) {
#ifdef TRACEBACK
														fprintf(stdout ,"line %d: sample bulge: k=%d, i=%d, j=%d, u=%d, v=%d, va=%e, anchor=%e%s\n",__LINE__,newK,newii,newjj,uu,vv,va,anchor,getJobID());
#endif
														fill_weighted_random_mutations(nb_mut_in_bulge,ii+1,ii+size_bulge,ss_sample);
														Esample = samplingHelix(newK,newii,newjj,uu,vv,ss_sample);
														return EBulge(lnt,uu,vv,rnt,size_bulge) + Esample;
													}
												}
											}
										}
									}
								}
							}
							
							/* bulge on the right */
							
							newii = ii+1;
							newjj = jj-1-size_bulge;
							
							if (((newii<cutpoint)&&(cutpoint<newjj))||((len>=__hairpin_min_length__+4+size_bulge)&&((newjj<cutpoint)||(newii>=cutpoint)))) { /* check that loop respects cutpoint */
								
								if (is_internal_loop(ii,newii,newjj,jj)) { /* lookup structure constraints */
									
									if (!(kronecker(newii,uu) && cst_tape[newii])||
										!(kronecker(newjj,vv) && cst_tape[newjj])) {
										
										nt2mut=size_bulge-cst_segment[jj-size_bulge][size_bulge-1];
										nb_max_mut_in_bulge = minimum(nt2mut,nb_mut_inside_boundary);
										
										for (nb_mut_in_bulge=0;nb_mut_in_bulge<=nb_max_mut_in_bulge;nb_mut_in_bulge++) { /* number of mutations in bulge */
											newK = nb_mut_inside_boundary  - nb_mut_in_bulge;
											if (newK>=0) {
												if (Zs[newK][newii][newjj][uu][vv].pf) {
													
													/* weight = genereMutant(nt2mut,nb_mut_in_bulge) */
													weight = sequence_bias(jj-size_bulge,jj-1,nb_mut_in_bulge);
													
													anchor += boltzmannBulge(lnt,uu,vv,rnt,size_bulge) * Zs[newK][newii][newjj][uu][vv].pf * weight * weight * nucleotide_bias(ii,lnt) * nucleotide_bias(jj,rnt);
#ifdef DEBUG_SAMPLING
													if ((ii==LEFT_CHECK)&&(jj==RIGHT_CHECK)) debug += boltzmannBulge(lnt,uu,vv,rnt,size_bulge) *
														Zs[newK][newii][newjj][uu][vv].pf * genereMutant(nt2mut,nb_mut_in_bulge);
#endif
													if (va<=anchor) {
#ifdef TRACEBACK
														fprintf(stdout ,"line %d: sample bulge: k=%d, i=%d, j=%d, u=%d, v=%d, va=%e, anchor=%e%s\n",__LINE__,newK,newii,newjj,uu,vv,va,anchor,getJobID());
#endif
														fill_weighted_random_mutations(nb_mut_in_bulge,jj-size_bulge,jj-1,ss_sample);
														Esample = samplingHelix(newK,newii,newjj,uu,vv,ss_sample);
														return EBulge(lnt,uu,vv,rnt,size_bulge) + Esample;
													}
												}
											}
										}
									}	      
								}
							}
						}
					}
					
					/* close multi-loop*/
					
					anchor += boltzmannEndMultiLoop() * Zm[nb_mut_inside_boundary][ii+1][jj-1][uu][vv].pf * nucleotide_bias(ii,lnt) * nucleotide_bias(jj,rnt);
#if 0
					if (Zm[nb_mut_inside_boundary][ii+1][jj-1][uu][vv].pf)
						printf("%d,%d: %e (%e * %e)\n",ii,jj,Zs[nb_mut_remaining][ii][jj][lnt][rnt].pf,boltzmannEndMultiLoop(),Zm[nb_mut_inside_boundary][ii+1][jj-1][uu][vv].pf);
#endif
					
#ifdef DEBUG_SAMPLING
					if ((ii==LEFT_CHECK)&&(jj==RIGHT_CHECK)) debug += boltzmannEndMultiLoop() * Zm[nb_mut_inside_boundary][ii+1][jj-1][uu][vv].pf;
#endif
					
					if (va<=anchor) {
#ifdef TRACEBACK
						fprintf(stdout ,"line %d: sample multi-loop: k=%d, i=%d, j=%d, u=%d, v=%d, va=%e, anchor=%e%s\n",__LINE__,nb_mut_inside_boundary,ii+1,jj-1,uu,vv,va,anchor,getJobID());
#endif
						Esample = samplingMultiLoop(nb_mut_inside_boundary,ii+1,jj-1,uu,vv,0,ss_sample);
						return EEndMultiLoop() + Esample;
					}	      	  
				}
			}
		}
	}
	
#ifdef DEBUG_SAMPLING
	if ((ii==LEFT_CHECK)&&(jj==RIGHT_CHECK)) printf("sample bulge & multi: %e\n",debug-pdebug);
	pdebug = debug;
#endif
	
	/*****************************************************************************************************/
	/** lookup special internal loops                                                                   **/
	/*****************************************************************************************************/
	
	for (cloop=0;cloop<4;cloop++) {
		
		nb_nt = 6 + (cloop>>1) + (cloop&1); /* number of nt in the internal loop */
		nb_max_mut=minimum(nb_nt,nb_mut_remaining);
		
		for (nb_mut_motif=0;nb_mut_motif<=nb_max_mut;nb_mut_motif++) {
			list_config = loop_tab[nb_mut_motif][cloop][ii][jj];
			ii_list = 0;
			while ((config=list_config[ii_list])) {
				switch (cloop) {
					case 0:
						newii=ii+2;
						newjj=jj-2;						
						if (!is_internal_loop(ii,newii,newjj,jj)) break; /* lookup structure constraints */
						if (((newii<cutpoint)&&(newjj>cutpoint))||((len>=__hairpin_min_length__+6)&&((newjj<cutpoint)||(newii>=cutpoint)))) { /* constrain loop to respect cutpoint */
							xx=config>>10;
							yy=config&3;
							if ((xx==lnt)&&(yy==rnt)) {
								aa=(config>>8)&3;
								uu=(config>>6)&3;
								vv=(config>>4)&3;
								bb=(config>>2)&3;
								newK=nb_mut_remaining-nb_mut_motif + kronecker(ii+2,uu) + kronecker(jj-2,vv);
								if ((newK<0)||(newK>nb_mut_remaining)) {
									fprintf(stderr,"Incoherency: newK CANNOT be negative!! (line %d)%s\n",__LINE__,getJobID());
									return 0;
								}
								anchor += boltzmannInternal_1x1(lnt,aa,uu,vv,bb,rnt) * Zs[newK][ii+2][jj-2][uu][vv].pf * nucleotide_bias(ii,lnt) * nucleotide_bias(ii+1,aa) * nucleotide_bias(jj-1,bb) * nucleotide_bias(jj,rnt);
#ifdef DEBUG_SAMPLING
								if ((ii==LEFT_CHECK)&&(jj==RIGHT_CHECK)) debug += boltzmannInternal_1x1(lnt,aa,uu,vv,bb,rnt) * Zs[newK][ii+2][jj-2][uu][vv].pf;
#endif
								if (va<=anchor) {
									ss_sample[0][ii+1] = index2char(ii+1,aa);
									ss_sample[1][ii+1] = '.';
									ss_sample[0][jj-1] = index2char(jj-1,bb);
									ss_sample[1][jj-1] = '.';
#ifdef TRACEBACK
									fprintf(stdout ,"line %d: sample internal-loop: k=%d, i=%d, j=%d, u=%d, v=%d, va=%e, anchor=%e%s\n",__LINE__,newK,ii+2,jj-2,uu,vv,va,anchor,getJobID());
#endif
									Esample = samplingHelix(newK,ii+2,jj-2,uu,vv,ss_sample);
									return EInternal_1x1(lnt,aa,uu,vv,bb,rnt) + Esample;
								}	      
							}
						}
						break;
					case 1:
						newii=ii+2;
						newjj=jj-3;
						if (!is_internal_loop(ii,newii,newjj,jj)) break; /* lookup structure constraints */
						if (((newii<cutpoint)&&(newjj>cutpoint))||((len>=__hairpin_min_length__+7)&&((newjj<cutpoint)||(newii>=cutpoint)))) { /* constrain loop to respect cutpoint */
							xx=config>>12;
							yy=config&3;
							if ((xx==lnt)&&(yy==rnt)) {
								aa=(config>>10)&3;
								uu=(config>>8)&3;
								vv=(config>>6)&3;
								cc=(config>>4)&3;
								bb=(config>>2)&3;
								newK=nb_mut_remaining-nb_mut_motif + kronecker(ii+2,uu) + kronecker(jj-3,vv);
								if ((newK<0)||(newK>nb_mut_remaining)) {
									fprintf(stderr,"Incoherency: newK CANNOT be negative!! (line %d)%s\n",__LINE__,getJobID());
									return 0;
								}
								anchor += boltzmannInternal_1x2(lnt,aa,uu,vv,cc,bb,rnt) * Zs[newK][ii+2][jj-3][uu][vv].pf * nucleotide_bias(ii,lnt) * nucleotide_bias(ii+1,aa) * nucleotide_bias(jj-2,cc) * nucleotide_bias(jj-1,bb) * nucleotide_bias(jj,rnt);
#ifdef DEBUG_SAMPLING
								if ((ii==LEFT_CHECK)&&(jj==RIGHT_CHECK)) debug += boltzmannInternal_1x2(lnt,aa,uu,vv,cc,bb,rnt) * Zs[newK][ii+2][jj-3][uu][vv].pf;
#endif
								if (va<=anchor) {
									ss_sample[0][ii+1] = index2char(ii+1,aa);
									ss_sample[1][ii+1] = '.';
									ss_sample[0][jj-2] = index2char(jj-2,cc);
									ss_sample[1][jj-2] = '.';
									ss_sample[0][jj-1] = index2char(jj-1,bb);
									ss_sample[1][jj-1] = '.';
#ifdef TRACEBACK
									fprintf(stdout ,"line %d: sample internal-loop: k=%d, i=%d, j=%d, u=%d, v=%d, va=%e, anchor=%e%s\n",__LINE__,newK,ii+2,jj-3,uu,vv,va,anchor,getJobID());
#endif
									Esample = samplingHelix(newK,ii+2,jj-3,uu,vv,ss_sample);
									return EInternal_1x2(lnt,aa,uu,vv,cc,bb,rnt) + Esample;
								}	      
							}
						}
						break;
					case 2:
						newii=ii+3;
						newjj=jj-2;
						if (!is_internal_loop(ii,newii,newjj,jj)) break; /* lookup structure constraints */
						if (((newii<cutpoint)&&(newjj>cutpoint))||((len>=__hairpin_min_length__+7)&&((newjj<cutpoint)||(newii>=cutpoint)))) { /* constrain loop to respect cutpoint */
							xx=config>>12;
							yy=config&3;
							if ((xx==lnt)&&(yy==rnt)) {
								aa=(config>>10)&3;
								cc=(config>>8)&3;
								uu=(config>>6)&3;
								vv=(config>>4)&3;
								bb=(config>>2)&3;
								newK=nb_mut_remaining-nb_mut_motif + kronecker(ii+3,uu) + kronecker(jj-2,vv);
								if ((newK<0)||(newK>nb_mut_remaining)) {
									fprintf(stderr,"Incoherency: newK CANNOT be negative!! (line %d)%s\n",__LINE__,getJobID());
									return 0;
								}
								anchor += boltzmannInternal_2x1(lnt,aa,cc,uu,vv,bb,rnt) * Zs[newK][ii+3][jj-2][uu][vv].pf * nucleotide_bias(ii,lnt) * nucleotide_bias(ii+1,aa) * nucleotide_bias(ii+2,cc) * nucleotide_bias(jj-1,bb) * nucleotide_bias(jj,rnt);
#ifdef DEBUG_SAMPLING
								if ((ii==LEFT_CHECK)&&(jj==RIGHT_CHECK)) debug += boltzmannInternal_2x1(lnt,aa,cc,uu,vv,bb,rnt) * Zs[newK][ii+3][jj-2][uu][vv].pf;
#endif
								if (va<=anchor) {
									ss_sample[0][ii+1] = index2char(ii+1,aa);
									ss_sample[1][ii+1] = '.';
									ss_sample[0][ii+2] = index2char(ii+2,cc);
									ss_sample[1][ii+2] = '.';
									ss_sample[0][jj-1] = index2char(jj-1,bb);
									ss_sample[1][jj-1] = '.';
#ifdef TRACEBACK
									fprintf(stdout ,"line %d: sample internal-loop: k=%d, i=%d, j=%d, u=%d, v=%d, va=%e, anchor=%e%s\n",__LINE__,newK,ii+3,jj-2,uu,vv,va,anchor,getJobID());
#endif
									Esample = samplingHelix(newK,ii+3,jj-2,uu,vv,ss_sample);
									return EInternal_2x1(lnt,aa,cc,uu,vv,bb,rnt) + Esample;
								}	      
							}
						}
						break;
					case 3:
						newii=ii+3;
						newjj=jj-3;
						if (!is_internal_loop(ii,newii,newjj,jj)) break; /* lookup structure constraints */
						if (((newii<cutpoint)&&(newjj>cutpoint))||((len>=__hairpin_min_length__+8)&&((newjj<cutpoint)||(newii>=cutpoint)))) { /* constrain loop to respect cutpoint */						
							xx=config>>14;
							yy=config&3;
							if ((xx==lnt)&&(yy==rnt)) {
								aa=(config>>12)&3;
								cc=(config>>10)&3;
								uu=(config>>8)&3;
								vv=(config>>6)&3;
								dd=(config>>4)&3;
								bb=(config>>2)&3;
								newK=nb_mut_remaining-nb_mut_motif + kronecker(ii+3,uu) + kronecker(jj-3,vv);
								if ((newK<0)||(newK>nb_mut_remaining)) {
									fprintf(stderr,"Incoherency: newK CANNOT be negative!! (line %d)%s\n",__LINE__,getJobID());
									return 0;
								}
								anchor += boltzmannInternal_2x2(lnt,aa,cc,uu,vv,dd,bb,rnt) * Zs[newK][ii+3][jj-3][uu][vv].pf * nucleotide_bias(ii,lnt) * nucleotide_bias(ii+1,aa) * nucleotide_bias(ii+2,cc) * nucleotide_bias(jj-2,dd) * nucleotide_bias(jj-1,bb) * nucleotide_bias(jj,rnt);
#ifdef DEBUG_SAMPLING
								if ((ii==LEFT_CHECK)&&(jj==RIGHT_CHECK)) debug += boltzmannInternal_2x2(lnt,aa,cc,uu,vv,dd,bb,rnt) * Zs[newK][ii+3][jj-3][uu][vv].pf;
#endif
								if (va<=anchor) {
									ss_sample[0][ii+1] = index2char(ii+1,aa);
									ss_sample[1][ii+1] = '.';
									ss_sample[0][ii+2] = index2char(ii+2,cc);
									ss_sample[1][ii+2] = '.';
									ss_sample[0][jj-2] = index2char(jj-2,dd);
									ss_sample[1][jj-2] = '.';
									ss_sample[0][jj-1] = index2char(jj-1,bb);
									ss_sample[1][jj-1] = '.';
#ifdef TRACEBACK
									fprintf(stdout ,"line %d: sample internal-loop: k=%d, i=%d, j=%d, u=%d, v=%d, va=%e, anchor=%e%s\n",__LINE__,newK,ii+3,jj-3,uu,vv,va,anchor,getJobID());
#endif
									Esample = samplingHelix(newK,ii+3,jj-3,uu,vv,ss_sample);
									return EInternal_2x2(lnt,aa,cc,uu,vv,dd,bb,rnt) + Esample;
								}
							}
						}
						break;
				}
				ii_list++;
			}
		}
	}
	
	
#ifdef DEBUG_SAMPLING
	if ((ii==LEFT_CHECK)&&(jj==RIGHT_CHECK)) printf("sample special: %e\n",debug-pdebug);
	pdebug = debug;
#endif
	
	/*****************************************************************************************************/
	/** 1xn and nx1 internal loops                                                                      **/
	/*****************************************************************************************************/
	
	if ((ii<cutpoint)&&(cutpoint<jj)) {
		nb_free_nt = minimum(__max_size_bulge__,len - 7);
	} else {
		nb_free_nt = minimum(__max_size_bulge__,len - __hairpin_min_length__ - 7);
	}
	
	
	/** enumerate **/
	for (nn=1;nn<=nb_free_nt;nn++) { /* number of undefined nucleotides in internal loop */
		
		nb_max_mut_close_bp=minimum(4,nb_mut_remaining);
		
		for (nb_mut_close_bp=0;nb_mut_close_bp<=nb_max_mut_close_bp;nb_mut_close_bp++) {
			
			list_config_close = stack_tab[nb_mut_close_bp][ii][jj][1];
			ii_list_close = 0;
			
			while ((config_close=list_config_close[ii_list_close])) {
				
				xx=config_close>>6;
				yy=config_close&3;
				
				if ((xx==lnt)&&(yy==rnt)) {
					aa=(config_close>>4)&3;
					bb=(config_close>>2)&3;
					
					/** 1xn **/
					
					newii=ii+2;
					newjj=jj-3-nn;
					
					if (is_internal_loop(ii,newii,newjj,jj)) { /* lookup structure constraints */
						
						if (((newii<cutpoint)&&(newjj>cutpoint))||((len>=__hairpin_min_length__+7+nn)&&((newjj<cutpoint)||(newii>=cutpoint)))) {
							
							nt2mut=nn-cst_segment[jj-1-nn][nn-1];
							nb_max_mut_unk_loop = minimum(nt2mut,nb_mut_remaining-nb_mut_close_bp); /* to be adjusted later */
							
							for (nb_mut_unk_loop=0;nb_mut_unk_loop<=nb_max_mut_unk_loop;nb_mut_unk_loop++) { /* number of mutations in unknown part of loops */
								
								/** dont forget to add kronecker because mutation on single loop is be counted twice **/
								nb_max_mut_open_bp=minimum(4,nb_mut_remaining-nb_mut_unk_loop-nb_mut_close_bp+kronecker(ii+1,aa)); 
								
								for (nb_mut_open_bp=0;nb_mut_open_bp<=nb_max_mut_open_bp;nb_mut_open_bp++) {
									
									list_config_open = stack_tab[nb_mut_open_bp][ii+1][jj-2-nn][2];
									ii_list_open = 0;
									while ((config_open=list_config_open[ii_list_open])) {
										
										cc=(config_open>>6)&3;
										
										if (aa==cc) {
											
											uu=(config_open>>4)&3;
											vv=(config_open>>2)&3;
											dd=config_open&3;
											newK = nb_mut_remaining - nb_mut_close_bp - nb_mut_open_bp - nb_mut_unk_loop + kronecker(ii+1,aa) + kronecker(ii+2,uu) + kronecker(jj-3-nn,vv);
											if ((newK<0)||(newK>nb_mut_remaining)) {
												fprintf(stderr,"Incoherency: newK CANNOT be negative!! (line %d)%s\n",__LINE__,getJobID());
												return 0;
											}
											
											/* weight = genereMutant(nt2mut,nb_mut_unk_loop) */
											weight = sequence_bias(jj-1-nn,jj-2,nb_mut_unk_loop);
											
											anchor += boltzmannInternal_generic(lnt,aa,cc,uu,vv,dd,bb,rnt,1,nn+2) *
											Zs[newK][ii+2][jj-3-nn][uu][vv].pf * weight * nucleotide_bias(ii,lnt) * nucleotide_bias(ii+1,aa) * nucleotide_bias(jj-2-nn,dd) * nucleotide_bias(jj-1,bb) * nucleotide_bias(jj,rnt);
											
#ifdef DEBUG_SAMPLING
											if ((ii==LEFT_CHECK)&&(jj==RIGHT_CHECK)) debug += boltzmannInternal_generic(lnt,aa,cc,uu,vv,dd,bb,rnt,nn+2,1) *
												Zs[newK][ii+2][jj-3-nn][uu][vv].pf * genereMutant(nt2mut,nb_mut_unk_loop);
#endif
											if (va<=anchor) {
#ifdef TRACEBACK
												fprintf(stdout ,"line %d: sample internal-loop: k=%d, i=%d, j=%d, u=%d, v=%d, va=%e, anchor=%e%s\n",__LINE__,newK,ii+2,jj-3-nn,uu,vv,va,anchor,getJobID());
#endif
												ss_sample[0][ii+1] = index2char(ii+1,aa);
												ss_sample[1][ii+1] = '.';
												ss_sample[0][jj-2-nn] = index2char(jj-2-nn,dd);
												ss_sample[1][jj-2-nn] = '.';
												ss_sample[0][jj-1] = index2char(jj-1,bb);
												ss_sample[1][jj-1] = '.';
												fill_weighted_random_mutations(nb_mut_unk_loop,jj-1-nn,jj-2,ss_sample);
												Esample = samplingHelix(newK,ii+2,jj-3-nn,uu,vv,ss_sample);
												return EInternal_generic(lnt,aa,cc,uu,vv,dd,bb,rnt,1,nn+2) + Esample;
											}
										}
										ii_list_open++;
									}
								} /* close nb_mut_open */
							} /* close nb_mut_unk_loop */
						}
					} /* close check structure constraints */
					
					/** nx1 **/
					
					newii=ii+3+nn;
					newjj=jj-2;
					
					if (is_internal_loop(ii,newii,newjj,jj)) { /* lookup structure constraints */
						
						if (((newii<cutpoint)&&(newjj>cutpoint))||((len>=__hairpin_min_length__+7+nn)&&((newjj<cutpoint)||(newii>=cutpoint)))) {
							
							nt2mut=nn-cst_segment[ii+2][nn-1];
							nb_max_mut_unk_loop = minimum(nt2mut,nb_mut_remaining-nb_mut_close_bp); /* to be adjusted later */
							
							for (nb_mut_unk_loop=0;nb_mut_unk_loop<=nb_max_mut_unk_loop;nb_mut_unk_loop++) { /* number of mutations in unknown part of loops */
								
								/** dont forget to add kronecker because mutation on single loop is be counted twice **/
								nb_max_mut_open_bp=minimum(4,nb_mut_remaining-nb_mut_unk_loop-nb_mut_close_bp+kronecker(jj-1,bb)); 
								
								for (nb_mut_open_bp=0;nb_mut_open_bp<=nb_max_mut_open_bp;nb_mut_open_bp++) {
									
									/** nx1 **/
									
									list_config_open = stack_tab[nb_mut_open_bp][ii+2+nn][jj-1][2];
									ii_list_open = 0;
									while ((config_open=list_config_open[ii_list_open])) {
										
										dd=config_open&3;
										
										if (bb==dd) {
											
											cc=(config_open>>6)&3;
											uu=(config_open>>4)&3;
											vv=(config_open>>2)&3;
											
											newK = nb_mut_remaining - nb_mut_close_bp - nb_mut_open_bp - nb_mut_unk_loop + kronecker(jj-1,bb) + kronecker(ii+3+nn,uu) + kronecker(jj-2,vv);
											if ((newK<0)||(newK>nb_mut_remaining)) {
												fprintf(stderr,"Incoherency: newK CANNOT be negative!! (line %d)%s\n",__LINE__,getJobID());
												return 0;
											}
											
											/* weight = genereMutant(nt2mut,nb_mut_unk_loop) */
											weight = sequence_bias(ii+2,ii+1+nn,nb_mut_unk_loop);
											
											anchor += boltzmannInternal_generic(lnt,aa,cc,uu,vv,dd,bb,rnt,nn+2,1) *
											Zs[newK][ii+3+nn][jj-2][uu][vv].pf * weight * nucleotide_bias(ii,lnt) * nucleotide_bias(ii+1,aa) * nucleotide_bias(ii+2+nn,cc) * nucleotide_bias(jj-1,bb) * nucleotide_bias(jj,rnt);
											
#ifdef DEBUG_SAMPLING
											if ((ii==LEFT_CHECK)&&(jj==RIGHT_CHECK)) debug += boltzmannInternal_generic(lnt,aa,cc,uu,vv,dd,bb,rnt,1,nn+2) *
												Zs[newK][ii+3+nn][jj-2][uu][vv].pf * genereMutant(nt2mut,nb_mut_unk_loop);
#endif
											if (va<=anchor) {
#ifdef TRACEBACK
												fprintf(stdout ,"line %d: sample internal-loop: k=%d, i=%d, j=%d, u=%d, v=%d, va=%e, anchor=%e%s\n",__LINE__,newK,ii+3+nn,jj-2,uu,vv,va,anchor,getJobID());
#endif
												ss_sample[0][ii+1] = index2char(ii+1,aa);
												ss_sample[1][ii+1] = '.';
												ss_sample[0][ii+2+nn] = index2char(ii+2+nn,cc);
												ss_sample[1][ii+2+nn] = '.';
												ss_sample[0][jj-1] = index2char(jj-1,bb);
												ss_sample[1][jj-1] = '.';
												fill_weighted_random_mutations(nb_mut_unk_loop,ii+2,ii+1+nn,ss_sample);
												Esample = samplingHelix(newK,ii+3+nn,jj-2,uu,vv,ss_sample);
												return EInternal_generic(lnt,aa,cc,uu,vv,dd,bb,rnt,nn+2,1) + Esample;
											}
										}
										ii_list_open++;
									}  /* close config_open */
								} /* close nb_mut_open_bp */
							} /* close nb_mut_unk_loop */
						}
					} /* close check structure constraints */
				}
				ii_list_close++;
			}
		}
	}
	
#ifdef DEBUG_SAMPLING
	if ((ii==LEFT_CHECK)&&(jj==RIGHT_CHECK)) printf("sample 1xn & nx1: %e\n",debug-pdebug);
	pdebug = debug;
#endif
	
	/*****************************************************************************************************/
	/** mxn internal loops (m+n>4)                                                                      **/
	/*****************************************************************************************************/
	
	if ((ii<cutpoint)&&(cutpoint<jj)) {
		nb_free_nt = minimum(2*__max_size_bulge__,len - 8);
	} else {
		nb_free_nt = minimum(2*__max_size_bulge__,len - __hairpin_min_length__ - 8);
	}
	
	/** lookup **/
	nb_min_mut_close_bp=kronecker(ii,lnt)+kronecker(jj,rnt);
	nb_max_mut_close_bp=minimum(4-cst_tape[ii]-cst_tape[ii+1]-cst_tape[jj-1]-cst_tape[jj],nb_mut_remaining);
	
	for (nb_mut_close_bp=nb_min_mut_close_bp;nb_mut_close_bp<=nb_max_mut_close_bp;nb_mut_close_bp++) {
		
		list_config_close = stack_tab[nb_mut_close_bp][ii][jj][1];
		ii_list_close = 0;
		
		while ((config_close=list_config_close[ii_list_close])) {
			
			xx=config_close>>6;
			aa=(config_close>>4)&3;
			bb=(config_close>>2)&3;
			yy=config_close&3;
			
			if ((xx==lnt)&&(yy==rnt)) {
				
				for (nn=1;nn<=nb_free_nt;nn++) { /* number of undefined nucleotides in internal loop */
					
					min_rr = maximum(0,nn-__max_size_bulge__);
					max_rr = minimum(nn,__max_size_bulge__);
					
					for (rr=min_rr;rr<=max_rr;rr++) {
						
						ll=nn-rr;
						newii=ii+3+ll;
						newjj=jj-3-rr;
						
						if (is_internal_loop(ii,newii,newjj,jj)) { /* lookup structure constraints */
							
							if (((newii<cutpoint)&&(newjj>cutpoint))||((len>=__hairpin_min_length__+8+nn)&&((newjj<cutpoint)||(newii>=cutpoint)))) { /* check that loop respects cutpoint */
								
								nt2mut=nn;
								if (ll>0) { nt2mut -= cst_segment[ii+2][ll-1]; }
								if (rr>0) { nt2mut -= cst_segment[jj-1-rr][rr-1]; }
								
								nb_max_mut_unk_loop = minimum(nt2mut,nb_mut_remaining-nb_mut_close_bp);
								
								for (nb_mut_unk_loop=0;nb_mut_unk_loop<=nb_max_mut_unk_loop;nb_mut_unk_loop++) { /* number of mutations in unknown part of loops */
									
									nb_max_mut_open_bp=minimum(4-cst_tape[ii+2+ll]-cst_tape[ii+3+ll]-cst_tape[jj-2-rr]-cst_tape[jj-3-rr],nb_mut_remaining-nb_mut_unk_loop-nb_mut_close_bp);
									
									for (nb_mut_open_bp=0;nb_mut_open_bp<=nb_max_mut_open_bp;nb_mut_open_bp++) {
										list_config_open = stack_tab[nb_mut_open_bp][ii+2+ll][jj-2-rr][2];
										ii_list_open = 0;
										while ((config_open=list_config_open[ii_list_open])) {
											
											cc=(config_open>>6)&3;
											uu=(config_open>>4)&3;
											vv=(config_open>>2)&3;
											dd=config_open&3;
											
											newK = nb_mut_remaining - nb_mut_close_bp - nb_mut_open_bp - nb_mut_unk_loop + kronecker(ii+3+ll,uu) + kronecker(jj-3-rr,vv);
											if ((newK<0)||(newK>nb_mut_remaining)) {
												fprintf(stderr,"Incoherency: newK CANNOT be negative!! (line %d)%s\n",__LINE__,getJobID());
												return 0;
											}
											
											/* weight = genereMutant(nt2mut,nb_mut_unk_loop) */
											weight = sequence_bias(ii+2,ii+1+ll,nb_mut_unk_loop) + sequence_bias(jj-1-rr,jj-2,nb_mut_unk_loop)
											+ sequence_bias(ii+2,jj-2,nb_mut_unk_loop)
											- sequence_bias(ii+2,jj-2-rr,nb_mut_unk_loop) - sequence_bias(ii+2+ll,jj-2,nb_mut_unk_loop)
											+ sequence_bias(ii+2+ll,jj-2-rr,nb_mut_unk_loop);
											
											anchor += boltzmannInternal_generic(lnt,aa,cc,uu,vv,dd,bb,rnt,ll+2,rr+2)
											* Zs[newK][ii+3+ll][jj-3-rr][uu][vv].pf * weight
											* nucleotide_bias(ii,lnt) * nucleotide_bias(ii+1,aa)
											* nucleotide_bias(ii+2+ll,cc) * nucleotide_bias(jj-2-rr,dd)
											* nucleotide_bias(jj-1,bb) * nucleotide_bias(jj,rnt);
											
#ifdef DEBUG_SAMPLING
											if ((ii==LEFT_CHECK)&&(jj==RIGHT_CHECK)) {
												debug += boltzmannInternal_generic(lnt,aa,cc,uu,vv,dd,bb,rnt,ll+2,rr+2)
												* Zs[newK][ii+3+ll][jj-3-rr][uu][vv].pf * weight * nucleotide_bias(ii,lnt)
												* nucleotide_bias(ii+1,aa) * nucleotide_bias(ii+2+nn,cc) * nucleotide_bias(jj-2-rr,dd)
												* nucleotide_bias(jj-1,bb) * nucleotide_bias(jj,rnt);
											}
#endif
											if (va<=anchor) {
												int mut_left, mut_right;
												
#ifdef TRACEBACK
												fprintf(stdout ,"line %d: sample internal-loop: k=%d, ll=%d, rr=%d, newi=%d, newj=%d, oldi=%d, oldj=%d, u=%d, v=%d, va=%e, anchor=%e%s\n",__LINE__,newK,ll,rr,ii+3+ll,jj-3-rr,ii,jj,uu,vv,va,anchor,getJobID());
												fprintf(stdout,"mutations open=%d, close=%d [%d,%d]\n",nb_mut_open_bp,nb_mut_close_bp,nb_min_mut_close_bp,nb_max_mut_close_bp);
#endif
												ss_sample[0][ii+1] = index2char(ii+1,aa);
												ss_sample[1][ii+1] = '.';
												ss_sample[0][ii+2+ll] = index2char(ii+2+ll,cc);
												ss_sample[1][ii+2+ll] = '.';
												ss_sample[0][jj-2-rr] = index2char(jj-2-rr,dd);
												ss_sample[1][jj-2-rr] = '.';
												ss_sample[0][jj-1] = index2char(jj-1,bb);
												ss_sample[1][jj-1] = '.';
												
												/** partitioning the number of mutations in both loops **/
												mut_left = partition_mutations(ii+2, ii+1+ll, jj-1-rr, jj-2, nb_mut_unk_loop);
												mut_right = nb_mut_unk_loop - mut_left;
#ifdef TRACEBACK
												fprintf(stdout ,"mutations in loop: %d, left: %d, right: %d.%s\n",
														nb_mut_unk_loop,mut_left,mut_right,getJobID());
#endif
												fill_weighted_random_mutations(mut_left,ii+2,ii+1+ll,ss_sample);
												fill_weighted_random_mutations(mut_right,jj-1-rr,jj-2,ss_sample);
												Esample = samplingHelix(newK,ii+3+ll,jj-3-rr,uu,vv,ss_sample);
												return EInternal_generic(lnt,aa,cc,uu,vv,dd,bb,rnt,ll+2,rr+2) + Esample;
											}
											ii_list_open++;
										}
									}
								}
								
							}
						}
					}
				}
			}
			ii_list_close++;
		}
	}
	
#ifdef DEBUG_SAMPLING
	if ((ii==LEFT_CHECK)&&(jj==RIGHT_CHECK)) printf("sample generic loop: %e\n",debug-pdebug);
	pdebug = debug;
#endif
	
	/*****************************************************************************************************/
	/* check. this line should NOT be reachable */
	/*****************************************************************************************************/
	
	fprintf(stderr,"helix sampling failed... Might be due to numerical precision (i=%d, j=%d, nb_mut_remaining=%d, va=%e, partition number=%e)%s\n",
			ii,jj,nb_mut_remaining,va,anchor,getJobID());
	fprintf(stderr,"%s\n%s\n",ss_sample[0],ss_sample[1]);
	return 0;
	
}

/***********************************************************************************************/
/* sample multi-loop                                                                           */
/***********************************************************************************************/

double samplingMultiLoop(int nb_mut_remaining,int ii,int jj, int lnt, int rnt, int last_helix, char **ss_sample) {
	
	int uu,vv,rr,ee,newK,err_bound,nb_mut_inside_boundary,nt2mut,min_gap;
	
	/* In that case, s and b are necessarily equals to 0 */
	
	double va,weight,Esample,Esample1,Esample2;
	
	/* initialization */
	
	double anchor = 0.0;
	
#ifdef TRACEBACK
	fprintf(stdout ,"MultiLoop sampling: mut=%d, ii=%d, jj=%d, lnt=%d, rnt=%d%s\n",nb_mut_remaining,ii,jj,lnt,rnt,getJobID());
#endif
	
	if (last_helix) {
		va = random_va() * Zms[nb_mut_remaining][ii][jj][lnt][rnt].pf;
	}
	else {
		va = random_va() * Zm[nb_mut_remaining][ii][jj][lnt][rnt].pf;
	}
	
	/* sample a single helix */
	
	if (((ss_cst[ii]==jj)||((ss_cst[ii]<0)&&(ss_cst[jj]<0)))&&(ii<jj)) {
		
		anchor += Zs[nb_mut_remaining][ii][jj][lnt][rnt].pf * boltzmann_penalty_close_helix(lnt,rnt);
		
		if (va<=anchor) {
#ifdef TRACEBACK
			fprintf(stdout ,"line %d: sample helix: k=%d, i=%d, j=%d, u=%d, v=%d, va=%e up2=%e, %e, %e%s\n",__LINE__,nb_mut_remaining,ii,jj,lnt,rnt,va,Zs[nb_mut_remaining][ii][jj][lnt][rnt].pf,
					Zms[nb_mut_remaining][ii][jj][lnt][rnt].pf,Zm[nb_mut_remaining][ii][jj][lnt][rnt].pf,getJobID());
#endif
			Esample = samplingHelix(nb_mut_remaining,ii,jj,lnt,rnt,ss_sample);
			return penalty_close_helix(lnt,rnt) + Esample;
		}  
	}
	
	/* try intermediate pairings in order to build multi-loop */
	
	nb_mut_inside_boundary = nb_mut_remaining - kronecker(ii,lnt);
	
	if ((ii<cutpoint)&&(cutpoint<jj)) { min_gap = 0; }
	else { min_gap = __hairpin_min_length__; }
	
	for (rr=ii+1;rr<=jj-min_gap-1;rr++) {
		for (uu=0;uu<4;uu++) {
			
			/* initiate array and extend on left-hand side */
			
			if ((last_helix)&&(single_strand_cst[rr-1]>=(rr-ii-1))&&((ss_cst[rr]==jj)||((ss_cst[rr]<0)&&(ss_cst[jj]<0)))) { /* sample the last helix */
				
				if ((!((cst_tape[ii])&&(kronecker(ii,lnt))))&&(validBasePair(uu,rnt))) {
					
					if ((jj-rr>__hairpin_min_length__)||((rr<cutpoint)&&(cutpoint<jj))) {
						if (rr-ii>1) {
							
							nt2mut = rr-ii-1 - cst_segment[ii+1][rr-ii-2];
							err_bound = minimum(nt2mut,nb_mut_inside_boundary);
							
							for (ee=0;ee<=err_bound;ee++) {
								
								newK = nb_mut_remaining - kronecker(ii,lnt) - ee;
								if ((newK<0)||(newK>nb_mut_remaining)) {
									fprintf(stderr,"Incoherency: newK CANNOT be negative!! (line %d)%s\n",__LINE__,getJobID());
									return 0;
								}
								/* weight = genereMutant(nt2mut,ee) */
								weight = sequence_bias(ii+1,rr-1,ee) * nucleotide_bias(ii,lnt);
								anchor += Zs[newK][rr][jj][uu][rnt].pf * boltzmann_penalty_close_helix(uu,rnt) * weight * boltzmannExtendMultiLoop(rr-ii);
								if (va<=anchor) {
									ss_sample[0][ii]=index2char(ii,lnt);
									ss_sample[1][ii]='.';
									fill_weighted_random_mutations(ee,ii+1,rr-1,ss_sample);
#ifdef TRACEBACK
									fprintf(stdout ,"line %d: sample helix: k=%d, i=%d, j=%d, u=%d, v=%d%s\n",__LINE__,newK,rr,jj,uu,rnt,getJobID());
#endif
									Esample = samplingHelix(newK,rr,jj,uu,rnt,ss_sample);
									return penalty_close_helix(uu,rnt) + EExtendMultiLoop(rr-ii) + Esample;
								}
							}
						}
						else { /* special case: only one residue on the left */
							newK = nb_mut_remaining - kronecker(ii,lnt);
							if ((newK<0)||(newK>nb_mut_remaining)) {
								fprintf(stderr,"Incoherency: newK CANNOT be negative!! (line %d)%s\n",__LINE__,getJobID());
								return 0;
							}
							anchor += Zs[newK][rr][jj][uu][rnt].pf * boltzmann_penalty_close_helix(uu,rnt) * boltzmannExtendMultiLoop(1) * nucleotide_bias(ii,lnt);
							if (va<=anchor) {
								ss_sample[0][ii]=index2char(ii,lnt);
								ss_sample[1][ii]='.';
#ifdef TRACEBACK
								fprintf(stdout ,"line %d: sample helix: k=%d, i=%d, j=%d, u=%d, v=%d%s\n",__LINE__,newK,rr,jj,uu,rnt,getJobID());
#endif
								Esample = samplingHelix(newK,rr,jj,uu,rnt,ss_sample);
								return penalty_close_helix(uu,rnt) + EExtendMultiLoop(1) + Esample;
							}
						}
					}
				}
			}
			else {
				
				if (((rr-ii>__hairpin_min_length__+1)||((ii<cutpoint)&&(cutpoint<rr-1)))
					&&((jj-rr>__hairpin_min_length__)||((rr<cutpoint)&&(cutpoint<jj)))) {
					if ((ss_cst[rr]==jj)||((ss_cst[rr]<0)&&(ss_cst[jj]<0))) {
						for (vv=0;vv<4;vv++) {
							if (validBasePair(vv,rnt)) {
								for (ee=0;ee<=nb_mut_remaining;ee++) {
									/* more than 2 helices */
									anchor += Zm[ee][ii][rr-1][lnt][uu].pf * Zs[nb_mut_remaining-ee][rr][jj][vv][rnt].pf * boltzmann_penalty_close_helix(vv,rnt) * boltzmannAddStemInMultiLoop(1);
									if (va<=anchor) {
#ifdef TRACEBACK
										fprintf(stdout ,"line %d: sample multi-loop: k=%d, i=%d, j=%d, u=%d, v=%d%s\n",__LINE__,ee,ii,rr-1,lnt,uu,getJobID());
										fprintf(stdout ,"AND sample helix: k=%d, i=%d, j=%d, u=%d, v=%d%s\n",nb_mut_remaining-ee,rr,jj,vv,rnt,getJobID());
#endif
										Esample1 = samplingMultiLoop(ee,ii,rr-1,lnt,uu,0,ss_sample);
										Esample2 = samplingHelix(nb_mut_remaining-ee,rr,jj,vv,rnt,ss_sample);
										return penalty_close_helix(vv,rnt) + EAddStemInMultiLoop(1) + Esample1 + Esample2;
									}
									/* exactly 2 helices */
									anchor += Zms[ee][ii][rr-1][lnt][uu].pf * Zs[nb_mut_remaining-ee][rr][jj][vv][rnt].pf * boltzmann_penalty_close_helix(vv,rnt) * boltzmannAddStemInMultiLoop(2);
									if (va<=anchor) {
#ifdef TRACEBACK
										fprintf(stdout ,"line %d: sample last multi-loop: k=%d, i=%d, j=%d, u=%d, v=%d%s\n",__LINE__,ee,ii,rr-1,lnt,uu,getJobID());
										fprintf(stdout ,"AND sample helix: k=%d, i=%d, j=%d, u=%d, v=%d%s\n",nb_mut_remaining-ee,rr,jj,vv,rnt,getJobID());
#endif
										Esample1 = samplingMultiLoop(ee,ii,rr-1,lnt,uu,1,ss_sample);
										Esample2 = samplingHelix(nb_mut_remaining-ee,rr,jj,vv,rnt,ss_sample);
										return penalty_close_helix(vv,rnt) + EAddStemInMultiLoop(2) + Esample1 + Esample2;
									}
								}
							}
						}
					}
				}
			}
		}
	}
	
	/* extend all tables with an unpaired nucleotide on right side */
	
	if ((ii<cutpoint)&&(cutpoint<jj)) { min_gap = 0; }
	else { min_gap = __hairpin_min_length__; }
	
	if ((single_strand_cst[jj]>0)&&(min_gap+2)) { /* lookup structure constraints */	
		if (!((cst_tape[jj])&&(kronecker(jj,rnt)))) {
			for (uu=0;uu<4;uu++) {
				/* right */
				newK = nb_mut_remaining - kronecker(jj,rnt);
				if ((newK<0)||(newK>nb_mut_remaining)) {
					fprintf(stderr,"Incoherency: newK CANNOT be negative!! (line %d)%s\n",__LINE__,getJobID());
					return 0;
				}
				anchor += boltzmannExtendMultiLoop(1) * Zms[newK][ii][jj-1][lnt][uu].pf * nucleotide_bias(jj,rnt);
				if (va<=anchor) {
					ss_sample[0][jj]=index2char(jj,rnt);
					ss_sample[1][jj]='.';
#ifdef TRACEBACK
					fprintf(stdout ,"line %d: sample multi-loop: k=%d, i=%d, j=%d, u=%d, v=%d%s\n",__LINE__,newK,ii,jj-1,lnt,uu,getJobID());
#endif
					Esample = samplingMultiLoop(newK,ii,jj-1,lnt,uu,last_helix,ss_sample);
					return EExtendMultiLoop(1) + Esample;
				}
				anchor += boltzmannExtendMultiLoop(1) * Zm[newK][ii][jj-1][lnt][uu].pf * nucleotide_bias(jj,rnt);
				if (va<=anchor) {
					ss_sample[0][jj]=index2char(jj,rnt);
					ss_sample[1][jj]='.';
#ifdef TRACEBACK
					fprintf(stdout ,"line %d: sample multi-loop: k=%d, i=%d, j=%d, u=%d, v=%d%s\n",__LINE__,newK,ii,jj-1,lnt,uu,getJobID());
#endif
					Esample = samplingMultiLoop(newK,ii,jj-1,lnt,uu,last_helix,ss_sample);
					return EExtendMultiLoop(1) + Esample;
				}
			}
		}
	}
	
	/* check. this line should NOT be reachable */
	
	fprintf(stderr,"Multi-loop sampling failed... Might be due to numerical precision (i=%d, j=%d, nb_mut_remaining=%d, va=%f, partition number=%f)%s\n",
			ii,jj,nb_mut_remaining,va,anchor,getJobID());
	
	return 0;
	
}  

/********************************/
/**** Sampling Exterior Loop ****/
/********************************/

double samplingExteriorLoop(int nb_mut_remaining,int ii,int jj, int lnt, int rnt, char **ss_sample) {
	
	int uu,vv,rr,ee,newK,err_bound,len=jj-ii+1,nb_mut_inside_boundary,nt2mut,min_gap;
	
	/* In that case, s and b are necessarily equals to 0 */
	
	double va, weight, Esample, Esample1, Esample2;
	
	/* initialization */
	
	double anchor = 0.0;
	
#ifdef TRACEBACK
	fprintf(stdout ,"Exterior helix sampling: mut=%d, ii=%d, jj=%d, lnt=%d, rnt=%d%s\n",nb_mut_remaining,ii,jj,lnt,rnt,getJobID());
#endif
	
	va = random_va() * (Zes[nb_mut_remaining][ii][jj][lnt][rnt].pf + Ze[nb_mut_remaining][ii][jj][lnt][rnt].pf);
	
	/* sample a single helix */
	
	if ((ii<jj)&&((ss_cst[ii]==jj)||((ss_cst[ii]<0)&&(ss_cst[jj]<0)))) { 
		
		anchor += Zs[nb_mut_remaining][ii][jj][lnt][rnt].pf * boltzmann_penalty_close_helix(lnt,rnt);
		
		if (va<=anchor) {
#ifdef TRACEBACK
			fprintf(stdout ,"line %d: sample helix: k=%d, i=%d, j=%d, u=%d, v=%d%s\n",__LINE__,nb_mut_remaining,ii,jj,lnt,rnt,getJobID());
#endif
			Esample = samplingHelix(nb_mut_remaining,ii,jj,lnt,rnt,ss_sample);
			return penalty_close_helix(lnt,rnt) + Esample;
		}  
	}
	
	/* try intermediate pairings in order to build multi-loop */
	
	nb_mut_inside_boundary = nb_mut_remaining - kronecker(ii,lnt);
	
	if ((ii<cutpoint)&&(cutpoint<jj)) { min_gap = 0; }
	else { min_gap = __hairpin_min_length__; }
	
	for (rr=ii+1;rr<=jj-min_gap-1;rr++) {
		for (uu=0;uu<4;uu++) {
			
			/* initiate array and extend on left-hand side */
			
			if ((!((cst_tape[ii])&&(kronecker(ii,lnt))))&&(validBasePair(uu,rnt))) {
				
				if ((single_strand_cst[rr-1]>=(rr-ii-1))&&((ss_cst[rr]==jj)||((ss_cst[rr]<0)&&(ss_cst[jj]<0)))) { /* sample the last helix and thus check appropriate constraints */
					
					if ((!((cst_tape[ii])&&(kronecker(ii,lnt))))&&(validBasePair(uu,rnt))) {
						
						if ((jj-rr>__hairpin_min_length__)||((rr<cutpoint)&&(cutpoint<jj))) {
							if ((rr-ii>1)&&(single_strand_cst[rr-1]>=(rr-ii-1))&&((ss_cst[rr]==jj)||((ss_cst[rr]<0)&&(ss_cst[jj]<0)))) {
								
								nt2mut = rr-ii-1 - cst_segment[ii+1][rr-ii-2];
								err_bound = minimum(nt2mut,nb_mut_inside_boundary);
								
								for (ee=0;ee<=err_bound;ee++) {
									
									newK = nb_mut_remaining - kronecker(ii,lnt) - ee;
									/* weight = genereMutant(nt2mut,ee) */
									weight = sequence_bias(ii+1,rr-1,ee) * nucleotide_bias(ii,lnt);
									anchor += Zs[newK][rr][jj][uu][rnt].pf * boltzmann_penalty_close_helix(uu,rnt) * weight;
									if (va<=anchor) {
										ss_sample[0][ii]=index2char(ii,lnt);
										ss_sample[1][ii]='.';
										fill_weighted_random_mutations(ee,ii+1, rr-1, ss_sample);
#ifdef TRACEBACK
										fprintf(stdout ,"line %d: sample helix: k=%d, i=%d, j=%d, u=%d, v=%d%s\n",__LINE__,newK,rr,jj,uu,rnt,getJobID());
#endif
										Esample = samplingHelix(newK,rr,jj,uu,rnt,ss_sample);
										return penalty_close_helix(uu,rnt) + Esample;
									}
								}
							}
							else { /* special case: only one residue on the left */
								newK = nb_mut_remaining - kronecker(ii,lnt);
								if ((newK<0)||(newK>nb_mut_remaining)) {
									fprintf(stderr,"Incoherency: newK CANNOT be negative!! (line %d)%s\n",__LINE__,getJobID());
									return 0;
								}
								anchor += Zs[newK][rr][jj][uu][rnt].pf * boltzmann_penalty_close_helix(uu,rnt) * nucleotide_bias(ii,lnt);
								if (va<=anchor) {
									ss_sample[0][ii]=index2char(ii,lnt);
									ss_sample[1][ii]='.';
#ifdef TRACEBACK
									fprintf(stdout ,"line %d: sample helix: k=%d, i=%d, j=%d, u=%d, v=%d%s\n",__LINE__,newK,rr,jj,uu,rnt,getJobID());
#endif
									Esample = samplingHelix(newK,rr,jj,uu,rnt,ss_sample);
									return penalty_close_helix(uu,rnt) + Esample;
								}
							}
						}
					}
				}
			}
			
			if (((rr-ii>__hairpin_min_length__+1)||((ii<cutpoint)&&(cutpoint<rr-1)))
				&&((jj-rr>__hairpin_min_length__)||((rr<cutpoint)&&(cutpoint<jj)))) {
				if ((ss_cst[rr]==jj)||((ss_cst[rr]<0)&&(ss_cst[jj]<0))) {
					for (vv=0;vv<4;vv++) {
						if (validBasePair(vv,rnt)) {
							for (ee=0;ee<=nb_mut_remaining;ee++) {
								/* more than 2 helices */
								anchor += Ze[ee][ii][rr-1][lnt][uu].pf * Zs[nb_mut_remaining-ee][rr][jj][vv][rnt].pf * boltzmann_penalty_close_helix(vv,rnt);
								if (va<=anchor) {
#ifdef TRACEBACK
									fprintf(stdout ,"line %d: sample exterior loop: k=%d, i=%d, j=%d, u=%d, v=%d%s\n",__LINE__,ee,ii,rr-1,lnt,uu,getJobID());
									fprintf(stdout ,"AND sample helix: k=%d, i=%d, j=%d, u=%d, v=%d%s\n",nb_mut_remaining-ee,rr,jj,vv,rnt,getJobID());
#endif
									Esample1 = samplingExteriorLoop(ee,ii,rr-1,lnt,uu,ss_sample);
									Esample2 = samplingHelix(nb_mut_remaining-ee,rr,jj,vv,rnt,ss_sample);
									return penalty_close_helix(vv,rnt) + Esample1 + Esample2;
								}
								/* exactly 2 helices */
								anchor += Zes[ee][ii][rr-1][lnt][uu].pf * Zs[nb_mut_remaining-ee][rr][jj][vv][rnt].pf * boltzmann_penalty_close_helix(vv,rnt);
								if (va<=anchor) {
#ifdef TRACEBACK
									fprintf(stdout ,"line %d: sample exterior loop: k=%d, i=%d, j=%d, u=%d, v=%d%s\n",__LINE__,ee,ii,rr-1,lnt,uu,getJobID());
									fprintf(stdout ,"AND sample helix: k=%d, i=%d, j=%d, u=%d, v=%d%s\n",nb_mut_remaining-ee,rr,jj,vv,rnt,getJobID());
#endif
									Esample1 = samplingExteriorLoop(ee,ii,rr-1,lnt,uu,ss_sample);
									Esample2 = samplingHelix(nb_mut_remaining-ee,rr,jj,vv,rnt,ss_sample);
									return penalty_close_helix(vv,rnt) + Esample1 + Esample2;
								}
							}
						}
					}
				}
			}
		}
	}
	
	/* extend all tables with an unpaired nucleotide on right side */
	
	if ((ii<cutpoint)&&(cutpoint<jj)) { min_gap = 0; }
	else { min_gap = __hairpin_min_length__; }
	
	if ((single_strand_cst[jj]>0)&&(len>min_gap+2)) { /* minimum length required */
		if (!((cst_tape[jj])&&(kronecker(jj,rnt)))) {
			for (uu=0;uu<4;uu++) {
				/* right */
				newK = nb_mut_remaining - kronecker(jj,rnt);
				if ((newK<0)||(newK>nb_mut_remaining)) {
					fprintf(stderr,"ERROR:%s:%d: newK should not be negative (newK=%d,nb_mut_remaining=%d).%s\n",__FILE__,__LINE__,newK,nb_mut_remaining,getJobID());
					return 0;
				}
				anchor += Zes[newK][ii][jj-1][lnt][uu].pf * nucleotide_bias(jj,rnt);
				if (va<=anchor) {
					ss_sample[0][jj]=index2char(jj,rnt);
					ss_sample[1][jj]='.';
#ifdef TRACEBACK
					fprintf(stdout ,"line %d: sample exterior loop: k=%d, i=%d, j=%d, u=%d, v=%d%s\n",__LINE__,newK,ii,jj-1,lnt,uu,getJobID());
#endif
					Esample = samplingExteriorLoop(newK,ii,jj-1,lnt,uu,ss_sample);
					return Esample;
				}
				anchor += Ze[newK][ii][jj-1][lnt][uu].pf * nucleotide_bias(jj,rnt);
				if (va<=anchor) {
					ss_sample[0][jj]=index2char(jj,rnt);
					ss_sample[1][jj]='.';
#ifdef TRACEBACK
					fprintf(stdout ,"line %d: sample exterior loop: k=%d, i=%d, j=%d, u=%d, v=%d%s\n",__LINE__,newK,ii,jj-1,lnt,uu,getJobID());
#endif
					Esample = samplingExteriorLoop(newK,ii,jj-1,lnt,uu,ss_sample);
					return Esample;
				}
			}
		}
	}
	
	/* check. this line should NOT be reachable */
	
	fprintf(stderr,"Exterior loop sampling failed... Might be due to numerical precision (i=%d, j=%d, nb_mut_remaining=%d, va=%e, partition number=%e)%s\n",
			ii,jj,nb_mut_remaining,va,anchor,getJobID());
	
	return 0;
	
}  

/****************************************************************************************************/
/* start sampling a sequence (i,j)                                                                  */
/****************************************************************************************************/

int startBasicSampling(int i, int j, char **ss_sample, double *Esample) {
	
	int k,xx,yy;
	double va;
	
	/* initialization */
	
	double anchor = 0.0;
	double partitionFunction = 0.0; /* init with empty structure */
	
	for (k=0;k<max_mutations;k++) {
		if (!ss_constraint_is_non_empty()) { /* check constraints */
			/* partitionFunction += genereMutant(j-i+1,k); */
			partitionFunction += sequence_bias(i,j,k);
		}
		for (xx=0;xx<4;xx++) {
			for (yy=0;yy<4;yy++) {
#if 0
				printf("(%d,%d,%d,%d,%d)---> %e, %e\n",k,i,j,xx,yy,Zes[k][i][j][xx][yy].mfe, Ze[k][i][j][xx][yy].mfe);
				printf("(%d,%d,%d,%d,%d)---> %e, %e\n",k,i,j,xx,yy,Zes[k][i][j][xx][yy].pf, Ze[k][i][j][xx][yy].pf);
#endif
				if (Zes[k][i][j][xx][yy].pf>DBL_EPSILON) {
					partitionFunction += Zes[k][i][j][xx][yy].pf;
				}
				if (Ze[k][i][j][xx][yy].pf>DBL_EPSILON) {
					partitionFunction += Ze[k][i][j][xx][yy].pf;
				}
			}
		}
	}
	
	/* va instanciation */
	
	va = random_va() * partitionFunction; 
	
	/* determine the config to sample as well as the number of mutations */
	
	for (k=0;k<max_mutations;k++) {
		
		/* empty structure */
		if (!ss_constraint_is_non_empty()) { /* check constraints */
			/* anchor += genereMutant(len,k); */
			anchor += sequence_bias(i,j,k);
			if (va<=anchor) {
				/* start to sample an exterior loop */
				fill_weighted_random_mutations(k,i, j, ss_sample);
				*Esample = 0;
				return k;
			}
		}
		
		/* general case */
		
		for (xx=0;xx<4;xx++) {
			for (yy=0;yy<4;yy++) {
				
				if (Zes[k][i][j][xx][yy].pf>DBL_EPSILON) {
					anchor += Zes[k][i][j][xx][yy].pf;
				}
				if (Ze[k][i][j][xx][yy].pf>DBL_EPSILON) {
					anchor += Ze[k][i][j][xx][yy].pf;
				}
				
				if (va<=anchor) {
					/* start to sample an exterior loop */
					*Esample = samplingExteriorLoop(k,i,j,xx,yy,ss_sample);
					return k;
				}
			}
		}
	}
	
	fprintf(stderr,"Basic sampling start failed... Anchor value of %e is found instead of %e.%s\n",anchor,partitionFunction,getJobID());
	
	return -1;
	
}

/****************************************************************************************************/
/* start sampling a sequence (i,j)                                                                  */
/****************************************************************************************************/

void startSamplingKmutant(int k, int i, int j, char **ss_sample, double *Esample) {
	
	int xx,yy;
	double va;
	
	/* initialization */
	
	double anchor = 0.0;
	double partitionFunction = 0.0; /* init with empty structure */
	
	if (!ss_constraint_is_non_empty()) { /* check constraints */
		/* partitionFunction += genereMutant(rna_len,k); */
		partitionFunction += sequence_bias(0,rna_len-1,k);
	}
	
	for (xx=0;xx<4;xx++) {
		for (yy=0;yy<4;yy++) {
#if 0
			printf("(%d,%d,%d,%d,%d)---> %e, %e\n",k,i,j,xx,yy,Zes[k][i][j][xx][yy] + Ze[k][i][j][xx][yy]);
#endif
			if (Zes[k][i][j][xx][yy].pf>DBL_EPSILON) {
				partitionFunction += Zes[k][i][j][xx][yy].pf;
			}
			if (Ze[k][i][j][xx][yy].pf>DBL_EPSILON) {
				partitionFunction += Ze[k][i][j][xx][yy].pf;
			}
		}
	}
	
	
	/* va instanciation */
	
	va = random_va() * partitionFunction; 
	
	/* determine the config to sample as well as the number of mutations */
	
	/* empty structure */
	
	if (!ss_constraint_is_non_empty()) { /* check constraints */
		/* anchor += genereMutant(rna_len,k); */
		anchor += sequence_bias(0,rna_len-1,k);
		if (va<=anchor) {
			/* start to sample an exterior loop */
			fill_weighted_random_mutations(k,i, j, ss_sample);
			*Esample = 0.0;
			return;
		}
	}
	
	/* general case */
	
	for (xx=0;xx<4;xx++) {
		for (yy=0;yy<4;yy++) {
			
			if (Zes[k][i][j][xx][yy].pf>DBL_EPSILON) {
				anchor += Zes[k][i][j][xx][yy].pf;
			}
			if (Ze[k][i][j][xx][yy].pf>DBL_EPSILON) {
				anchor += Ze[k][i][j][xx][yy].pf;
			}
			
			if (va<=anchor) {
				/* start to sample an exterior loop */
				*Esample = samplingExteriorLoop(k,i,j,xx,yy,ss_sample);
				return;
			}
		}
	}
	
	fprintf(stderr,"Sampling start failed... Anchor value of %e is found instead of %e.%s\n",anchor,partitionFunction,getJobID());
	
}

/************************************************************************************************************/
/*                                                                                                          */
/* Sample k structures                                                                                      */
/*                                                                                                          */
/************************************************************************************************************/


int basicSamplingEngine(FILE *sample_file_flux, int nos, int stat_flag, int warning_flag, int compatible_neighbors) {
	
	int nbSample = 0, ii,jj;
	char **ss_sample;
	char *buffer_seq1=NULL,*buffer_ss1=NULL;
	int length_seq1=0;
	int nb_mut_in_sample,nb_corrupted=0,error_rate=0;
	double Esample;
	
	/* init array if cutpoint is used */
	
	if (cutpoint) {
		length_seq1=(int)(cutpoint+0.5);
		buffer_seq1=(char *)xmalloc((length_seq1+1)*sizeof(char));
		buffer_ss1=(char *)xmalloc((length_seq1+1)*sizeof(char));
	}
	
	/* initialization */
	
	ss_sample = (char **)xmalloc(2*sizeof(char *));
	ss_sample[0] = emptyRNAss();
	ss_sample[1] = emptyRNAss();
	
	if (stat_flag) {
		
		/* initialize base pair count table */
		
		bp_stat_table=(long int **)xmalloc(rna_len*sizeof(long int *));
		for (ii=0;ii<rna_len;ii++) {
			bp_stat_table[ii]=(long int *)xmalloc(rna_len*sizeof(long int));
			memset(bp_stat_table[ii],0,rna_len*sizeof(long int));
		}
		
		/* initialize mutation count table */
		
		mutation_stat_table=(long int ***)xmalloc(rna_len*sizeof(long int **));
		for (ii=0;ii<rna_len;ii++) {
			mutation_stat_table[ii]=(long int **)xmalloc(4*sizeof(long int *));  /* fields for nucleotides */
			for (jj=0;jj<4;jj++) {
				mutation_stat_table[ii][jj]=(long int *)xmalloc(rna_len*sizeof(long int));
				memset(mutation_stat_table[ii][jj],0,rna_len*sizeof(long int));
			}
		}
		
	}
	
	/* check */
	
	if ((ss_constraint_is_non_empty())&&(!compatible_neighbors)) {
		fprintf(stderr,"> Cannot sample %d-mutants. No sequence meets the structural constraints.\n",nos);
		return 0;
	}
	
	/* sampling */
	
	fprintf(sample_file_flux,"> sampling %d mutant(s) and secondary structure(s)\n",nos);
	error_rate = (int)((double)(nos)/1000.0) + 1;
	ii=0;
	while (ii<nos) {
		cleanRNAss(ss_sample[0]);
		cleanRNAss(ss_sample[1]);
		nb_mut_in_sample = startBasicSampling(0,rna_len-1,ss_sample,&Esample);
		if (check_sample(ss_sample[0],nb_mut_in_sample)) {
			if (stat_flag) {
				fill_stat_tables(ss_sample[0],ss_sample[1]);
			}
			else {
				if (cutpoint) {
					strncpy(buffer_seq1,ss_sample[0],(length_seq1)*sizeof(char));
					strncpy(buffer_ss1,ss_sample[1],(length_seq1)*sizeof(char));
					buffer_seq1[length_seq1]='\0';
					buffer_ss1[length_seq1]='\0';
					fprintf(sample_file_flux,"%s&%s\t(%.2f)\n",buffer_seq1,ss_sample[0]+length_seq1,Esample);
					fprintf(sample_file_flux,"%s&%s\n",buffer_ss1,ss_sample[1]+length_seq1);
				}
				else {
					fprintf(sample_file_flux,"%s\t(%.2f)\n",ss_sample[0],Esample);
					fprintf(sample_file_flux,"%s\n",ss_sample[1]);
				}
			}
			ii++; /* increment only if sample is valid */
		}
		else {
			if (nb_corrupted>error_rate) {
				fprintf(stderr,"%s:%d:ERROR: Sampling engine failed while sampling %d-mutants. Please report bug.%s\n",__FILE__,__LINE__,nb_mut_in_sample,getJobID());
				fprintf(stderr,"%s\t(%.2f)\n",ss_sample[0],Esample);
				fprintf(stderr,"%s\n",ss_sample[1]);
				fprintf(stderr,"exit\n");
				return ii;
			}
			nb_corrupted++;
		}
		nbSample++;
	} 
	
	/* print statistics*/
	
	if (stat_flag) {
		print_stat_tables(max_mutations,nos);
	}
	
	/* clean tables */
	
	if (stat_flag) {
		
		/* initialize base pair count table */
		
		for (ii=0;ii<rna_len;ii++) {
			free(bp_stat_table[ii]);
		}
		free(bp_stat_table);
		
		/* initialize mutation count table */
		
		for (ii=0;ii<rna_len;ii++) {
			for (jj=0;jj<4;jj++) {
				free(mutation_stat_table[ii][jj]);
			}
			free(mutation_stat_table[ii]);
		}
		free(mutation_stat_table);
		
	}
	
	free(ss_sample[1]);
	free(ss_sample[0]);
	free(ss_sample);	
	
	if (cutpoint) {
		free(buffer_seq1);
		free(buffer_ss1);
	}
	
	/* print collisions */
	
	if (warning_flag) {
		if (nb_corrupted>0) {
			fprintf(stderr,"WARNING: %d corrupted samples generated and re-sampled.%s\n",nb_corrupted,getJobID());
		}
	}
	/* exit */
	
	return nbSample;
}

/************************************************************************************************************/
/*                                                                                                          */
/* Read a command file storing all sampling request, and execute it                                         */
/*                                                                                                          */
/************************************************************************************************************/

int sampleFromFile(FILE *sample_file_flux, const char *commandfile, int stat_flag, int warning_flag, int *compatible_neighbors) {
	int nbSample = 0, nb_mut_in_sample, howManySample, ii, jj, line=1;
	FILE *commands;
	char **ss_sample,buffer[1000];
	char *buffer_seq1=NULL,*buffer_ss1=NULL;
	int length_seq1=0, nb_corrupted[1000], indcorr=0, error_rate=0;
	double Esample;
	
	/* init array if cutpoint is used */
	
	if (cutpoint) {
		length_seq1=(int)(cutpoint+0.5);
		buffer_seq1=(char *)xmalloc((length_seq1+1)*sizeof(char));
		buffer_ss1=(char *)xmalloc((length_seq1+1)*sizeof(char));
	}
	
	/* open file */
	
	if ((commands=fopen(commandfile,"r"))==NULL) {
		fprintf(stderr,"Cannot open command file %s%s\n",commandfile,getJobID());
		return 0;
	}
	
	/* initialization */
	
	ss_sample = (char **)xmalloc(2*sizeof(char *));
	ss_sample[0] = emptyRNAss();
	ss_sample[1] = emptyRNAss();
	
	if (stat_flag) {
		
		/* initialize base pair count table */
		
		bp_stat_table=(long int **)xmalloc(rna_len*sizeof(long int *));
		for (ii=0;ii<rna_len;ii++) {
			bp_stat_table[ii]=(long int *)xmalloc(rna_len*sizeof(long int));
		}
		
		/* initialize mutation count table */
		
		mutation_stat_table=(long int ***)xmalloc(rna_len*sizeof(long int **));
		for (ii=0;ii<rna_len;ii++) {
			mutation_stat_table[ii]=(long int **)xmalloc(4*sizeof(long int *));  /* fields for nucleotides */
			for (jj=0;jj<4;jj++) {
				mutation_stat_table[ii][jj]=(long int *)xmalloc(rna_len*sizeof(long int));
			}
		}
		
	}
	
	/* sampling */
	
	while (fgets(&buffer[0],1000,commands) != NULL) {
		
		if ((sscanf(buffer,"%d%d",&nb_mut_in_sample,&howManySample)) != 2) {
			fprintf(stderr,"Corrupted sampling command file %s.%s\n",commandfile,getJobID());
		}
		
		/* checks */ 
		
		if ((ss_constraint_is_non_empty())&&(!compatible_neighbors[nb_mut_in_sample])) {
			fprintf(stderr,"> Cannot sample %d-mutants. No sequence meets the structural constraints.\n",nb_mut_in_sample);
			continue;
		}
		
		if ((nb_mut_in_sample>=max_mutations)||(nb_mut_in_sample<0)) {
			fprintf(stderr,"%s\n: line %d: Cannot sample mutants with more than %d mutations. Re-run with \"-m %d\".%s\n",commandfile,line,max_mutations,nb_mut_in_sample,getJobID());
			continue;
		}
		
		/* init values */ 	  
		
		indcorr++;
		nb_corrupted[indcorr]=nb_mut_in_sample;
		indcorr++;
		nb_corrupted[indcorr]=0;
		
		/* reset tables */
		
		if (stat_flag) {
			
			/* initialize base pair count table */
			
			for (ii=0;ii<rna_len;ii++) {
				memset(bp_stat_table[ii],0,rna_len*sizeof(long int));
			}
			
			/* initialize mutation count table */
			
			for (ii=0;ii<rna_len;ii++) {
				for (jj=0;jj<4;jj++) {
					memset(mutation_stat_table[ii][jj],0,rna_len*sizeof(long int));
				}
			}
			
		}
		
		/* sample <howManySample> sequence and structures with <nb_mut_in_sample> mutations */
		
		fprintf(sample_file_flux,"> sampling %d sequence and secondary structure(s) with %d mutations\n",howManySample,nb_mut_in_sample);
		
		error_rate = (int)((double)(howManySample)/1000.0) + 1;
		ii=0;
		while (ii<howManySample) {
			cleanRNAss(ss_sample[0]);
			cleanRNAss(ss_sample[1]);
			startSamplingKmutant(nb_mut_in_sample,0,rna_len-1,ss_sample,&Esample);
			if (check_sample(ss_sample[0],nb_mut_in_sample)) {
				if (stat_flag) {
					fill_stat_tables(ss_sample[0],ss_sample[1]);
				}
				else {
					if (cutpoint) {
						strncpy(buffer_seq1,ss_sample[0],(length_seq1)*sizeof(char));
						strncpy(buffer_ss1,ss_sample[1],(length_seq1)*sizeof(char));
						buffer_seq1[length_seq1]='\0';
						buffer_ss1[length_seq1]='\0';
						fprintf(sample_file_flux,"%s&%s\t(%.2f)\n",buffer_seq1,ss_sample[0]+length_seq1,Esample);
						fprintf(sample_file_flux,"%s&%s\n",buffer_ss1,ss_sample[1]+length_seq1);
					}
					else {
						fprintf(sample_file_flux,"%s\t(%.2f)\n",ss_sample[0],Esample);
						fprintf(sample_file_flux,"%s\n",ss_sample[1]);
					}
				}
				ii++; /* increment only if sample is valid */
			}
			else {
				if (nb_corrupted[indcorr]>error_rate) {
					fprintf(stderr,"%s:%d:ERROR: Sampling engine failed while sampling %d-mutants. Please report bug.%s\n",__FILE__,__LINE__,nb_mut_in_sample,getJobID());
					fprintf(stderr,"%s\t(%.2f)\n",ss_sample[0],Esample);
					fprintf(stderr,"%s\n",ss_sample[1]);
					fprintf(stderr,"exit\n");
					return ii;
				}
				nb_corrupted[indcorr]++;					
			}
		}
		
		/* print statistics*/
		
		if (stat_flag) {
			print_stat_tables(nb_mut_in_sample,howManySample);
		}
		
		nbSample += howManySample;
		line ++;
		
	}
	
	/* clean tables */
	
	if (stat_flag) {
		
		/* initialize base pair count table */
		
		for (ii=0;ii<rna_len;ii++) {
			free(bp_stat_table[ii]);
		}
		free(bp_stat_table);
		
		/* initialize mutation count table */
		
		for (ii=0;ii<rna_len;ii++) {
			for (jj=0;jj<4;jj++) {
				free(mutation_stat_table[ii][jj]);
			}
			free(mutation_stat_table[ii]);
		}
		free(mutation_stat_table);
		
	}
	
	free(ss_sample[1]);
	free(ss_sample[0]);
	free(ss_sample);
	
	if (cutpoint) {
		free(buffer_seq1);
		free(buffer_ss1);
	}
	
	/* print collisions */
	
	if (warning_flag) {
		for (ii=0;ii<=indcorr;ii+=2) {
			if (nb_corrupted[ii+1]>0) {
				fprintf(stderr,"WARNING: %d corrupted samples with %d mutations generated and re-sampled.%s\n",nb_corrupted[ii+1],nb_corrupted[ii],getJobID());
			}
		}
	}
	
	/* exit */
	
	fclose(commands);
	
	return nbSample;
}
