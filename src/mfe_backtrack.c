#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <string.h>
#include "reader_energy_params.h"
#include "energy_functions.h"
#include "util.h"
#include "mfe_backtrack.h"
#include "energy_params_tables.h"
#include "sampling.h"
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
extern int nb_value_triloop_table;
extern double *triloop_cmpt_table;
extern int nb_value_tetraloop_table;
extern double *tetraloop_cmpt_table;
extern int uncst_mutation[2][4][4];

#define NO_TRACEBACK
#define NO_DEBUG_SAMPLING
#define NO_DEBUG_MFE
#define LEFT_CHECK 7
#define RIGHT_CHECK 26


const double __backtrack_double_precision__ = 0.01;

inline int are_equal(double a, double b) {
	
	if (fabs(b-a)<__backtrack_double_precision__ ) {
		return 1;
	}
	return 0; 
}

/**** variables for hybrid ****/

extern double mfe_hybrid_initialization_energy;
extern double pf_hybrid_initialization_energy;
extern int include_intermolecular_interactions;
extern double cutpoint;

/***********************************************************************************************/
/* sampling random mutation in a segment                                                       */
/***********************************************************************************************/


void mfe_fill_random_mutations(int nb_mutations, int i, int j, char **ss_sample) {
	int nt, ilist, *nt_list;
	double va,anchor;
	
#ifdef TRACEBACK
	fprintf(stderr,"fill mfe random: mut=%d, ii=%d, jj=%d: %s\n",nb_mutations,i,j,getJobID());
	fflush(stderr);
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
			ss_sample[0][i]=input_tape[i];
			ss_sample[1][i]='.';
			mfe_fill_random_mutations(nb_mutations, i+1, j, ss_sample);
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
				ss_sample[0][i]=index2char(i,nt);
				ss_sample[1][i]='.';
				mfe_fill_random_mutations(nb_mutations-1, i+1, j, ss_sample);
				return;			
			}
			ilist++;
		}
	}
	
	fprintf(stderr,"%s:line %d: sampling failed (i=%d,j=%d,k=%d,va=%e,anchor=%e,max=%f).%s\n",__FILE__,__LINE__,j,j,nb_mutations,va,anchor,sequence_bias(i,j,nb_mutations),getJobID());
	
#ifdef TRACEBACK
	fprintf(stderr,"done%s\n",getJobID());
#endif
	
}


/***********************************************************************************************/
/* sampling hairpins                                                                           */
/*                                                                                             */
/* WARNING: We assume that there is no overlap between triloop, tetraloop and GGG loop bonuses */
/*                                                                                             */
/***********************************************************************************************/

void backtrack_triloop(int nb_mut_in_unk,int ii,int jj, int lnt, int uu, int vv, int rnt, char **ss_sample, double bonus_value) {
	int ilist1,base_hcode,hcode;
	int *nt_list1,nt1;
	int ind_table,final_rank,current_rank;
	
	triloop_counter(nb_mut_in_unk,ii,jj,lnt,uu,vv,rnt);
	ind_table=0;
	while (triloop_cmpt_table[2*ind_table]!=bonus_value) {
		ind_table++;
		if (ind_table>=nb_value_triloop_table) {
			fprintf(stderr,"%s:line %d: Cannot find bonus value %e into triloop table.%s\n",__FILE__,__LINE__,bonus_value,getJobID());
			exit(EXIT_FAILURE);
		}
	}
	
	/* rank of sample triloop */
	final_rank = (int)(random_va() * triloop_cmpt_table[2*ind_table+1]);
	
	base_hcode = (lnt<<8) | (uu<<6) | (vv<<2) | rnt;
	ilist1=0;
	nt_list1 = uncst_mutation[nb_mut_in_unk][char2index(ii+2)];
	current_rank = 0;
	while ((nt1=nt_list1[ilist1])!=INDEX_STOP) {
		hcode = base_hcode | (nt1<<4);
		if (refTableTriLoop[hcode]==bonus_value) {
			if (current_rank == final_rank) {
				ss_sample[0][ii+2]=index2char(ii+2,nt1);
				ss_sample[1][ii+2]='.';
				return;
			}
			else {
				current_rank++;
			}
		}
		ilist1++;
	} 
	fprintf(stderr,"%s:line %d: Triloop sampling failed.%s\n",__FILE__,__LINE__,getJobID());
}

void backtrack_tetraloop(int nb_mut_in_unk,int ii,int jj, int lnt, int uu, int vv, int rnt, char **ss_sample, double bonus_value) {
	int ilist1,ilist2,base_hcode,hcode;
	int *nt_list1,nt1,*nt_list2,nt2;
	int ind_table,final_rank,current_rank;
	int tab_mut_remaining[3][2][2] = {{{0,0}},{{0,1},{1,0}},{{1,1}}};
	int kk,maxk;
	
	tetraloop_counter(nb_mut_in_unk,ii,jj,lnt,uu,vv,rnt);
	ind_table=0;
	while (tetraloop_cmpt_table[2*ind_table]!=bonus_value) {
		ind_table++;
		if (ind_table>=nb_value_tetraloop_table) {
			fprintf(stderr,"%s:line %d: Cannot find bonus value %e into tetraloop table.%s\n",__FILE__,__LINE__,bonus_value,getJobID());
			exit(EXIT_FAILURE);
		}
	}
	
	/* rank of sample tetraloop */
	final_rank = (int)(random_va() * tetraloop_cmpt_table[2*ind_table+1]);
	
	
	base_hcode = (lnt<<10) | (uu<<8) | (vv<<2) | rnt;
	current_rank = 0;
	
	if (nb_mut_in_unk==1) { maxk=2; }
	else { maxk=1; }
	
	for (kk=0;kk<maxk;kk++) {
		ilist1=0;
		nt_list1 = uncst_mutation[tab_mut_remaining[nb_mut_in_unk][kk][0]][char2index(ii+2)];
		while ((nt1=nt_list1[ilist1])!=INDEX_STOP) {
			ilist2=0;
			nt_list2 = uncst_mutation[tab_mut_remaining[nb_mut_in_unk][kk][1]][char2index(ii+3)];
			while ((nt2=nt_list2[ilist2])!=INDEX_STOP) {
				hcode = base_hcode | (nt1<<6) | (nt2<<4);
				
				if (refTableTetraLoop[hcode]==bonus_value) {
					if (current_rank == final_rank) {
						ss_sample[0][ii+2]=index2char(ii+2,nt1);
						ss_sample[1][ii+2]='.';
						ss_sample[0][ii+3]=index2char(ii+3,nt2);
						ss_sample[1][ii+3]='.';
						return;
					}
					else {
						current_rank++;
					}
				}
				ilist2++;
			}
			ilist1++;
		} 
	}  
	fprintf(stderr,"%s:line %d: Tetraloop sampling failed.%s\n",__FILE__,__LINE__,getJobID());
}

/***********************************************************************************************/
/* sampling when (i,j) base pair                                                               */
/***********************************************************************************************/

void backtrackHelix(int nb_mut_remaining,int ii,int jj, int lnt, int rnt, char **ss_sample, double mfe) {
	
	int cloop,nb_mut_motif,nb_nt,nb_max_mut,ii_list,ii_list_open,ii_list_close,config,config_close,config_open;
	int nb_max_mut_unk_loop,nb_mut_unk_loop,nb_min_mut_close_bp,nb_max_mut_close_bp,nb_mut_close_bp,nb_max_mut_open_bp,nb_mut_open_bp;
	int xx,aa,cc,uu,vv,dd,bb,yy,newK,nb_free_nt,len,nn,ll,rr,nb_mut_ext_aa,nb_max_mut_in_bulge;
	int size_bulge,nb_mut_in_bulge,nb_mut_inside_boundary,nt2mut,min_rr,max_rr,ind;
	unsigned int *list_config, *list_config_close, *list_config_open;		
	double base_hairpin_energy,bonus,n_total,n_partial,n_remaining;	
	double anchor = 0;
	int newii,newjj,min_length;
	
#ifdef DEBUG_SAMPLING
	double debug=0,pdebug=0;
#endif
	
#ifdef TRACEBACK
	fprintf(stderr,"Helix backtracking: mut=%d, ii=%d, jj=%d, lnt=%d, rnt=%d, mfe=%f\n",nb_mut_remaining,ii,jj,lnt,rnt,mfe);
#endif
	
	/* check */
	
	if ((!Zs[nb_mut_remaining][ii][jj][lnt][rnt].pf)&&(jj-ii>1)) {
		fprintf(stderr,"%s:line %d:helix sampling for [%d,%d] with %d base pair failed.%s\n%e\n",__FILE__,__LINE__,ii,jj,nb_mut_remaining,getJobID(),mfe);
	}
	
	if (((jj<cutpoint)||(ii>=cutpoint))&&(jj-ii-1 < __hairpin_min_length__)) {
		fprintf(stderr,"%s:line %d:helix sampling failed... Hairpin threshold is not respected at position (%d,%d)%s\n",__FILE__,__LINE__,ii,jj,getJobID());
	}
	
	if ((cst_tape[ii])&&(kronecker(ii,lnt))) {
		fprintf(stderr,"%s:line %d:%s: %d: helix sampling failed... index %d cannot mutate%s\n",__FILE__,__LINE__,__FILE__,__LINE__,ii,getJobID());
		return;
	}
	
	if ((cst_tape[jj])&&(kronecker(jj,rnt))) {
		fprintf(stderr,"%s:line %d: helix sampling failed... index %d cannot mutate%s\n",__FILE__,__LINE__,jj,getJobID());
		return;
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
	
	if ((ii<cutpoint)&&(jj>cutpoint)) { /* inter molecular hairpin */
		
		if (is_hairpin(ii,jj)) {	
			
			if (include_intermolecular_interactions) {
				newK = nb_mut_remaining-kronecker(ii,lnt)-kronecker(jj,rnt);
				if (ii==5 && jj==10) printf(">>>>>>>> newK=%d\n",newK);
				if (len>2) {
					nt2mut = len-2 - cst_segment[ii+1][jj-ii-2];
				}
				else {
					nt2mut = 0;
				}
				
				if ((newK>=0)||(newK<=nt2mut)) { /* check configurations */
					
					int nmut_left,nmut_right;
					double mfe_left,mfe_right;
					int rightmostseq1=(int)(cutpoint-0.5);
					int leftmostseq2=(int)(cutpoint+0.5);
					int valid_index_left,valid_index_right,valid_index;
					
					/* single stranded intermolecular hairpin */
					
					n_remaining = genereMutant(nt2mut,newK);	
					anchor =  mfe_hybrid_initialization_energy;
					
					if ((n_remaining)&&(are_equal(mfe,anchor))) {
						
#ifdef TRACEBACK
						fprintf(stderr,"line %d: sample intermolecular hairpin: k=%d, i=%d, j=%d, mfe=%e, anchor=%e%s\n",__LINE__,newK,ii,jj,mfe,anchor,getJobID());
#endif
						if (len>2) {
							mfe_fill_random_mutations(newK,ii+1, jj-1, ss_sample);
						}
						return;
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
									mfe_left=0;
									if (valid_index_left) {
										mfe_left=minimum_double(Zes[nmut_left][ii+1][rightmostseq1][uu][ss].mfe,mfe_left);
										mfe_left=minimum_double(Ze[nmut_left][ii+1][rightmostseq1][uu][ss].mfe,mfe_left);
									}
									for(ww=0;ww<4;ww++) {
										for(vv=0;vv<4;vv++) {
											mfe_right=0;
											if (valid_index_right) {
												mfe_right=minimum_double(Zes[nmut_right][leftmostseq2][jj-1][ww][vv].mfe,mfe_right);
												mfe_right=minimum_double(Ze[nmut_right][leftmostseq2][jj-1][ww][vv].mfe,mfe_right);
											}
											
											anchor=mfe_left+mfe_right+mfe_hybrid_initialization_energy;
											
											if (((mfe_left)||(mfe_right))&&(are_equal(mfe,anchor))) {
												
#ifdef TRACEBACK
												fprintf(stderr,"line %d: sample intermolecular hairpin: k=%d, i=%d, j=%d, u=%d, v=%d, mfe=%e, anchor=%e%s\n",__LINE__,newK,ii,jj,uu,vv,mfe,anchor,getJobID());
#endif
												if (valid_index_left) {
													backtrackExteriorLoop(nmut_left,ii+1,rightmostseq1,uu,ss,ss_sample,mfe_left);
												}
												if (valid_index_right) {
													backtrackExteriorLoop(nmut_right,leftmostseq2,jj-1,ww,vv,ss_sample,mfe_right);
												}
												return;
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
	else { /* single strand hairpin */
		if (is_hairpin(ii,jj)) {	
			for (uu=0;uu<4;uu++) {
				if (!((cst_tape[ii+1])&&(kronecker(ii+1,uu)))) {
					for (vv=0;vv<4;vv++) {				
						if (!((cst_tape[jj-1])&&(kronecker(jj-1,vv)))) {
							
							base_hairpin_energy = EHairpin(lnt,uu,vv,rnt,len-2);
							newK = nb_mut_remaining-kronecker(ii,lnt)-kronecker(ii+1,uu)-kronecker(jj-1,vv)-kronecker(jj,rnt);
							nt2mut = len-4 - cst_segment[ii+2][jj-ii-4];
							
							if ((newK<0)||(newK>nt2mut)) continue; /* check configurations */
							
							n_total = n_partial = 0;
							
							if (len==5) { /*** triloop ***/
								triloop_counter(newK,ii,jj,lnt,uu,vv,rnt);
								for (ind=0;ind<nb_value_triloop_table;ind++) {
									if ((n_partial = triloop_cmpt_table[2*ind+1])) {
										n_total += n_partial;
										bonus = triloop_cmpt_table[2*ind];
										anchor =  base_hairpin_energy + bonus;
										if (are_equal(mfe,anchor)) {
#ifdef TRACEBACK
											fprintf(stderr,"line %d: bactrack triloop: k=%d, i=%d, j=%d, u=%d, v=%d, mfe=%e, anchor=%e%s\n",__LINE__,newK,ii,jj,uu,vv,mfe,anchor,getJobID());
#endif
											ss_sample[0][ii+1]=index2char(ii+1,uu);
											ss_sample[1][ii+1]='.';
											ss_sample[0][jj-1]=index2char(jj-1,vv);
											ss_sample[1][jj-1]='.';
											backtrack_triloop(newK, ii, jj, lnt, uu, vv, rnt, ss_sample, bonus);
											return;
										}
									}
								}
							}
							else if (len==6) { /*** tetraloop ***/
								tetraloop_counter(newK,ii,jj,lnt,uu,vv,rnt);
								for (ind=0;ind<nb_value_tetraloop_table;ind++) {
									if ((n_partial = tetraloop_cmpt_table[2*ind+1])) {
										n_total += n_partial;
										bonus = tetraloop_cmpt_table[2*ind];
										anchor =  base_hairpin_energy + bonus;
										if (are_equal(mfe,anchor)) {
#ifdef TRACEBACK
											fprintf(stderr,"line %d: sample tetraloop: k=%d, i=%d, j=%d, u=%d, v=%d, mfe=%e, anchor=%e%s\n",__LINE__,newK,ii,jj,uu,vv,mfe,anchor,getJobID());
#endif
											ss_sample[0][ii+1]=index2char(ii+1,uu);
											ss_sample[1][ii+1]='.';
											ss_sample[0][jj-1]=index2char(jj-1,vv);
											ss_sample[1][jj-1]='.';
											backtrack_tetraloop(newK, ii, jj, lnt, uu, vv, rnt, ss_sample, bonus);
											return;
										}
									}
								}
							}
							
							/*** special case ***/
#ifdef OLD_GGG_HAIRPIN							
							if((lnt==INDEX_G)&&(uu==INDEX_G)&&(rnt==INDEX_U)) { /*** GGG loop ***/
								double ggg_newK = newK;
								n_partial = 0;
								if (kronecker(ii+2,INDEX_G)) { /*** 3rd nt needs to mutate ***/
									if ((newK)&&(!(cst_tape[ii+2]))) {
										ggg_newK--; /** remove the mutation for GGG loop **/
										n_partial = genereMutant(nt2mut-1,ggg_newK); 
										n_total += n_partial;
										anchor =  base_hairpin_energy + ggg_hairpin_bonus;
									}
								}
								else { /*** already a GGG loop ***/
									n_partial = genereMutant(nt2mut-1,ggg_newK); 
									n_total += n_partial;
									anchor =  base_hairpin_energy + ggg_hairpin_bonus;
								}
								if ((n_partial)&&(are_equal(mfe,anchor))) {
#ifdef TRACEBACK
									fprintf(stderr,"line %d: sample GGG loop: k=%d, i=%d, j=%d, u=%d, v=%d, mfe=%e, anchor=%e%s\n",__LINE__,newK,ii,jj,uu,vv,mfe,anchor,getJobID());
#endif
									if ((len-2<__hairpin_min_length__)||(len-2>__hairpin_max_length__)) {
										fprintf(stderr,"Incoherency: Hairpin does not have the correct length!! (line %d)%s\n",__LINE__,getJobID());
										return;
									}
									ss_sample[0][ii+1]=index2char(ii+1,uu);
									ss_sample[1][ii+1]='.';
									ss_sample[0][ii+2]=index2char(ii+2,INDEX_G);
									ss_sample[1][ii+2]='.';
									ss_sample[0][jj-1]=index2char(jj-1,vv);
									ss_sample[1][jj-1]='.';
									if (len>5) { /*** fill region if the size of the hairpin is more than 3 ***/
										mfe_fill_random_mutations(ggg_newK,ii+3, jj-2, ss_sample);
									}
									return;
								}	  
							}
#endif
							
							/*** generic case ***/
							
							n_remaining = genereMutant(nt2mut,newK) - n_total;	
							anchor =  base_hairpin_energy;
							if ((n_remaining)&&(are_equal(mfe,anchor))) {
								
#ifdef TRACEBACK
								fprintf(stderr,"line %d: sample hairpin: k=%d, i=%d, j=%d, u=%d, v=%d, mfe=%e, anchor=%e%s\n",__LINE__,newK,ii,jj,uu,vv,mfe,anchor,getJobID());
#endif
								if ((len-2<__hairpin_min_length__)||(len-2>__hairpin_max_length__)) {
									fprintf(stderr,"Incoherency: Hairpin does not have the correct length!! (line %d)%s\n",__LINE__,getJobID());
									return;
								}
								ss_sample[0][ii+1]=index2char(ii+1,uu);
								ss_sample[1][ii+1]='.';
								ss_sample[0][jj-1]=index2char(jj-1,vv);
								ss_sample[1][jj-1]='.';
								mfe_fill_random_mutations(newK,ii+2, jj-2, ss_sample);
								return;
							}	  
						}
					}
				}
			}
		}
	}
	
	/*****************************************************************************************************/
	/** Stacking pairs                                                                                  **/
	/*****************************************************************************************************/
	
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
					return;
				}
				anchor = EStack(lnt,uu,vv,rnt) + Zs[newK][ii+1][jj-1][uu][vv].mfe;
#ifdef DEBUG_SAMPLING
				if ((ii==LEFT_CHECK)&&(jj==RIGHT_CHECK)) debug = EStack(lnt,uu,vv,rnt) + Zs[newK][ii+1][jj-1][uu][vv].mfe;
#endif
				if (are_equal(mfe,anchor)) {
#ifdef TRACEBACK
					fprintf(stderr,"line %d: sample stack: k=%d, i=%d, j=%d, u=%d, v=%d, mfe=%e, anchor=%e%s\n",__LINE__,newK,ii+1,jj-1,uu,vv,mfe,anchor,getJobID());
#endif
					backtrackHelix(newK,ii+1,jj-1,uu,vv,ss_sample,Zs[newK][ii+1][jj-1][uu][vv].mfe);
					return;
				}
			}
			ii_list++;
		}
	}
	
#ifdef DEBUG_SAMPLING
	if ((ii==LEFT_CHECK)&&(jj==RIGHT_CHECK)) printf("sample stack: %e\n",debug-pdebug);
	pdebug = debug;
#endif
	
	/*****************************************************************************************************/
	/** bulge AND multi-loop                                                                            **/
	/*****************************************************************************************************/
	
	nb_mut_inside_boundary = nb_mut_remaining - kronecker(ii,lnt) - kronecker(jj,rnt);
	
	if (nb_mut_inside_boundary>=0) {
		
		if ((ii<cutpoint)&&(cutpoint<jj)) {
			min_length = 4;
			nb_free_nt = minimum(__max_size_bulge__,len - 4);
		} else {
			min_length = __hairpin_min_length__ +4;
			nb_free_nt = minimum(__max_size_bulge__,len - __hairpin_min_length__ - 4);
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
							if (((newii<cutpoint)&&(cutpoint<newjj))||((len>=__hairpin_min_length__+4+size_bulge)&&((newjj<cutpoint)||(newii>=cutpoint)))) {  /* check that loop respects cutpoint */
								
								if (!(kronecker(newii,uu) && cst_tape[newii])||
									!(kronecker(newjj,vv) && cst_tape[newjj])) {
									
									nt2mut=size_bulge-cst_segment[ii+1][size_bulge-1];
									nb_max_mut_in_bulge = minimum(nt2mut,nb_mut_inside_boundary);
									
									for (nb_mut_in_bulge=0;nb_mut_in_bulge<=nb_max_mut_in_bulge;nb_mut_in_bulge++) { /* number of mutations in bulge */
										newK = nb_mut_inside_boundary  - nb_mut_in_bulge;
										if (newK>=0) {
											if ((Zs[newK][newii][newjj][uu][vv].pf)&&(genereMutant(nt2mut,nb_mut_in_bulge))) {
												anchor = EBulge(lnt,uu,vv,rnt,size_bulge) + Zs[newK][newii][newjj][uu][vv].mfe;
#ifdef DEBUG_SAMPLING
												if ((ii==LEFT_CHECK)&&(jj==RIGHT_CHECK)) debug = EBulge(lnt,uu,vv,rnt,size_bulge) + Zs[newK][newii][newjj][uu][vv].mfe;
#endif
												if (are_equal(mfe,anchor)) {
#ifdef TRACEBACK
													fprintf(stderr,"line %d: sample bulge: k=%d, i=%d, j=%d, u=%d, v=%d, mfe=%e, anchor=%e%s\n",__LINE__,newK,newii,newjj,uu,vv,mfe,anchor,getJobID());
#endif
													mfe_fill_random_mutations(nb_mut_in_bulge,ii+1,ii+size_bulge,ss_sample);
													backtrackHelix(newK,newii,newjj,uu,vv,ss_sample,Zs[newK][ii+1+size_bulge][jj-1][uu][vv].mfe);
													return;
												}
											}
										}
									}
								}
								
							}
							
							/* bulge on the right */
							
							newii = ii+1;
							newjj = jj-1-size_bulge;
							if (((newii<cutpoint)&&(cutpoint<newjj))||((len>=__hairpin_min_length__+4+size_bulge)&&((newjj<cutpoint)||(newii>=cutpoint)))) {  /* check that loop respects cutpoint */
								
								if (!(kronecker(newii,uu) && cst_tape[newii])||
									!(kronecker(newjj,vv) && cst_tape[newjj])) {
									
									nt2mut=size_bulge-cst_segment[jj-size_bulge][size_bulge-1];
									nb_max_mut_in_bulge = minimum(nt2mut,nb_mut_inside_boundary);
									
									for (nb_mut_in_bulge=0;nb_mut_in_bulge<=nb_max_mut_in_bulge;nb_mut_in_bulge++) { /* number of mutations in bulge */
										newK = nb_mut_inside_boundary  - nb_mut_in_bulge;
										if (newK>=0) {
											if ((Zs[newK][newii][newjj][uu][vv].pf)&&(genereMutant(nt2mut,nb_mut_in_bulge))) {
												anchor = EBulge(lnt,uu,vv,rnt,size_bulge) + Zs[newK][newii][newjj][uu][vv].mfe;
#ifdef DEBUG_SAMPLING
												if ((ii==LEFT_CHECK)&&(jj==RIGHT_CHECK)) debug = EBulge(lnt,uu,vv,rnt,size_bulge) +
													Zs[newK][newii][newjj][uu][vv].mfe;
#endif
												if (are_equal(mfe,anchor)) {
#ifdef TRACEBACK
													fprintf(stderr,"line %d: sample bulge: k=%d, i=%d, j=%d, u=%d, v=%d, mfe=%e, anchor=%e%s\n",__LINE__,newK,newii,newjj,uu,vv,mfe,anchor,getJobID());
#endif
													mfe_fill_random_mutations(nb_mut_in_bulge,jj-size_bulge,jj-1,ss_sample);
													backtrackHelix(newK,ii+1,jj-1-size_bulge,uu,vv,ss_sample,Zs[newK][newii][newjj][uu][vv].mfe);
													return;
												}
											}
										}
									}
								}
							}
						}
					}
					
					/* close multi-loop*/
					
					anchor = EEndMultiLoop() + Zm[nb_mut_inside_boundary][ii+1][jj-1][uu][vv].mfe;
#ifdef DEBUG_SAMPLING
					if ((ii==LEFT_CHECK)&&(jj==RIGHT_CHECK)) debug = EEndMultiLoop() + Zm[nb_mut_inside_boundary][ii+1][jj-1][uu][vv].mfe;
#endif
					if (are_equal(mfe,anchor)) {
#ifdef TRACEBACK
						fprintf(stderr,"line %d: sample multi-loop: k=%d, i=%d, j=%d, u=%d, v=%d, mfe=%e, anchor=%e%s\n",__LINE__,nb_mut_inside_boundary,ii+1,jj-1,uu,vv,mfe,anchor,getJobID());
#endif
						backtrackMultiLoop(nb_mut_inside_boundary,ii+1,jj-1,uu,vv,0,ss_sample,Zm[nb_mut_inside_boundary][ii+1][jj-1][uu][vv].mfe);
						return;
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
									return;
								}
								anchor = EInternal_1x1(lnt,aa,uu,vv,bb,rnt) + Zs[newK][ii+2][jj-2][uu][vv].mfe;
#ifdef DEBUG_SAMPLING
								if ((ii==LEFT_CHECK)&&(jj==RIGHT_CHECK)) debug = EInternal_1x1(lnt,aa,uu,vv,bb,rnt) + Zs[newK][ii+2][jj-2][uu][vv].mfe;
#endif
								if (are_equal(mfe,anchor)) {
									ss_sample[0][ii+1] = index2char(ii+1,aa);
									ss_sample[1][ii+1] = '.';
									ss_sample[0][jj-1] = index2char(jj-1,bb);
									ss_sample[1][jj-1] = '.';
#ifdef TRACEBACK
									fprintf(stderr,"line %d: sample internal-loop: k=%d, i=%d, j=%d, u=%d, v=%d, mfe=%e, anchor=%e%s\n",__LINE__,newK,ii+2,jj-2,uu,vv,mfe,anchor,getJobID());
#endif
									backtrackHelix(newK,ii+2,jj-2,uu,vv,ss_sample,Zs[newK][ii+2][jj-2][uu][vv].mfe);
									return;
								}	      
							}
						}
						break;
					case 1:
						newii=ii+2;
						newjj=jj-3;
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
									return;
								}
								anchor = EInternal_1x2(lnt,aa,uu,vv,cc,bb,rnt) + Zs[newK][ii+2][jj-3][uu][vv].mfe;
#ifdef DEBUG_SAMPLING
								if ((ii==LEFT_CHECK)&&(jj==RIGHT_CHECK)) debug = EInternal_1x2(lnt,aa,uu,vv,cc,bb,rnt) + Zs[newK][ii+2][jj-3][uu][vv].mfe;
#endif
								if (are_equal(mfe,anchor)) {
									ss_sample[0][ii+1] = index2char(ii+1,aa);
									ss_sample[1][ii+1] = '.';
									ss_sample[0][jj-2] = index2char(jj-2,cc);
									ss_sample[1][jj-2] = '.';
									ss_sample[0][jj-1] = index2char(jj-1,bb);
									ss_sample[1][jj-1] = '.';
#ifdef TRACEBACK
									fprintf(stderr,"line %d: sample internal-loop: k=%d, i=%d, j=%d, u=%d, v=%d, mfe=%e, anchor=%e%s\n",__LINE__,newK,ii+2,jj-3,uu,vv,mfe,anchor,getJobID());
#endif
									backtrackHelix(newK,ii+2,jj-3,uu,vv,ss_sample,Zs[newK][ii+2][jj-3][uu][vv].mfe);
									return;
								}	      
							}
						}
						break;
					case 2:
						newii=ii+3;
						newjj=jj-2;
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
									return;
								}
								anchor = EInternal_2x1(lnt,aa,cc,uu,vv,bb,rnt) + Zs[newK][ii+3][jj-2][uu][vv].mfe;
#ifdef DEBUG_SAMPLING
								if ((ii==LEFT_CHECK)&&(jj==RIGHT_CHECK)) debug = EInternal_2x1(lnt,aa,cc,uu,vv,bb,rnt) + Zs[newK][ii+3][jj-2][uu][vv].mfe;
#endif
								if (are_equal(mfe,anchor)) {
									ss_sample[0][ii+1] = index2char(ii+1,aa);
									ss_sample[1][ii+1] = '.';
									ss_sample[0][ii+2] = index2char(ii+2,cc);
									ss_sample[1][ii+2] = '.';
									ss_sample[0][jj-1] = index2char(jj-1,bb);
									ss_sample[1][jj-1] = '.';
#ifdef TRACEBACK
									fprintf(stderr,"line %d: sample internal-loop: k=%d, i=%d, j=%d, u=%d, v=%d, mfe=%e, anchor=%e%s\n",__LINE__,newK,ii+3,jj-2,uu,vv,mfe,anchor,getJobID());
#endif
									backtrackHelix(newK,ii+3,jj-2,uu,vv,ss_sample,Zs[newK][ii+3][jj-2][uu][vv].mfe);
									return;
								}	      
							}
						}
						break;
					case 3:
						newii=ii+3;
						newjj=jj-3;
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
									return;
								}
								anchor = EInternal_2x2(lnt,aa,cc,uu,vv,dd,bb,rnt) + Zs[newK][ii+3][jj-3][uu][vv].mfe;
#ifdef DEBUG_SAMPLING
								if ((ii==LEFT_CHECK)&&(jj==RIGHT_CHECK)) debug = EInternal_2x2(lnt,aa,cc,uu,vv,dd,bb,rnt) + Zs[newK][ii+3][jj-3][uu][vv].mfe;
#endif
								if (are_equal(mfe,anchor)) {
									ss_sample[0][ii+1] = index2char(ii+1,aa);
									ss_sample[1][ii+1] = '.';
									ss_sample[0][ii+2] = index2char(ii+2,cc);
									ss_sample[1][ii+2] = '.';
									ss_sample[0][jj-2] = index2char(jj-2,dd);
									ss_sample[1][jj-2] = '.';
									ss_sample[0][jj-1] = index2char(jj-1,bb);
									ss_sample[1][jj-1] = '.';
#ifdef TRACEBACK
									fprintf(stderr,"line %d: sample internal-loop: k=%d, i=%d, j=%d, u=%d, v=%d, mfe=%e, anchor=%e%s\n",__LINE__,newK,ii+3,jj-3,uu,vv,mfe,anchor,getJobID());
#endif
									backtrackHelix(newK,ii+3,jj-3,uu,vv,ss_sample,Zs[newK][ii+3][jj-3][uu][vv].mfe);
									return;
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
											return;
										}
										
										anchor = EInternal_generic(lnt,aa,cc,uu,vv,dd,bb,rnt,1,nn+2) + Zs[newK][ii+2][jj-3-nn][uu][vv].mfe;
										
#ifdef DEBUG_SAMPLING
										if ((ii==LEFT_CHECK)&&(jj==RIGHT_CHECK)) debug += boltzmannInternal_generic(lnt,aa,cc,uu,vv,dd,bb,rnt,1,nn+2) *
											Zs[newK][ii+2][jj-3-nn][uu][vv].pf * genereMutant(nt2mut,nb_mut_unk_loop);
#endif
										if ((are_equal(mfe,anchor))&&(genereMutant(nt2mut,nb_mut_unk_loop))) {
#ifdef TRACEBACK
											fprintf(stderr,"line %d: sample internal-loop: k=%d, i=%d, j=%d, u=%d, v=%d, mfe=%e, anchor=%e%s\n",__LINE__,newK,ii+2,jj-3-nn,uu,vv,mfe,anchor,getJobID());
#endif
											ss_sample[0][ii+1] = index2char(ii+1,aa);
											ss_sample[1][ii+1] = '.';
											ss_sample[0][jj-2-nn] = index2char(jj-2-nn,dd);
											ss_sample[1][jj-2-nn] = '.';
											ss_sample[0][jj-1] = index2char(jj-1,bb);
											ss_sample[1][jj-1] = '.';
											mfe_fill_random_mutations(nb_mut_unk_loop,jj-1-nn,jj-2,ss_sample);
											backtrackHelix(newK,ii+2,jj-3-nn,uu,vv,ss_sample,Zs[newK][ii+2][jj-3-nn][uu][vv].mfe);
											return;
										}
									}
									ii_list_open++;
								}
							}
						}
					}
					
					/** nx1 **/
					
					newii=ii+3+nn;
					newjj=jj-2;

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
											return;
										}
										anchor = EInternal_generic(lnt,aa,cc,uu,vv,dd,bb,rnt,nn+2,1) + Zs[newK][ii+3+nn][jj-2][uu][vv].mfe;
										
#ifdef DEBUG_SAMPLING
										if ((ii==LEFT_CHECK)&&(jj==RIGHT_CHECK)) debug = EInternal_generic(lnt,aa,cc,uu,vv,dd,bb,rnt,1,nn+2) + Zs[newK][ii+3+nn][jj-2][uu][vv].mfe;
#endif
										if ((are_equal(mfe,anchor))&&(genereMutant(nt2mut,nb_mut_unk_loop))) {
#ifdef TRACEBACK
											fprintf(stderr,"line %d: sample internal-loop: k=%d, i=%d, j=%d, u=%d, v=%d, mfe=%e, anchor=%e%s\n",__LINE__,newK,ii+3+nn,jj-2,uu,vv,mfe,anchor,getJobID());
#endif
											ss_sample[0][ii+1] = index2char(ii+1,aa);
											ss_sample[1][ii+1] = '.';
											ss_sample[0][ii+2+nn] = index2char(ii+2+nn,cc);
											ss_sample[1][ii+2+nn] = '.';
											ss_sample[0][jj-1] = index2char(jj-1,bb);
											ss_sample[1][jj-1] = '.';
											mfe_fill_random_mutations(nb_mut_unk_loop,ii+2,ii+1+nn,ss_sample);
											backtrackHelix(newK,ii+3+nn,jj-2,uu,vv,ss_sample,Zs[newK][ii+3+nn][jj-2][uu][vv].mfe);
											return;
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
	nb_min_mut_close_bp=0;
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
											return;
										}
										
										anchor = EInternal_generic(lnt,aa,cc,uu,vv,dd,bb,rnt,ll+2,rr+2) + Zs[newK][ii+3+ll][jj-3-rr][uu][vv].mfe;
										
#ifdef DEBUG_SAMPLING
										if ((ii==LEFT_CHECK)&&(jj==RIGHT_CHECK)) debug = EInternal_generic(lnt,aa,cc,uu,vv,dd,bb,rnt,ll+2,rr+2) +
											Zs[newK][ii+3+ll][jj-3-rr][uu][vv].mfe;
#endif
#ifdef DEBUG_MFE
										if ((ii==LEFT_CHECK)&&(jj==RIGHT_CHECK)) {
											fprintf(stderr,"internal loop :(%d,%d,%d,%d,%d,%d,%d,%d):[%d,%d]:%e\n",lnt,aa,cc,uu,vv,dd,bb,rnt,ii+3+ll,jj-3-rr,EInternal_generic(lnt,aa,cc,uu,vv,dd,bb,rnt,ll+2,rr+2)+Zs[newK][ii+3+ll][jj-3-rr][uu][vv].mfe);
										}
#endif
										if ((are_equal(mfe,anchor))&&(genereMutant(nt2mut,nb_mut_unk_loop))) {
											int mut_left, mut_right;
#ifdef TRACEBACK
											fprintf(stderr,"line %d: sample internal-loop: k=%d, ll=%d, rr=%d, i=%d, j=%d, u=%d, v=%d, mfe=%e, anchor=%e%s\n",__LINE__,newK,ll,rr,ii+3+ll,jj-3-rr,uu,vv,mfe,anchor,getJobID());
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
											mut_left = partition_mutations_uniform(ii+2, ii+1+ll, jj-1-rr, jj-2, nb_mut_unk_loop);
											mut_right = nb_mut_unk_loop - mut_left;
#ifdef TRACEBACK
											fprintf(stdout ,"mutations in loop: %d, left: %d, right: %d.%s\n",
													nb_mut_unk_loop,mut_left,mut_right,getJobID());
#endif
											mfe_fill_random_mutations(mut_left,ii+2,ii+1+ll,ss_sample);
											mfe_fill_random_mutations(mut_right,jj-1-rr,jj-2,ss_sample);
											backtrackHelix(newK,ii+3+ll,jj-3-rr,uu,vv,ss_sample,Zs[newK][ii+3+ll][jj-3-rr][uu][vv].mfe);
											return;
										}
										ii_list_open++;
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
	
	fprintf(stderr,"helix backtracking failed... Might be due to numerical precision (i=%d, j=%d, nb_mut_remaining=%d, mfe=%e, partition number=%e)%s\n",
			ii,jj,nb_mut_remaining,mfe,anchor,getJobID());
	
}

/***********************************************************************************************/
/* sample multi-loop                                                                           */
/***********************************************************************************************/

void backtrackMultiLoop(int nb_mut_remaining,int ii,int jj, int lnt, int rnt, int last_helix, char **ss_sample, double mfe) {
	
	int uu,vv,rr,ee,newK,err_bound,nb_mut_inside_boundary,nt2mut,min_gap;
	double anchor;
	
#ifdef TRACEBACK
	fprintf(stderr,"MultiLoop backtrack: mut=%d, ii=%d, jj=%d, lnt=%d, rnt=%d, mfe=%f%s\n",nb_mut_remaining,ii,jj,lnt,rnt,mfe,getJobID());
#endif
	
	/* sample a single helix */
	
	anchor = Zs[nb_mut_remaining][ii][jj][lnt][rnt].mfe  + penalty_close_helix(lnt,rnt);
	if (are_equal(mfe,anchor)) {
#ifdef TRACEBACK
		fprintf(stderr,"line %d: sample helix: k=%d, i=%d, j=%d, u=%d, v=%d, mfe=%e up2=%e, %e, %e%s\n",
				__LINE__,nb_mut_remaining,ii,jj,lnt,rnt,mfe,Zs[nb_mut_remaining][ii][jj][lnt][rnt].mfe,
				Zms[nb_mut_remaining][ii][jj][lnt][rnt].mfe,Zm[nb_mut_remaining][ii][jj][lnt][rnt].mfe, getJobID());
#endif
		backtrackHelix(nb_mut_remaining,ii,jj,lnt,rnt,ss_sample,Zs[nb_mut_remaining][ii][jj][lnt][rnt].mfe);
		return;
	}  
	
	/* try intermediate pairings in order to build multi-loop */
	
	nb_mut_inside_boundary = nb_mut_remaining - kronecker(ii,lnt);
	
	if ((ii<cutpoint)&&(cutpoint<jj)) { min_gap = 0; }
	else { min_gap = __hairpin_min_length__; }
	
	for (rr=ii+1;rr<=jj-min_gap-1;rr++) {
		for (uu=0;uu<4;uu++) {
			
			/* initiate array and extend on left-hand side */
			
			if (last_helix) { /* sample the last helix */
				
				if ((!((cst_tape[ii])&&(kronecker(ii,lnt))))&&(validBasePair(uu,rnt))) {
					
					if ((jj-rr>__hairpin_min_length__)||((rr<cutpoint)&&(cutpoint<jj))) {
						if (rr-ii>1) {
							
							nt2mut = rr-ii-1 - cst_segment[ii+1][rr-ii-2];
							err_bound = minimum(nt2mut,nb_mut_inside_boundary);
							
							for (ee=0;ee<=err_bound;ee++) {
								
								newK = nb_mut_remaining - kronecker(ii,lnt) - ee;
								if ((newK<0)||(newK>nb_mut_remaining)) {
									fprintf(stderr,"Incoherency: newK CANNOT be negative!! (line %d)%s\n",__LINE__,getJobID());
									return;
								}
								
								anchor = Zs[newK][rr][jj][uu][rnt].mfe + penalty_close_helix(uu,rnt) + EExtendMultiLoop(rr-ii);
								if ((are_equal(mfe,anchor))&&(genereMutant(nt2mut,ee))) {
									ss_sample[0][ii]=index2char(ii,lnt);
									ss_sample[1][ii]='.';
									mfe_fill_random_mutations(ee,ii+1,rr-1,ss_sample);
#ifdef TRACEBACK
									fprintf(stderr,"line %d: sample helix: k=%d, i=%d, j=%d, u=%d, v=%d%s\n",__LINE__,newK,rr,jj,uu,rnt,getJobID());
#endif
									backtrackHelix(newK,rr,jj,uu,rnt,ss_sample,Zs[newK][rr][jj][uu][rnt].mfe);
									return;
								}
							}
						}
						else { /* special case: only one residue on the left */
							newK = nb_mut_remaining - kronecker(ii,lnt);
							if ((newK<0)||(newK>nb_mut_remaining)) {
								fprintf(stderr,"Incoherency: newK CANNOT be negative!! (line %d)%s\n",__LINE__,getJobID());
								return;
							}
							anchor = Zs[newK][rr][jj][uu][rnt].mfe + penalty_close_helix(uu,rnt) + EExtendMultiLoop(1);
							if (are_equal(mfe,anchor)) {
								ss_sample[0][ii]=index2char(ii,lnt);
								ss_sample[1][ii]='.';
#ifdef TRACEBACK
								fprintf(stderr,"line %d: sample helix: k=%d, i=%d, j=%d, u=%d, v=%d%s\n",__LINE__,newK,rr,jj,uu,rnt,getJobID());
#endif
								backtrackHelix(newK,rr,jj,uu,rnt,ss_sample,Zs[newK][rr][jj][uu][rnt].mfe);
								return;
							}
						}
					}
				}
			}
			else {
				
				if (((rr-ii>__hairpin_min_length__+1)||((ii<cutpoint)&&(cutpoint<rr-1)))
					&&((jj-rr>__hairpin_min_length__)||((rr<cutpoint)&&(cutpoint<jj)))) {
					for (vv=0;vv<4;vv++) {
						if (validBasePair(vv,rnt)) {
							for (ee=0;ee<=nb_mut_remaining;ee++) {
								/* more than 2 helices */
								anchor = Zm[ee][ii][rr-1][lnt][uu].mfe + Zs[nb_mut_remaining-ee][rr][jj][vv][rnt].mfe + penalty_close_helix(vv,rnt) + EAddStemInMultiLoop(1);
								if (are_equal(mfe,anchor)) {
#ifdef TRACEBACK
									fprintf(stderr,"line %d: sample multi-loop: k=%d, i=%d, j=%d, u=%d, v=%d%s\n",__LINE__,ee,ii,rr-1,lnt,uu,getJobID());
									fprintf(stderr,"AND sample helix: k=%d, i=%d, j=%d, u=%d, v=%d%s\n",nb_mut_remaining-ee,rr,jj,vv,rnt,getJobID());
#endif
									backtrackMultiLoop(ee,ii,rr-1,lnt,uu,0,ss_sample,Zm[ee][ii][rr-1][lnt][uu].mfe);
									backtrackHelix(nb_mut_remaining-ee,rr,jj,vv,rnt,ss_sample,Zs[nb_mut_remaining-ee][rr][jj][vv][rnt].mfe);
									return;
								}
								/* exactly 2 helices */
								anchor = Zms[ee][ii][rr-1][lnt][uu].mfe + Zs[nb_mut_remaining-ee][rr][jj][vv][rnt].mfe + penalty_close_helix(vv,rnt) + EAddStemInMultiLoop(2);
								if (are_equal(mfe,anchor)) {
#ifdef TRACEBACK
									fprintf(stderr,"line %d: sample last multi-loop: k=%d, i=%d, j=%d, u=%d, v=%d%s\n",__LINE__,ee,ii,rr-1,lnt,uu,getJobID());
									fprintf(stderr,"AND sample helix: k=%d, i=%d, j=%d, u=%d, v=%d%s\n",nb_mut_remaining-ee,rr,jj,vv,rnt,getJobID());
#endif
									backtrackMultiLoop(ee,ii,rr-1,lnt,uu,1,ss_sample,Zms[ee][ii][rr-1][lnt][uu].mfe);
									backtrackHelix(nb_mut_remaining-ee,rr,jj,vv,rnt,ss_sample,Zs[nb_mut_remaining-ee][rr][jj][vv][rnt].mfe);
									return;
								}
							}
						}
					}
				}
			}
		}
	}
	
	/* extend all tables with an unpaired nucleotide on right side */
	
	if (!((cst_tape[jj])&&(kronecker(jj,rnt)))) {
		for (uu=0;uu<4;uu++) {
			/* right */
			newK = nb_mut_remaining - kronecker(jj,rnt);
			if ((newK<0)||(newK>nb_mut_remaining)) {
				fprintf(stderr,"Incoherency: newK CANNOT be negative!! (line %d)%s\n",__LINE__,getJobID());
				return;
			}
			anchor = EExtendMultiLoop(1) + Zms[newK][ii][jj-1][lnt][uu].mfe;
			if (are_equal(mfe,anchor)) {
				ss_sample[0][jj]=index2char(jj,rnt);
				ss_sample[1][jj]='.';
#ifdef TRACEBACK
				fprintf(stderr,"line %d: sample multi-loop: k=%d, i=%d, j=%d, u=%d, v=%d%s\n",__LINE__,newK,ii,jj-1,lnt,uu,getJobID());
#endif
				backtrackMultiLoop(newK,ii,jj-1,lnt,uu,last_helix,ss_sample,Zms[newK][ii][jj-1][lnt][uu].mfe);
				return;
			}
			anchor = EExtendMultiLoop(1) + Zm[newK][ii][jj-1][lnt][uu].mfe;
			if (are_equal(mfe,anchor)) {
				ss_sample[0][jj]=index2char(jj,rnt);
				ss_sample[1][jj]='.';
#ifdef TRACEBACK
				fprintf(stderr,"line %d: sample multi-loop: k=%d, i=%d, j=%d, u=%d, v=%d%s\n",__LINE__,newK,ii,jj-1,lnt,uu,getJobID());
#endif
				backtrackMultiLoop(newK,ii,jj-1,lnt,uu,last_helix,ss_sample,Zm[newK][ii][jj-1][lnt][uu].mfe);
				return;
			}
		}
	}
	
	/* check. this line should NOT be reachable */
	
	fprintf(stderr,"Multi-loop bactracking failed... Might be due to numerical precision (i=%d, j=%d, nb_mut_remaining=%d, mfe=%f, partition number=%f)%s\n",
			ii,jj,nb_mut_remaining,mfe,anchor,getJobID());
}  

/******************************************************************************************************************/

void backtrackExteriorLoop(int nb_mut_remaining,int ii,int jj, int lnt, int rnt, char **ss_sample, double mfe) {
	
	int uu,vv,rr,ee,newK,err_bound,len=jj-ii+1,nb_mut_inside_boundary,nt2mut,min_gap;
	double anchor;
	
#ifdef TRACEBACK
	fprintf(stderr,"Exterior helix bactracking: mut=%d, ii=%d, jj=%d, lnt=%d, rnt=%d, mfe=%f%s\n",nb_mut_remaining,ii,jj,lnt,rnt,mfe,getJobID());
#endif
	
	/* sample a single helix */
	
	anchor = Zs[nb_mut_remaining][ii][jj][lnt][rnt].mfe + penalty_close_helix(lnt,rnt);
	if (are_equal(mfe,anchor)) {
#ifdef TRACEBACK
		fprintf(stderr,"line %d: sample helix: k=%d, i=%d, j=%d, u=%d, v=%d%s\n",__LINE__,nb_mut_remaining,ii,jj,lnt,rnt,getJobID());
#endif
		backtrackHelix(nb_mut_remaining,ii,jj,lnt,rnt,ss_sample,Zs[nb_mut_remaining][ii][jj][lnt][rnt].mfe);
		return;
	}  
	
	/* try intermediate pairings in order to build multi-loop */
	
	
	nb_mut_inside_boundary = nb_mut_remaining - kronecker(ii,lnt);
	
	if ((ii<cutpoint)&&(cutpoint<jj)) { min_gap = 0; }
	else { min_gap = __hairpin_min_length__; }
	
	for (rr=ii+1;rr<=jj-min_gap-1;rr++) {
		for (uu=0;uu<4;uu++) {
			
			/* initiate array and extend on left-hand side */
			
			if ((!((cst_tape[ii])&&(kronecker(ii,lnt))))&&(validBasePair(uu,rnt))) {
				
				if ((jj-rr>__hairpin_min_length__)||((rr<cutpoint)&&(cutpoint<jj))) {
					if (rr-ii>1) {
						
						nt2mut = rr-ii-1 - cst_segment[ii+1][rr-ii-2];
						err_bound = minimum(nt2mut,nb_mut_inside_boundary);
						
						for (ee=0;ee<=err_bound;ee++) {
							newK = nb_mut_remaining - kronecker(ii,lnt) - ee;
							anchor = Zs[newK][rr][jj][uu][rnt].mfe + penalty_close_helix(uu,rnt);
							if ((are_equal(mfe,anchor))&&(genereMutant(nt2mut,ee))) {
								ss_sample[0][ii]=index2char(ii,lnt);
								ss_sample[1][ii]='.';
								mfe_fill_random_mutations(ee,ii+1, rr, ss_sample);
#ifdef TRACEBACK
								fprintf(stderr,"line %d: sample helix: k=%d, i=%d, j=%d, u=%d, v=%d%s\n",__LINE__,newK,rr,jj,uu,rnt,getJobID());
#endif
								backtrackHelix(newK,rr,jj,uu,rnt,ss_sample,Zs[newK][rr][jj][uu][rnt].mfe);
								return;
							}
						}
					}
					else { /* special case: only one residue on the left */
						newK = nb_mut_remaining - kronecker(ii,lnt);
						if ((newK<0)||(newK>nb_mut_remaining)) {
							fprintf(stderr,"Incoherency: newK (=%d) CANNOT be negative or larger than %d!! (line %d)%s\n",newK,nb_mut_remaining,__LINE__,getJobID());
							return;
						}
						anchor = Zs[newK][rr][jj][uu][rnt].mfe + penalty_close_helix(uu,rnt);
						if (are_equal(mfe,anchor)) {
							ss_sample[0][ii]=index2char(ii,lnt);
							ss_sample[1][ii]='.';
#ifdef TRACEBACK
							fprintf(stderr,"line %d: sample helix: k=%d, i=%d, j=%d, u=%d, v=%d%s\n",__LINE__,newK,rr,jj,uu,rnt,getJobID());
#endif
							backtrackHelix(newK,rr,jj,uu,rnt,ss_sample,Zs[newK][rr][jj][uu][rnt].mfe);
							return;
						}
					}
				}
			}
			
			if (((rr-ii>__hairpin_min_length__+1)||((ii<cutpoint)&&(cutpoint<rr-1)))
				&&((jj-rr>__hairpin_min_length__)||((rr<cutpoint)&&(cutpoint<jj)))) {
				for (vv=0;vv<4;vv++) {
					if (validBasePair(vv,rnt)) {
						for (ee=0;ee<=nb_mut_remaining;ee++) {
							/* more than 2 helices */
							anchor = Ze[ee][ii][rr-1][lnt][uu].mfe + Zs[nb_mut_remaining-ee][rr][jj][vv][rnt].mfe + penalty_close_helix(vv,rnt);
							if (are_equal(mfe,anchor)) {
#ifdef TRACEBACK
								fprintf(stderr,"line %d: sample exterior loop: k=%d, i=%d, j=%d, u=%d, v=%d%s\n",__LINE__,ee,ii,rr-1,lnt,uu,getJobID());
								fprintf(stderr,"AND sample helix: k=%d, i=%d, j=%d, u=%d, v=%d%s\n",nb_mut_remaining-ee,rr,jj,vv,rnt,getJobID());
#endif
								backtrackExteriorLoop(ee,ii,rr-1,lnt,uu,ss_sample,Ze[ee][ii][rr-1][lnt][uu].mfe);
								backtrackHelix(nb_mut_remaining-ee,rr,jj,vv,rnt,ss_sample,Zs[nb_mut_remaining-ee][rr][jj][vv][rnt].mfe);
								return;
							}
							/* exactly 2 helices */
							anchor = Zes[ee][ii][rr-1][lnt][uu].mfe + Zs[nb_mut_remaining-ee][rr][jj][vv][rnt].mfe + penalty_close_helix(vv,rnt);
							if (are_equal(mfe,anchor)) {
#ifdef TRACEBACK
								fprintf(stderr,"line %d: sample exterior loop: k=%d, i=%d, j=%d, u=%d, v=%d, E=%e%s\n",__LINE__,ee,ii,rr-1,lnt,uu,Zes[ee][ii][rr-1][lnt][uu].mfe,getJobID());
								fprintf(stderr,"AND sample helix: k=%d, i=%d, j=%d, u=%d, v=%d, E=%e%s\n",nb_mut_remaining-ee,rr,jj,vv,rnt,Zs[nb_mut_remaining-ee][rr][jj][vv][rnt].mfe,getJobID());
#endif
								backtrackExteriorLoop(ee,ii,rr-1,lnt,uu,ss_sample,Zes[ee][ii][rr-1][lnt][uu].mfe);
								backtrackHelix(nb_mut_remaining-ee,rr,jj,vv,rnt,ss_sample,Zs[nb_mut_remaining-ee][rr][jj][vv][rnt].mfe);
								return;
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
	
	if (len>min_gap+2) { /* minimum length required */
		if (!((cst_tape[jj])&&(kronecker(jj,rnt)))) {
			for (uu=0;uu<4;uu++) {
				/* right */
				newK = nb_mut_remaining - kronecker(jj,rnt);
				if ((newK<0)||(newK>nb_mut_remaining)) {
					fprintf(stderr,"Incoherency: newK CANNOT be negative!! (line %d)%s\n",__LINE__,getJobID());
					return;
				}
				anchor = Zes[newK][ii][jj-1][lnt][uu].mfe;
				if (are_equal(mfe,anchor)) {
					ss_sample[0][jj]=index2char(jj,rnt);
					ss_sample[1][jj]='.';
#ifdef TRACEBACK
					fprintf(stderr,"line %d: sample exterior loop: k=%d, i=%d, j=%d, u=%d, v=%d%s\n",__LINE__,newK,ii,jj-1,lnt,uu,getJobID());
#endif
					backtrackExteriorLoop(newK,ii,jj-1,lnt,uu,ss_sample,Zes[newK][ii][jj-1][lnt][uu].mfe);
					return;
				}
				anchor = Ze[newK][ii][jj-1][lnt][uu].mfe;
				if (are_equal(mfe,anchor)) {
					ss_sample[0][jj]=index2char(jj,rnt);
					ss_sample[1][jj]='.';
#ifdef TRACEBACK
					fprintf(stderr,"line %d: sample exterior loop: k=%d, i=%d, j=%d, u=%d, v=%d%s\n",__LINE__,newK,ii,jj-1,lnt,uu,getJobID());
#endif
					backtrackExteriorLoop(newK,ii,jj-1,lnt,uu,ss_sample,Ze[newK][ii][jj-1][lnt][uu].mfe);
					return;
				}
			}
		}
	}
	
	/* check. this line should NOT be reachable */
	
	fprintf(stderr,"Exterior loop bactracking failed... Might be due to numerical precision (i=%d, j=%d, nb_mut_remaining=%d, mfe=%f, partition number=%f)%s\n",
			ii,jj,nb_mut_remaining,mfe,anchor,getJobID());
}  

/****************************************************************************************************/
/* start sampling a sequence (i,j)                                                                  */
/****************************************************************************************************/

void startBacktrackKSuperOptimal(FILE *mfe_file_flux,int k, double mfe) {
	char **ss_sample;
	int uu,vv,xx=-1,yy=-1,iLast=rna_len-1;
	char *buffer_seq1=NULL,*buffer_ss1=NULL;
	int length_seq1=0;
	
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
	
	for (uu=0;uu<4;uu++) {
		for (vv=0;vv<4;vv++) {
#if 0
			printf("uu=%d, vv=%d, mfe=%f, Zes=%f, Ze=%f\n",uu,vv,mfe,Zes[k][0][iLast][uu][vv].mfe,Ze[k][0][iLast][uu][vv].mfe);
#endif
			if (are_equal(mfe,Zes[k][0][iLast][uu][vv].mfe)) {
				xx = uu;
				yy = vv;
			}
			if (are_equal(mfe,Ze[k][0][iLast][uu][vv].mfe)) {
				xx = uu;
				yy = vv;
			}
		} 
	}
	
	if ((xx<0)||(yy<0)) {
		if (ss_constraint_is_non_empty()) {
			fprintf(stderr,"Cannot initiate m.f.e. bactracking. Check constraints.%s\n",getJobID());
			return;
		}
		else {
			mfe_fill_random_mutations(k, 0, rna_len-1, ss_sample);
			if (cutpoint) {
				strncpy(buffer_seq1,ss_sample[0],(length_seq1)*sizeof(char));
				strncpy(buffer_ss1,ss_sample[1],(length_seq1)*sizeof(char));
				buffer_seq1[length_seq1]='\0';
				buffer_ss1[length_seq1]='\0';
				fprintf(mfe_file_flux,"> %d-superoptimal structure\n",k);
				fprintf(mfe_file_flux,"%s&%s\t(0.0)\n",buffer_seq1,ss_sample[0]+length_seq1);
				fprintf(mfe_file_flux,"%s&%s\n",buffer_ss1,ss_sample[1]+length_seq1);
				return;
			}
			else {
				fprintf(mfe_file_flux,"> %d-superoptimal structure\n",k);
				fprintf(mfe_file_flux,"%s\t(0.0)\n",ss_sample[0]);
				fprintf(mfe_file_flux,"%s\n",ss_sample[1]);
				return;
			}
		}
	}
	
	/* backtrack */
	
	printf("> %d-superoptimal structure\n",k);
	backtrackExteriorLoop(k,0,iLast,xx,yy,ss_sample,mfe);
	if (cutpoint) {
		strncpy(buffer_seq1,ss_sample[0],(length_seq1)*sizeof(char));
		strncpy(buffer_ss1,ss_sample[1],(length_seq1)*sizeof(char));
		buffer_seq1[length_seq1]='\0';
		buffer_ss1[length_seq1]='\0';
		fprintf(mfe_file_flux,"%s&%s\t(%.2f)\n",buffer_seq1,ss_sample[0]+length_seq1,mfe);
		fprintf(mfe_file_flux,"%s&%s\n",buffer_ss1,ss_sample[1]+length_seq1);
	}
	else {
		fprintf(mfe_file_flux,"%s\t(%.2f)\n",ss_sample[0],mfe);
		fprintf(mfe_file_flux,"%s\n",ss_sample[1]);
	}
	
	/* free tables */
	
	free(ss_sample[1]);
	free(ss_sample[0]);
	free(ss_sample);
		
	if (cutpoint) {
		free(buffer_seq1);
		free(buffer_ss1);
	}
	
}
