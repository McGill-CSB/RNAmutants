#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <string.h>
#include <ctype.h>
#include <time.h>
#include "energy_functions.h"
#include "constraints.h"
#include "util.h"

#define INDEX_A 0
#define INDEX_C 1
#define INDEX_G 2
#define INDEX_U 3
#define INDEX_STOP -1


#define NO_DEBUG_MFE
#define NO_DEBUG_SAMPLING
#define K_CHECK 0
#define LEFT_CHECK 16
#define RIGHT_CHECK 20

extern const double __infinity__;

/* sequence constraints (structure constraints are declared in a separate file) */

extern int *cst_tape;
extern int **cst_segment;

/* hairpin parameters */
extern double au_penalty;
extern double ggg_hairpin_bonus;
extern double c_hairpin_slope;
extern double c_hairpin_intersect;
extern double c_hairpin_3;

/** table size values **/

extern int nb_value_triloop_table;
extern double *triloop_cmpt_table;
extern int nb_value_tetraloop_table;
extern double *tetraloop_cmpt_table;
extern double *triloop_weight_table;
extern double *tetraloop_weight_table;

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
extern double __RT__;
extern double cutpoint;
extern int include_intermolecular_interactions;
extern double inter_molecular_init_energy;

/* internal loop lookup table */

extern unsigned int *****loop_tab;
extern unsigned int *****stack_tab;

/*** add-on for main algo ***/

/**** new variables ****/

double mfe_hybrid_initialization_energy = 0;
double pf_hybrid_initialization_energy = 1.0;

/************** update functions ****************/

void update_init_hybrid_param(double value) {
	
	mfe_hybrid_initialization_energy = value;
	pf_hybrid_initialization_energy =	exp(-value/__RT__);
	
}

/*********************************************************************/
/*                                                                   */
/* auxiliar functions                                                */
/*                                                                   */
/*********************************************************************/

/* init function used for Mutant algorithm */

void init_Mutant_arrays_hybrid() { /* data structure can be optimized */
	int nn,ii,jj,kk,ll,local_len;
	
	Zs=(struct Zcell *****)xmalloc(max_mutations*sizeof(struct Zcell ****));
	Zes=(struct Zcell *****)xmalloc(max_mutations*sizeof(struct Zcell *****));
	Zms=(struct Zcell *****)xmalloc(max_mutations*sizeof(struct Zcell ****));
	Zm=(struct Zcell *****)xmalloc(max_mutations*sizeof(struct Zcell ****));
	Ze=(struct Zcell *****)xmalloc(max_mutations*sizeof(struct Zcell ****));
	
	for (nn=0;nn<max_mutations;nn++) {
		
		Zs[nn]=(struct Zcell ****)xmalloc(rna_len*sizeof(struct Zcell ***));
		Zes[nn]=(struct Zcell ****)xmalloc(rna_len*sizeof(struct Zcell ***));
		Zms[nn]=(struct Zcell ****)xmalloc(rna_len*sizeof(struct Zcell ***));
		Zm[nn]=(struct Zcell ****)xmalloc(rna_len*sizeof(struct Zcell ***));
		Ze[nn]=(struct Zcell ****)xmalloc(rna_len*sizeof(struct Zcell ***));
		
		for (ii=0;ii<rna_len;ii++) {
			
			Zs[nn][ii]=(struct Zcell ***)xmalloc(rna_len*sizeof(struct Zcell **));
			Zes[nn][ii]=(struct Zcell ***)xmalloc(rna_len*sizeof(struct Zcell **));
			Zms[nn][ii]=(struct Zcell ***)xmalloc(rna_len*sizeof(struct Zcell **));
			Zm[nn][ii]=(struct Zcell ***)xmalloc(rna_len*sizeof(struct Zcell **));
			Ze[nn][ii]=(struct Zcell ***)xmalloc(rna_len*sizeof(struct Zcell **));
			
			for (jj=0;jj<rna_len;jj++) {
				local_len = jj - ii + 1;
				
				if ((ii>=jj)||
					(((jj<cutpoint)||(ii>=cutpoint))&&(local_len < __hairpin_min_length__))) { /* do not allocate memory for unexpected subsequences */
					Zs[nn][ii][jj]=0;
					Zes[nn][ii][jj]=0;
					Zms[nn][ii][jj]=0;
					Zm[nn][ii][jj]=0;
					Ze[nn][ii][jj]=0;
					continue;
				}
				
				Zs[nn][ii][jj]=(struct Zcell **)xmalloc(4*sizeof(struct Zcell *));
				Zes[nn][ii][jj]=(struct Zcell **)xmalloc(4*sizeof(struct Zcell *));
				Zms[nn][ii][jj]=(struct Zcell **)xmalloc(4*sizeof(struct Zcell *));
				Zm[nn][ii][jj]=(struct Zcell **)xmalloc(4*sizeof(struct Zcell *));
				Ze[nn][ii][jj]=(struct Zcell **)xmalloc(4*sizeof(struct Zcell *));
				
				/* allocate memory for the visiblity fields */
				
				for (kk=0;kk<4;kk++) {
					
					Zs[nn][ii][jj][kk]=(struct Zcell *)xmalloc(4*sizeof(struct Zcell));
					Zes[nn][ii][jj][kk]=(struct Zcell *)xmalloc(4*sizeof(struct Zcell));
					Zms[nn][ii][jj][kk]=(struct Zcell *)xmalloc(4*sizeof(struct Zcell));
					Zm[nn][ii][jj][kk]=(struct Zcell *)xmalloc(4*sizeof(struct Zcell));
					Ze[nn][ii][jj][kk]=(struct Zcell *)xmalloc(4*sizeof(struct Zcell));
					
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
					}
				}
			}
		}
	}
}

/*************** precompute hairpin contribution (can also be integrate in main algo) ********************/

void precomputeHairpin_hybrid(int cutoff, int include_intermolecular_interactions) {
	int kk,ii,jj,xx,yy,uu,vv,ll,maxii,newK,nt2mut,ind;
	double bonus,n_total,n_partial,n_remaining;
	double hairpin_base_energy,mfe_hairpin_base_energy,bound_bias;
	
	update_init_hybrid_param(inter_molecular_init_energy);
	
	for (kk=0;kk<max_mutations;kk++) {
		for (jj=1;jj<rna_len;jj++) {
			maxii=jj;
			for (ii=0;ii<maxii;ii++) {
				if (is_hairpin(ii,jj)) { /* lookup structure constraints */
					
					ll = jj-ii-1;
					
					if ((ii<cutpoint)&&(jj>cutpoint)) {
						int last1,first2,ncst1=0,ncst2=0;
						last1 = (int)(cutpoint);
						first2 = (int)(cutpoint)+1;
						if (ii<last1) {
							ncst1 = cst_segment[ii+1][last1-ii-1];
						}
						if (jj>first2) {
							ncst2 = cst_segment[first2][jj-first2-1];
						}
						nt2mut = ll - ncst1 - ncst2;
					}
					else if ((ll>=__hairpin_min_length__)&&(ll<=__hairpin_max_length__)) {
						nt2mut = ll-2 - cst_segment[ii+2][ll-3];
					}
					else {
						continue;
					}

					for (xx=0;xx<4;xx++) {
						for (yy=0;yy<4;yy++) {
							if (validBasePair(xx,yy)) { /* special case when xx and yy base pair */
								
								/**** enumerate all interior base pair cases ****/
								
								Zs[kk][ii][jj][xx][yy].pf = 0.0;
								Zs[kk][ii][jj][xx][yy].mfe = __infinity__;
								
								if ((ii<cutpoint)&&(jj>cutpoint)) {	/* (ii,jj) is an intermolecular interaction */
									
									if (include_intermolecular_interactions) {
										
										newK=kk-kronecker(ii,xx)-kronecker(jj,yy);
										
										if ((newK>=0)&&(newK<=nt2mut)) {
											/* check cst */
											if (!(((cst_tape[ii])&&(kronecker(ii,xx)))||
												  ((cst_tape[jj])&&(kronecker(jj,yy))))) {

												/*** generic cases (interior bases uu and vv are still considered for further extensions) ***/
												
												/* n_remaining = genereMutant(nt2mut,newK); */
												if (jj-ii>1) {
													n_remaining = sequence_bias(ii+1,jj-1,newK) * nucleotide_bias(ii,xx) * nucleotide_bias(jj,yy);
												}
												else {
													n_remaining = nucleotide_bias(ii,xx) * nucleotide_bias(jj,yy);
												}
												
												Zs[kk][ii][jj][xx][yy].pf += pf_hybrid_initialization_energy * n_remaining;
												Zs[kk][ii][jj][xx][yy].mfe =  minimum_double(Zs[kk][ii][jj][xx][yy].mfe, mfe_hybrid_initialization_energy);
											}	
										}
									}
								}
								else { /* (ii,jj) is an hairpin. length check has already been done */
									
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
																/* n_partial = triloop_cmpt_table[2*ind+1];
																bonus = triloop_cmpt_table[2*ind]; */
																n_partial = triloop_weight_table[2*ind+1] * bound_bias;
																bonus = triloop_weight_table[2*ind];
																n_total += n_partial;
																if (n_partial>0) {
																	Zs[kk][ii][jj][xx][yy].pf +=  hairpin_base_energy * exp(-(bonus/__RT__)) * n_partial;
																	Zs[kk][ii][jj][xx][yy].mfe =  minimum_double(Zs[kk][ii][jj][xx][yy].mfe,mfe_hairpin_base_energy + bonus);
																}																
															}
														} 
														else if (ll==4) { /** tetraloop **/ 
															tetraloop_counter(newK,ii,jj,xx,uu,vv,yy);
															n_total = 0;
															for (ind=0;ind<nb_value_tetraloop_table;ind++) {
																/* n_partial = tetraloop_cmpt_table[2*ind+1];
																bonus = tetraloop_cmpt_table[2*ind]; */
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
															if (kronecker(ii+2,INDEX_G)) { /*** 3rd nt needs to mutate ***/
																if ((newK)&&(!(cst_tape[ii+2]))&&(nt2mut>0)) {
																	/* n_partial = genereMutant(nt2mut-1,newK-1); */
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
																	/* n_partial = genereMutant(nt2mut-1,newK); */
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
														
														/* n_remaining = genereMutant(nt2mut,newK) - n_total; */
														n_remaining = sequence_bias(ii+2,jj-2,newK) * bound_bias - n_total;
														Zs[kk][ii][jj][xx][yy].pf +=  hairpin_base_energy  * n_remaining;
														Zs[kk][ii][jj][xx][yy].mfe =  minimum_double(Zs[kk][ii][jj][xx][yy].mfe,mfe_hairpin_base_energy);

														/* initialize related arrays */
														
														Zms[kk][ii][jj][xx][yy].pf  = Zes[kk][ii][jj][xx][yy].pf  = Zs[kk][ii][jj][xx][yy].pf * boltzmann_penalty_close_helix(xx,yy);
														Zms[kk][ii][jj][xx][yy].mfe = Zes[kk][ii][jj][xx][yy].mfe = Zs[kk][ii][jj][xx][yy].mfe + penalty_close_helix(xx,yy);
														
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
}


/*******************************************************************************/
/*                                                                             */
/* Mutant algorithm for computing of the boltzmann distribution                */
/*                                                                             */
/*******************************************************************************/

void RNAmutants_algorithm_hybrid(int verbose) {
	
	int len, ii, jj, newii, newjj, rr, kk, xx, yy, nn, uu, vv, aa, bb, cc, dd, ee, ll, newK,nb_free_nt,nt2mut;
	unsigned int *list_config, *list_config_close, *list_config_open;
	int cloop,nb_mut_motif,nb_nt,nb_max_mut,ii_list,ii_list_open,ii_list_close,config,config_close,config_open;
	int nb_max_mut_unk_loop,nb_mut_unk_loop,nb_max_mut_close_bp,nb_mut_close_bp,nb_max_mut_open_bp,nb_mut_open_bp;
	int nb_max_mut_in_bulge, nb_mut_in_bulge, size_bulge,nb_mut_inside_boundary,min_rr,max_rr;
	double cmpt_mutants;
	clock_t start=0, current=0;
	double cpu_time_used;
	int leftmostseq2,rightmostseq1;
	int valid_index, min_gap;
	int cut_left = (int)(cutpoint);
	int cut_right = (int)(cutpoint)+1;
	
	rightmostseq1=(int)(cutpoint-0.5);
	leftmostseq2=(int)(cutpoint+0.5);
	
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
			
			if ((ii<cutpoint)&&(jj>cutpoint)) continue; /* only intramolecular interactions */
			
			for (xx=0;xx<4;xx++) {
				
				if (kronecker(ii,xx) && cst_tape[ii]) continue;
				
				for (yy=0;yy<4;yy++) {
					
					if (kronecker(jj,yy) && cst_tape[jj]) continue;
					
					/* right */
					if (!((cst_tape[jj])&&(kronecker(jj,yy)))) {
						for (uu=0;uu<4;uu++) {
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
					
					/* left */
					if (!((cst_tape[ii])&&(kronecker(ii,xx)))) {
						for (uu=0;uu<4;uu++) {										
							newK = kk - kronecker(ii,xx);
							if (newK>=0) {
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
		
		for (len=4;len<=rna_len;len++) { /* the number of nucleic acids in the subsequence */
			
			for (ii=0;ii<=rna_len-len;ii++) { /* ii is the index of the first nucleic acid of the subsequence (indices start at 0) */
				jj=ii+len-1;                    /* jj is the index of the last nucleic acid of the subsequence */
				
				if (((jj<cutpoint)||(ii>=cutpoint))&&(len<__hairpin_min_length__+4)) continue; /* remove short intramolecular interactions */
				
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
				
				
				/* Add potential secondary structure in intermolecular hairpins. This has not been done by precomputation b/c Ze was not available */
				
				
				if (is_basepair(ii,jj)){ /* lookup structure constraints */
					if ((ii<cutpoint)&&(cutpoint<jj)) { /* only intermolecular interaction are concerned */
						
						int valid_index_left=0, valid_index_right=0;
						
						valid_index=0;
						
						if (rightmostseq1-ii>=__hairpin_min_length__+2) {
							valid_index_left=1;
						}
						if (jj-leftmostseq2>=__hairpin_min_length__+2) {
							valid_index_right=1;
						}
						
						for(xx=0;xx<4;xx++) {
							for(yy=0;yy<4;yy++) {
								if (validBasePair(xx,yy)) {
									int	nb_mut_in_internal_region = kk - kronecker(ii,xx) - kronecker(jj,yy);
									int nmut_left,nmut_right;
									double pf_left,pf_right,mfe_right,mfe_left;
									
									
									if (nb_mut_in_internal_region>=0) {
										
										for (nmut_left=0;nmut_left<=nb_mut_in_internal_region;nmut_left++) {
											
											nmut_right=nb_mut_in_internal_region-nmut_left;
											
											pf_left = sequence_bias(ii+1,cut_left,nmut_left)  * nucleotide_bias(ii,xx);
											mfe_left = 0.0;
											pf_right = sequence_bias(cut_right,jj-1,nmut_right)  * nucleotide_bias(jj,yy);
											mfe_right = 0.0;
											
											for(uu=0;uu<4;uu++) {
												for(vv=0;vv<4;vv++) {
													if (valid_index_left) {
														pf_left+=Zes[nmut_left][ii+1][rightmostseq1][uu][vv].pf+Ze[nmut_left][ii+1][rightmostseq1][uu][vv].pf;
														mfe_left=minimum_double(Zes[nmut_left][ii+1][rightmostseq1][uu][vv].mfe,mfe_left);
														mfe_left=minimum_double(Ze[nmut_left][ii+1][rightmostseq1][uu][vv].mfe,mfe_left);
													}
													if (valid_index_right) {
														pf_right+=Zes[nmut_right][leftmostseq2][jj-1][uu][vv].pf+Ze[nmut_right][leftmostseq2][jj-1][uu][vv].pf;
														mfe_right=minimum_double(Zes[nmut_right][leftmostseq2][jj-1][uu][vv].mfe,mfe_right);
														mfe_right=minimum_double(Zes[nmut_right][leftmostseq2][jj-1][uu][vv].mfe,mfe_right);
													}
												}
											}
											
											if ((mfe_left!=0.0)||(mfe_right!=0.0)) { /* cannot efficiently use pf value b/c of bias */
												Zs[kk][ii][jj][xx][yy].pf += pf_left * pf_right * pf_hybrid_initialization_energy * nucleotide_bias(ii,xx) * nucleotide_bias(jj,yy);
												Zs[kk][ii][jj][xx][yy].mfe = minimum_double(mfe_left+mfe_right+mfe_hybrid_initialization_energy,Zs[kk][ii][jj][xx][yy].mfe);
#ifdef DEBUG_MFE
												if ((kk==K_CHECK)&&(ii==LEFT_CHECK)&&(jj==RIGHT_CHECK)) printf("hairpin:(%d,%d)(%d,%d): %e\n",xx,yy,uu,vv,mfe_left+mfe_right+mfe_hybrid_initialization_energy);
#endif
											}
										}
									}
								}
							}
						}
					}
				}
				
				/* lookup stacks */
				
				nb_max_mut=minimum(4-cst_tape[ii]-cst_tape[ii+1]-cst_tape[jj-1]-cst_tape[jj],kk);
				
				if (is_stack(ii,jj)){ /* lookup structure constraints */
					
					valid_index=0;
					
					if ((jj<cutpoint)||(ii>=cutpoint)) {
						if (len>=__hairpin_min_length__+4) {
							valid_index=1;
						}
					}
					else if ((ii<rightmostseq1)&&(leftmostseq2<jj)) {
						valid_index=1;
					}
					
					if (valid_index){ /* just to avoid useless checks */
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
#ifdef DEBUG_MFE
								if ((kk==K_CHECK)&&(ii==LEFT_CHECK)&&(jj==RIGHT_CHECK)) printf("stack:(%d,%d)(%d,%d): %e\n",xx,yy,uu,vv,EStack(xx,uu,vv,yy) + Zs[newK][ii+1][jj-1][uu][vv].mfe);
#endif
								ii_list++;
							}
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
									newii=ii+2;
									newjj=jj-2;
									if (!is_internal_loop(ii,newii,newjj,jj)) break; /* lookup structure constraints */
									if (((newii<cutpoint)&&(newjj>cutpoint))||((len>=__hairpin_min_length__+6)&&((newjj<cutpoint)||(newii>=cutpoint)))) { /* constrain loop to respect cutpoint */
										xx=config>>10;
										aa=(config>>8)&3;
										uu=(config>>6)&3;
										vv=(config>>4)&3;
										bb=(config>>2)&3;
										yy=config&3;
										newK = kk - nb_mut_motif + kronecker(ii+2,uu) + kronecker(jj-2,vv);
										Zs[kk][ii][jj][xx][yy].pf += boltzmannInternal_1x1(xx,aa,uu,vv,bb,yy) * Zs[newK][ii+2][jj-2][uu][vv].pf * nucleotide_bias(ii,xx) * nucleotide_bias(ii+1,aa) * nucleotide_bias(jj-1,bb) * nucleotide_bias(jj,yy);
										Zs[kk][ii][jj][xx][yy].mfe = minimum_double(Zs[kk][ii][jj][xx][yy].mfe,EInternal_1x1(xx,aa,uu,vv,bb,yy) + Zs[newK][ii+2][jj-2][uu][vv].mfe);
 
#ifdef DEBUG_SAMPLING
										if ((ii==LEFT_CHECK)&&(jj==RIGHT_CHECK)) debug += boltzmannInternal_1x1(xx,aa,uu,vv,bb,yy) * Zs[newK][ii+2][jj-2][uu][vv].pf;
#endif
#ifdef DEBUG_MFE
										if ((kk==K_CHECK)&&(ii==LEFT_CHECK)&&(jj==RIGHT_CHECK)) printf("1x1:(%d,%d)(%d,%d): %e\n",xx,yy,uu,vv,EInternal_1x1(xx,aa,uu,vv,bb,yy) + Zs[newK][ii+2][jj-2][uu][vv].mfe);
#endif
									}
									break;
								case 1:
									newii=ii+2;
									newjj=jj-3;
									if (!is_internal_loop(ii,newii,newjj,jj)) break; /* lookup structure constraints */
									if (((newii<cutpoint)&&(newjj>cutpoint))||((len>=__hairpin_min_length__+7)&&((newjj<cutpoint)||(newii>=cutpoint)))) { /* constrain loop to respect cutpoint */
										xx=config>>12;
										aa=(config>>10)&3;
										uu=(config>>8)&3;
										vv=(config>>6)&3;
										cc=(config>>4)&3;
										bb=(config>>2)&3;
										yy=config&3;
										newK = kk - nb_mut_motif + kronecker(ii+2,uu) + kronecker(jj-3,vv);
										Zs[kk][ii][jj][xx][yy].pf += boltzmannInternal_1x2(xx,aa,uu,vv,cc,bb,yy) * Zs[newK][ii+2][jj-3][uu][vv].pf * nucleotide_bias(ii,xx) * nucleotide_bias(ii+1,aa) * nucleotide_bias(jj-2,cc) * nucleotide_bias(jj-1,bb) * nucleotide_bias(jj,yy);
										Zs[kk][ii][jj][xx][yy].mfe = minimum_double(Zs[kk][ii][jj][xx][yy].mfe,EInternal_1x2(xx,aa,uu,vv,cc,bb,yy) + Zs[newK][ii+2][jj-3][uu][vv].mfe);
#ifdef DEBUG_SAMPLING
										if ((ii==LEFT_CHECK)&&(jj==RIGHT_CHECK)) debug += boltzmannInternal_1x2(xx,aa,uu,vv,cc,bb,yy) * Zs[newK][ii+2][jj-3][uu][vv].pf;
#endif
#ifdef DEBUG_MFE
										if ((kk==K_CHECK)&&(ii==LEFT_CHECK)&&(jj==RIGHT_CHECK)) printf("1x2:(%d,%d)(%d,%d): %e\n",xx,yy,uu,vv,EInternal_1x2(xx,aa,uu,vv,cc,bb,yy) + Zs[newK][ii+2][jj-3][uu][vv].mfe);
#endif
									}
									break;
								case 2:
									newii=ii+3;
									newjj=jj-2;
									if (!is_internal_loop(ii,newii,newjj,jj)) break; /* lookup structure constraints */
									if (((newii<cutpoint)&&(newjj>cutpoint))||((len>=__hairpin_min_length__+7)&&((newjj<cutpoint)||(newii>=cutpoint)))) { /* constrain loop to respect cutpoint */
										xx=config>>12;
										aa=(config>>10)&3;
										cc=(config>>8)&3;
										uu=(config>>6)&3;
										vv=(config>>4)&3;
										bb=(config>>2)&3;
										yy=config&3;
										newK = kk - nb_mut_motif + kronecker(ii+3,uu) + kronecker(jj-2,vv);
										Zs[kk][ii][jj][xx][yy].pf += boltzmannInternal_2x1(xx,aa,cc,uu,vv,bb,yy) * Zs[newK][ii+3][jj-2][uu][vv].pf * nucleotide_bias(ii,xx) * nucleotide_bias(ii+1,aa) * nucleotide_bias(ii+2,cc) * nucleotide_bias(jj-1,bb) * nucleotide_bias(jj,yy);
										Zs[kk][ii][jj][xx][yy].mfe = minimum_double(Zs[kk][ii][jj][xx][yy].mfe,EInternal_2x1(xx,aa,cc,uu,vv,bb,yy) + Zs[newK][ii+3][jj-2][uu][vv].mfe);
#ifdef DEBUG_SAMPLING
										if ((ii==LEFT_CHECK)&&(jj==RIGHT_CHECK)) debug += boltzmannInternal_2x1(xx,aa,cc,uu,vv,bb,yy) * Zs[newK][ii+3][jj-2][uu][vv].pf;
#endif
#ifdef DEBUG_MFE
										if ((kk==K_CHECK)&&(ii==LEFT_CHECK)&&(jj==RIGHT_CHECK)) printf("2x1:(%d,%d)(%d,%d): %e\n",xx,yy,uu,vv,EInternal_2x1(xx,aa,cc,uu,vv,bb,yy) + Zs[newK][ii+3][jj-2][uu][vv].mfe);
#endif
									}
									break;
								case 3:
									newii=ii+3;
									newjj=jj-3;
									if (!is_internal_loop(ii,newii,newjj,jj)) break; /* lookup structure constraints */
									if (((newii<cutpoint)&&(newjj>cutpoint))||((len>=__hairpin_min_length__+8)&&((newjj<cutpoint)||(newii>=cutpoint)))) { /* constrain loop to respect cutpoint */
										xx=config>>14;
										aa=(config>>12)&3;
										cc=(config>>10)&3;
										uu=(config>>8)&3;
										vv=(config>>6)&3;
										dd=(config>>4)&3;
										bb=(config>>2)&3;
										yy=config&3;
										newK = kk - nb_mut_motif + kronecker(ii+3,uu) + kronecker(jj-3,vv);
										Zs[kk][ii][jj][xx][yy].pf += boltzmannInternal_2x2(xx,aa,cc,uu,vv,dd,bb,yy) * Zs[newK][ii+3][jj-3][uu][vv].pf * nucleotide_bias(ii,xx) * nucleotide_bias(ii+1,aa) * nucleotide_bias(ii+2,cc) * nucleotide_bias(jj-2,dd) * nucleotide_bias(jj-1,bb) * nucleotide_bias(jj,yy);
										Zs[kk][ii][jj][xx][yy].mfe = minimum_double(Zs[kk][ii][jj][xx][yy].mfe,EInternal_2x2(xx,aa,cc,uu,vv,dd,bb,yy) + Zs[newK][ii+3][jj-3][uu][vv].mfe);
#ifdef DEBUG_SAMPLING
										if ((ii==LEFT_CHECK)&&(jj==RIGHT_CHECK)) debug += boltzmannInternal_2x2(xx,aa,cc,uu,vv,dd,bb,yy) * Zs[newK][ii+3][jj-3][uu][vv].pf;
#endif
#ifdef DEBUG_MFE
										if ((kk==K_CHECK)&&(ii==LEFT_CHECK)&&(jj==RIGHT_CHECK)) printf("2x2:(%d,%d)(%d,%d): %e\n",xx,yy,uu,vv,EInternal_2x2(xx,aa,cc,uu,vv,dd,bb,yy) + Zs[newK][ii+3][jj-3][uu][vv].mfe);
#endif
									}
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
				
				/* partial check further requirement will be conducted later */
				valid_index=1;
				if ((ii<cutpoint)&&(cutpoint<jj)) {
					if ((cutpoint<ii+2)||(jj-2<cutpoint)) {
						valid_index=0;
					}
				}
				else if (jj-ii-5<__hairpin_min_length__) {
					valid_index=0;
				}
				
				if (valid_index) {
					
					if ((ii<cutpoint)&&(cutpoint<jj)) {
						nb_free_nt = minimum(__max_size_bulge__,len - 7);
					} else {
						nb_free_nt = minimum(__max_size_bulge__,len - __hairpin_min_length__ - 7);
					}
					
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
								
								newii=ii+2;
								newjj=jj-3-nn;
								
								if (is_internal_loop(ii,newii,newjj,jj)) { /* lookup structure constraints */
									
									if (((newii<cutpoint)&&(newjj>cutpoint))||((len>=__hairpin_min_length__+7+nn)&&((newjj<cutpoint)||(newii>=cutpoint)))) {
										
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
													
#if 0
													if (kk==2) {
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
																Zs[kk][ii][jj][xx][yy].pf += boltzmannInternal_generic(xx,aa,cc,uu,vv,dd,bb,yy,1,nn+2) * Zs[newK][ii+2][jj-3-nn][uu][vv].pf
																* cmpt_mutants * nucleotide_bias(ii,xx) * nucleotide_bias(ii+1,aa) * nucleotide_bias(jj-2-nn,dd) * nucleotide_bias(jj-1,bb)
																* nucleotide_bias(jj,yy);
																Zs[kk][ii][jj][xx][yy].mfe = minimum_double(Zs[kk][ii][jj][xx][yy].mfe,EInternal_generic(xx,aa,cc,uu,vv,dd,bb,yy,1,nn+2) +
																											Zs[newK][ii+2][jj-3-nn][uu][vv].mfe);
#ifdef DEBUG_SAMPLING
																if ((ii==LEFT_CHECK)&&(jj==RIGHT_CHECK)) debug += boltzmannInternal_generic(xx,aa,cc,uu,vv,dd,bb,yy,1,nn+2) *
																	Zs[newK][ii+2][jj-3-nn][uu][vv].pf * cmpt_mutants;
#endif
#ifdef DEBUG_MFE
																if ((kk==K_CHECK)&&(ii==LEFT_CHECK)&&(jj==RIGHT_CHECK)) printf("1xn:(%d,%d)(%d,%d): %e\n",xx,yy,uu,vv,EInternal_generic(xx,aa,cc,uu,vv,dd,bb,yy,1,nn+2) + Zs[newK][ii+2][jj-3-nn][uu][vv].mfe);
#endif
															}
														}
													}
													ii_list_open++;
												}
											}
										}
									}
								}
								
								
								/** nx1 **/
								
								newii=ii+3+nn;
								newjj=jj-2;
								
								if (is_internal_loop(ii,newii,newjj,jj)) { /* lookup structure constraints */
									
									if (((newii<cutpoint)&&(newjj>cutpoint))||((len>=__hairpin_min_length__+7+nn)&&((newjj<cutpoint)||(newii>=cutpoint)))) {
										
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
																Zs[kk][ii][jj][xx][yy].mfe = minimum_double(Zs[kk][ii][jj][xx][yy].mfe,EInternal_generic(xx,aa,cc,uu,vv,dd,bb,yy,nn+2,1) +
																											Zs[newK][ii+3+nn][jj-2][uu][vv].mfe);
#ifdef DEBUG_SAMPLING
																if ((ii==LEFT_CHECK)&&(jj==RIGHT_CHECK)) debug += boltzmannInternal_generic(xx,aa,cc,uu,vv,dd,bb,yy,nn+2,1) *
																	Zs[newK][ii+3+nn][jj-2][uu][vv].pf * cmpt_mutants;
#endif
#ifdef DEBUG_MFE
																if ((kk==K_CHECK)&&(ii==LEFT_CHECK)&&(jj==RIGHT_CHECK)) printf("nx1:(%d,%d)(%d,%d): %e\n",xx,yy,uu,vv,EInternal_generic(xx,aa,cc,uu,vv,dd,bb,yy,nn+2,1) + Zs[newK][ii+3+nn][jj-2][uu][vv].mfe);
#endif
															}
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
				}
				
#ifdef DEBUG_SAMPLING
				if ((ii==LEFT_CHECK)&&(jj==RIGHT_CHECK)) printf("1xn & nx1: %e\n",debug-pdebug);
				pdebug = debug;
#endif
				
				/** mxn internal loops (m+n>4) **/
				
				valid_index=1;
				
				if ((ii<cutpoint)&&(jj>cutpoint)) {
					nb_free_nt = minimum(2*__max_size_bulge__,len - 8);
					if ((ii+3>cutpoint)||(jj-3<cutpoint)) {
						valid_index=0;
					}
				}
				else {
					nb_free_nt = minimum(2*__max_size_bulge__,len - __hairpin_min_length__ - 8);
					if (jj-ii-5<__hairpin_min_length__) {
						valid_index=0;
					}
				}
				
				
				if (valid_index) {
					
					
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
									
									newii=ii+3+ll;
									newjj=jj-3-rr;
									
									if (!is_internal_loop(ii,newii,newjj,jj)) continue; /* lookup structure constraints */
									
									if (((newii<cutpoint)&&(newjj>cutpoint))||((len>=__hairpin_min_length__+8+nn)&&((newjj<cutpoint)||(newii>=cutpoint)))) { /* check that loop respects cutpoint */
										
										
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
															Zs[kk][ii][jj][xx][yy].mfe = minimum_double(Zs[kk][ii][jj][xx][yy].mfe,EInternal_generic(xx,aa,cc,uu,vv,dd,bb,yy,ll+2,rr+2) +
																										Zs[newK][ii+3+ll][jj-3-rr][uu][vv].mfe);
#ifdef DEBUG_SAMPLING
															if ((ii==LEFT_CHECK)&&(jj==RIGHT_CHECK)) debug += boltzmannInternal_generic(xx,aa,cc,uu,vv,dd,bb,yy,ll+2,rr+2) *
																Zs[newK][ii+3+ll][jj-3-rr][uu][vv].pf * cmpt_mutants;
#endif
#ifdef DEBUG_MFE
															if ((kk==K_CHECK)&&(ii==LEFT_CHECK)&&(jj==RIGHT_CHECK)) printf("mxn:(%d,%d)(%d,%d):[%d,%d]: %e\n",xx,yy,uu,vv,ii+3+ll,jj-3-rr,EInternal_generic(xx,aa,cc,uu,vv,dd,bb,yy,ll+2,rr+2) + Zs[newK][ii+3+ll][jj-3-rr][uu][vv].mfe);
#endif
														}
													}
													ii_list_open++;
												}
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
						
						if (validBasePair(xx,yy)) { /* special case when xx and yy base pair */
							
							nb_mut_inside_boundary = kk - kronecker(ii,xx) - kronecker(jj,yy);
							
							if (nb_mut_inside_boundary>=0) {
								
								if ((len > __hairpin_min_length__ +4)||((len > 4)&&(ii<cutpoint)&&(jj>cutpoint))) { /* min length required for a stack (complete check is completed later) */
									
									for (uu=0;uu<4;uu++) {
										
										for (vv=0;vv<4;vv++) {
											
											if (validBasePair(uu,vv)) {
												
												if ((ii<cutpoint)&&(cutpoint<jj)) {
													nb_free_nt = minimum(__max_size_bulge__,len - 4);
												} else {
													nb_free_nt = minimum(__max_size_bulge__,len - __hairpin_min_length__ - 4);
												}

												/* bulge */
												
												for (size_bulge=1;size_bulge<=nb_free_nt;size_bulge++) { /* size of the bulge */
													
													/* we need to distinguish the bulge because the asymetry of constraints */
													
													/* bulge on the left */
													
													newii = ii+1+size_bulge;
													newjj = jj-1;
													
													if (((newii<cutpoint)&&(cutpoint<newjj))||((len>=__hairpin_min_length__+4+size_bulge)&&((newjj<cutpoint)||(newii>=cutpoint)))) {  /* check that loop respects cutpoint */
														
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
																				Zs[kk][ii][jj][xx][yy].mfe = minimum_double(Zs[kk][ii][jj][xx][yy].mfe,EBulge(xx,uu,vv,yy,size_bulge) + Zs[newK][newii][newjj][uu][vv].mfe);
#ifdef DEBUG_SAMPLING
																				if ((ii==LEFT_CHECK)&&(jj==RIGHT_CHECK)) debug += boltzmannBulge(xx,uu,vv,yy,size_bulge) *
																					Zs[newK][newii][newjj][uu][vv].pf * cmpt_mutants;
#endif
#ifdef DEBUG_MFE
																				if ((kk==K_CHECK)&&(ii==LEFT_CHECK)&&(jj==RIGHT_CHECK)) printf("bulgeL:(%d,%d)(%d,%d): %e\n",xx,yy,uu,vv,EBulge(xx,uu,vv,yy,size_bulge) + Zs[newK][newii][newjj][uu][vv].mfe);
#endif
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

													if (((newii<cutpoint)&&(cutpoint<newjj))||((len>=__hairpin_min_length__+4+size_bulge)&&((newjj<cutpoint)||(newii>=cutpoint)))) {  /* check that loop respects cutpoint */
														
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
																				Zs[kk][ii][jj][xx][yy].mfe = minimum_double(Zs[kk][ii][jj][xx][yy].mfe, EBulge(xx,uu,vv,yy,size_bulge) + Zs[newK][newii][newjj][uu][vv].mfe);
#ifdef DEBUG_SAMPLING
																				if ((ii==LEFT_CHECK)&&(jj==RIGHT_CHECK)) debug += boltzmannBulge(xx,uu,vv,yy,size_bulge) *
																					Zs[newK][newii][newjj][uu][vv].pf * cmpt_mutants;
#endif
#ifdef DEBUG_MFE
																				if ((kk==K_CHECK)&&(ii==LEFT_CHECK)&&(jj==RIGHT_CHECK)) printf("bulgeR:(%d,%d)(%d,%d): %e\n",xx,yy,uu,vv,EBulge(xx,uu,vv,yy,size_bulge) + Zs[newK][newii][newjj][uu][vv].mfe);
#endif
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
											
											if ((ii<jj)&&(is_basepair(ii,jj))) { /* lookup structure constraints */
												
												Zs[kk][ii][jj][xx][yy].pf += boltzmannEndMultiLoop() * Zm[nb_mut_inside_boundary][ii+1][jj-1][uu][vv].pf * nucleotide_bias(ii,xx) * nucleotide_bias(jj,yy);
												Zs[kk][ii][jj][xx][yy].mfe = minimum_double(Zs[kk][ii][jj][xx][yy].mfe,EEndMultiLoop() + Zm[nb_mut_inside_boundary][ii+1][jj-1][uu][vv].mfe);
												
#if 0
												if (Zm[nb_mut_inside_boundary][ii+1][jj-1][uu][vv].pf)
													printf("%d,%d: %e (%e + %e)\n",ii,jj,Zs[kk][ii][jj][xx][yy].mfe,boltzmannEndMultiLoop(),Zm[nb_mut_inside_boundary][ii+1][jj-1][uu][vv].mfe);
#endif
#ifdef DEBUG_SAMPLING
												if ((ii==LEFT_CHECK)&&(jj==RIGHT_CHECK)) debug += boltzmannEndMultiLoop() * Zm[nb_mut_inside_boundary][ii+1][jj-1][uu][vv].pf;
#endif
#ifdef DEBUG_MFE
												if ((kk==K_CHECK)&&(ii==LEFT_CHECK)&&(jj==RIGHT_CHECK)) printf("multi:(%d,%d)(%d,%d): %e\n",xx,yy,uu,vv,EEndMultiLoop() + Zm[nb_mut_inside_boundary][ii+1][jj-1][uu][vv].mfe);
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
						
						if ((ii<cutpoint)&&(cutpoint<jj)) { min_gap = 0; }
						else { min_gap = __hairpin_min_length__; }
						
						for (rr=ii+1;rr<=jj-min_gap-1;rr++) {
							
							for (uu=0;uu<4;uu++) {
								
								/* initiate array and extend on left-hand side */
								
								if (single_strand_length(rr-1)>=(rr-ii-1)) { /* lookup structure constraints */
									
									if ((!((cst_tape[ii])&&(kronecker(ii,xx))))&&(validBasePair(uu,yy))) {
										
										if ((jj-rr>__hairpin_min_length__)||((rr<cutpoint)&&(cutpoint<jj))) {
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
								}
								
								if (((rr-ii>__hairpin_min_length__+1)||((ii<cutpoint)&&(cutpoint<rr-1)))
									&&((jj-rr>__hairpin_min_length__)||((rr<cutpoint)&&(cutpoint<jj)))) {
									for (vv=0;vv<4;vv++) {
										if (validBasePair(vv,yy)) {
											for (ee=0;ee<=kk;ee++) {
												//if ((ii==1)&&(jj==12)) printf("rr=%d, le=%d, re=%d, %e - %e \n",rr,ee,kk-ee,Zms[ee][ii][rr-1][xx][uu].pf,Zs[kk-ee][rr][jj][vv][yy].pf);	
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
							
							if ((ii<cutpoint)&&(cutpoint<jj)) { min_gap = 0; }
							else { min_gap = __hairpin_min_length__; }

							if (len>min_gap+2) { /* minimum length required */
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
	
#ifdef INFO
	printf("\n");
#endif
	
}


