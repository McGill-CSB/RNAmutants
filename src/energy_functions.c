#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "energy_functions.h"
#include "energy_params_tables.h"
#include "util.h"
 
/* For using Nussinov-Jacobson energy model you need to recompile the whole program.        */
/* You need before to declare a preprocessor variable such that:                            */
/*                                                                                          */
/* #define NUSSINOV_JACOBSON                                                                */

//#define NUSSINOV_JACOBSON

/* physical constants */

/* formal temperature used in Boltzmann function */
/* other constants defined in zuker_functions.h  */

double __temperature=310.15; /* 37 Celsius degrees in Kelvin (273.15 + 37 = 310.15) */

double __RT__ = (double)(DEFAULT_TEMPERATURE * GAS_CONSTANT)/1000.0; /* RT in kCal */
                                                                     /* Constant values are declared in zuker_function.h */

/* By default GAIL rule and C loop are not used even if they are declared in the parameter files */
/* If you want to use then you need to declare the macro:                                        */
/*                                                                                               */
/* #define C_LOOP_ALLOWED                                                                        */
/* #define GAIL_RULE_ALLOWED                                                                     */
/*                                                                                               */
/* and then recompile the software. Obviously, these rules will be used if and only they are     */
/* specified in the parameter files!                                                             */


/* PLEASE, for safety do not declare or edit anything after this line */

#define INDEX_A 0
#define INDEX_C 1
#define INDEX_G 2
#define INDEX_U 3
#define INDEX_STOP -1

/* global variables */

/* from RNAmutants.h */
extern char *input_tape;
extern int rna_len;
extern const int __hairpin_min_length__;
extern const int __hairpin_max_length__;
extern int dangle_used;
extern int uncst_mutation[2][4][4];
extern int nb_value_triloop_table;
extern double *triloop_cmpt_table;
extern double *triloop_weight_table;
extern int nb_value_tetraloop_table;
extern double *tetraloop_cmpt_table;
extern double *tetraloop_weight_table;
extern int *cst_tape;

/* routines used in this files */

inline int min3_int(int a, int b, int c) {
  if (c<b) {
    if (c<a) {
    return c;
    }
  }
  else {
    if (b<a) {
      return b;
    }
  }
  return a;
}

inline double min_double(double a, double b) {
  if (b<a) {
    return b;
  }
  return a;
}

/******************************************************************/
/* compute the Boltzmann contribution of non GC closing base pair */
/******************************************************************/

inline double penalty_close_helix(int a, int b) {
  if (((a==INDEX_G)&&(b==INDEX_C))||
      ((a==INDEX_C)&&(b==INDEX_G))) { /* GC base pair */
    return 0.0;
  } 
  else { /* non GC base pair */
    return au_penalty;
  }  
}

inline double boltzmann_penalty_close_helix(int a, int b) {
	if (((a==INDEX_G)&&(b==INDEX_C))||
		((a==INDEX_C)&&(b==INDEX_G))) { /* GC base pair */
		return 1.0;
	} 
	else { /* non GC base pair */
		return exp(- au_penalty / __RT__);
	}  
}


/******************************************************************/
/* compute the hairpin energy                                     */
/******************************************************************/


double boltzmannHairpin(int xx, int uu, int vv, int yy, int len) {
#ifdef NUSSINOV_JACOBSON
  return 1.0;
#else
  double baseEnergy=0, misStacking = 0;

  /* base energy computed from length */

  if (len<=30) {
    baseEnergy = hairpin[len];

    /* dont use bonus for tri, tetra loop or GGG loop yet */

  }
  else {
    baseEnergy = hairpin[30] + param_large_loop * log((double)(len)/30);
  }

  /* misStacking energy */

  if (len>3) {
    misStacking = tstackh[xx][yy][uu][vv];
  }
  else {
    misStacking =  penalty_close_helix(xx,yy);
  }

  /* return value of the Boltzmann contribution */

  return exp(- (baseEnergy + misStacking) / __RT__);

#endif
}

double EHairpin(int xx, int uu, int vv, int yy, int len) {
#ifdef NUSSINOV_JACOBSON
  return 0.0;
#else
  double baseEnergy=0, misStacking = 0;

  /* base energy computed from length */

  if (len<=30) {
    baseEnergy = hairpin[len];

    /* dont use bonus for tri, tetra loop or GGG loop yet */

  }
  else {
    baseEnergy = hairpin[30] + param_large_loop * log((double)(len)/30);
  }

  /* misStacking energy */

  if (len>3) {
    misStacking = tstackh[xx][yy][uu][vv];
  }
  else {
    misStacking =  penalty_close_helix(xx,yy);
  }

  /* return value of the Boltzmann contribution */

  return baseEnergy + misStacking;

#endif
}

/** special function for tri and tetra loop **/

/** return the number if ggg loop **/

double ggg_loop_counter(int nb_mut, int ileft, int iright, int lnt, int int_lnt, int int_rnt, int rnt) {
  double cmpt = 0;
  int unklen=iright-ileft-4; /* length of unknown region of the hairpin -1 for the G at 3rd position */
  int nb_mut_unk;
  return 0.0;
  if (nb_mut<0) {
    fprintf(stderr,"WARNING:%s:line %d: invalid call of GGG loop cmpt.%s\n",__FILE__,__LINE__,getJobID());
    return 0.0; /** cancel this config **/
  }

  if ((lnt!=INDEX_G)||(int_lnt!=INDEX_G)||(rnt!=INDEX_U)) return 1.0;

  /* number of mutations in unknown region */

  nb_mut_unk = nb_mut - kronecker(ileft+2,INDEX_G);

  if (nb_mut_unk>unklen) {
    fprintf(stderr,"WARNING:%s:line %d: invalid call of GGG loop cmpt.%s\n",__FILE__,__LINE__,getJobID());
    return 0.0; /** cancel this config **/
  }

  if (nb_mut_unk<0) return 0.0; /* cannot build a ggg loop */

  if (unklen) {
    cmpt = genereMutant(nb_mut_unk,unklen); /* number of config */
  }
  else {
    cmpt = 1.0;
  }

  return cmpt;
}

/* the mutation constraint does not need to check because this function is called with nb_mut==1 iff cst_tape==0                              */
/* fill a table such that triloop_cmpt_table[2*ii] store the triloop bonus and triloop_cmpt_table[2*ii+1] the number/weight of configurations */
/* WARNING: not optimized */


int triloop_counter(int nb_mut, int ileft, int iright, int lnt, int int_lnt, int int_rnt, int rnt) {
#ifdef NUSSINOV_JACOBSON
	return 1.0;
#else
	int hcode, ilist=0,nt;
	int *list_nt=NULL;
	int ii, cmpt = 0;
	double val;
	
	if ((nb_mut<0)||(nb_mut>1)) {
		fprintf(stderr,"WARNING:%s:line %d: invalid call of triloop bonus.%s\n",__FILE__,__LINE__,getJobID());
		return 0.0; /** cancel this config **/
	}
	
	/* reset counter */
	for (ii=0;ii<nb_value_triloop_table;ii++) {
		triloop_cmpt_table[2*ii+1]=0.0;
		triloop_weight_table[2*ii+1]=0.0;
	}
	
	/* fill tables */
	if ((nb_mut==0)||(!cst_tape[ileft+2])) {
		ilist=0;
		list_nt=uncst_mutation[nb_mut][char2index(ileft+2)];
		while ((nt=list_nt[ilist])!=INDEX_STOP) {
			hcode = (lnt<<8) | (int_lnt<<6) | (nt<<4) | (int_rnt<<2) | rnt;
			if ((val=refTableTriLoop[hcode])) { /* continue if bonus exists */
				for (ii=0;ii<nb_value_triloop_table;ii++) {
					if (triloop_cmpt_table[2*ii]==val) {
						triloop_cmpt_table[2*ii+1]++;
						triloop_weight_table[2*ii+1]+=nucleotide_bias(ileft+2,nt);
						break;
					}
				}
				if (ii==nb_value_triloop_table) {
					fprintf(stderr,"%s:line %d: Corrupted table. Cannot find index for value %e.%s\n",__FILE__,__LINE__,val,getJobID());
				}
				cmpt++;
			}
			ilist++;
		}
	}
	
	return cmpt;
#endif
}

/* Although, mutational constraints must be checked in tetraloop because when the function is called we do not know which position is frozen. */

int tetraloop_counter(int nb_mut, int ileft, int iright, int lnt, int int_lnt, int int_rnt, int rnt) {
#ifdef NUSSINOV_JACOBSON
	return 1.0;
#else
	int hcode, ilist1=0, ilist2=0, nt1, nt2;
	int *list_nt1=NULL, *list_nt2=NULL;
	double val;
	int tab_mut_remaining[3][2][2] = {{{0,0}},{{0,1},{1,0}},{{1,1}}};
	int kk,maxk,ii,cmpt=0;
	
	if ((nb_mut<0)||(nb_mut>2)) {
		fprintf(stderr,"WARNING:%s:line %d: invalid call of tetraloop bonus.%s\n",__FILE__,__LINE__,getJobID());
		return 0.0; /** cancel this config **/
	}
	
	/* reset counter */
	for (ii=0;ii<nb_value_tetraloop_table;ii++) {
		tetraloop_cmpt_table[2*ii+1]=0.0;
		tetraloop_weight_table[2*ii+1]=0.0;
	}
	
	if (nb_mut==1) { maxk=2; } /* used to decide the number of configurations to try */
	else { maxk=1; }
	
	for (kk=0;kk<maxk;kk++) {
		if ((!tab_mut_remaining[nb_mut][kk][0])||(!cst_tape[ileft+2])) { /* continue if no mutation or unconstrained position */
			ilist1=0;
			list_nt1=uncst_mutation[tab_mut_remaining[nb_mut][kk][0]][char2index(ileft+2)];
			while ((nt1=list_nt1[ilist1])!=INDEX_STOP) {
				if ((!tab_mut_remaining[nb_mut][kk][1])||(!cst_tape[ileft+3])) { /* continue if no mutation or unconstrained position */
					ilist2=0;
					list_nt2=uncst_mutation[tab_mut_remaining[nb_mut][kk][1]][char2index(ileft+3)];
					while ((nt2=list_nt2[ilist2])!=INDEX_STOP) {
						hcode = (lnt<<10) | (int_lnt<<8) | (nt1<<6) | (nt2<<4) | (int_rnt<<2) | rnt;
						if ((val=refTableTetraLoop[hcode])) {
							for (ii=0;ii<nb_value_tetraloop_table;ii++) {
								if (tetraloop_cmpt_table[2*ii]==val) {
									tetraloop_cmpt_table[2*ii+1]++;
									tetraloop_weight_table[2*ii+1]+=nucleotide_bias(ileft+2,nt1)*nucleotide_bias(ileft+3,nt2);
									break;
								}
							}
							if (ii==nb_value_tetraloop_table) {
								fprintf(stderr,"%s:line %d: Corrupted table. Cannot find index for value %e.%s\n",__FILE__,__LINE__,val,getJobID());
							}
							cmpt++;
						}
						ilist2++;
					}
				}
				ilist1++;
			}
		}
	}
	
	return cmpt;
#endif
}

/******************************************************************/
/* compute the stacking energy                                    */
/******************************************************************/

double boltzmannStack(int xx, int uu, int vv, int yy) {
#ifdef NUSSINOV_JACOBSON
  return 1.0;
#else
  return exp(- (stack[xx][yy][uu][vv] / __RT__));
#endif
}

double EStack(int xx, int uu, int vv, int yy) {
#ifdef NUSSINOV_JACOBSON
  return 1.0;
#else
  return stack[xx][yy][uu][vv];
#endif
}

/******************************************************************/
/* compute the bulge energy                                       */
/******************************************************************/

double boltzmannBulge(int xx, int uu, int vv, int yy, int len) {
#ifdef NUSSINOV_JACOBSON
  return 1.0;
#else
  double baseEnergy=0, bonus=0;

  /* base energy computed from length */

  if (len<=30) {
    baseEnergy = bulge[len];

    /* bonus */
    
    if (len==1) {
      bonus=stack[xx][yy][uu][vv];

    }
    else { /* non-GC closing base pair penalties must be added if size > 1 */
      bonus = penalty_close_helix(xx,yy) + penalty_close_helix(uu,vv);
    }
  }
  else { /* generic computation of penalty if size greater than 30 */

    baseEnergy = bulge[30] + param_large_loop * log((double)(len)/30);

    /* non-GC closing base pair penalties must be added if size > 1 */
    bonus = penalty_close_helix(xx,yy) + penalty_close_helix(uu,vv);
  }

  return exp(-(baseEnergy + bonus)/__RT__);
#endif
}

double EBulge(int xx, int uu, int vv, int yy, int len) {
#ifdef NUSSINOV_JACOBSON
  return 0.0;
#else
  double baseEnergy=0, bonus=0;

  /* base energy computed from length */

  if (len<=30) {
    baseEnergy = bulge[len];

    /* bonus */
    
    if (len==1) {
      bonus=stack[xx][yy][uu][vv];

    }
    else { /* non-GC closing base pair penalties must be added if size > 1 */
      bonus = penalty_close_helix(xx,yy) + penalty_close_helix(uu,vv);
    }
  }
  else { /* generic computation of penalty if size greater than 30 */

    baseEnergy = bulge[30] + param_large_loop * log((double)(len)/30);

    /* non-GC closing base pair penalties must be added if size > 1 */
    bonus = penalty_close_helix(xx,yy) + penalty_close_helix(uu,vv);
  }

  return baseEnergy + bonus;
#endif
}

/*********************** generic internal loop **********************/

double boltzmannInternal_generic(int xx, int aa, int cc, int uu, int vv, int dd, int bb, int yy, int llen, int rlen) {
#ifdef NUSSINOV_JACOBSON
  return 1.0;
#else

  int len=llen+rlen,lenDiff;
  double baseEnergy,misStacking,asymetry=0.0;

  /* size dependant penalty */
  
  if (len<=30) {
    baseEnergy = internal[len];
  }
  else {
    baseEnergy = internal[30] + param_large_loop * log((double)(len)/30);
  }

  /* misstacking energy */

  /* general case */
  misStacking = tstacki[xx][yy][aa][bb] + tstacki[vv][uu][dd][cc];
  
  lenDiff = abs(llen-rlen);
  if (lenDiff>0) {
    asymetry = min_double(max_correction_asym_internal_loop,
			  lenDiff * array_ninio_internal_loop[min3_int(4,llen,rlen)]);
  }

  return exp(- (baseEnergy + misStacking + asymetry) / __RT__);
  
#endif
}

double EInternal_generic(int xx, int aa, int cc, int uu, int vv, int dd, int bb, int yy, int llen, int rlen) {
#ifdef NUSSINOV_JACOBSON
  return 0.0;
#else

  int len=llen+rlen,lenDiff;
  double baseEnergy,misStacking,asymetry=0.0;

  /* size dependant penalty */
  
  if (len<=30) {
    baseEnergy = internal[len];
  }
  else {
    baseEnergy = internal[30] + param_large_loop * log((double)(len)/30);
  }

  /* misstacking energy */

  /* general case */
  misStacking = tstacki[xx][yy][aa][bb] + tstacki[vv][uu][dd][cc];
  
  lenDiff = abs(llen-rlen);
  if (lenDiff>0) {
    asymetry = min_double(max_correction_asym_internal_loop,
			  lenDiff * array_ninio_internal_loop[min3_int(4,llen,rlen)]);
  }

  return baseEnergy + misStacking + asymetry;
  
#endif
}

/*******************************************************************************************/

double boltzmannInternal_1x1(int xx, int aa, int uu, int vv, int bb, int yy) {
#ifdef NUSSINOV_JACOBSON
  return 1.0;
#else
  return exp(- sint2[xx][yy][uu][vv][aa][bb] / __RT__);
#endif
}

double EInternal_1x1(int xx, int aa, int uu, int vv, int bb, int yy) {
#ifdef NUSSINOV_JACOBSON
  return 0.0;
#else
  return sint2[xx][yy][uu][vv][aa][bb];
#endif
}

/********/

double boltzmannInternal_1x2(int xx, int aa, int uu, int vv, int cc, int bb, int yy) {
#ifdef NUSSINOV_JACOBSON
  return 1.0;
#else
  return exp(- asint1x2[xx][yy][uu][vv][aa][bb][cc] / __RT__);
#endif
}

double EInternal_1x2(int xx, int aa, int uu, int vv, int cc, int bb, int yy) {
#ifdef NUSSINOV_JACOBSON
  return 0.0;
#else
  return asint1x2[xx][yy][uu][vv][aa][bb][cc];
#endif
}

/********/

double boltzmannInternal_2x1(int xx, int aa, int cc, int uu, int vv, int bb, int yy) {
#ifdef NUSSINOV_JACOBSON
  return 1.0;
#else
  return  exp(- asint1x2[vv][uu][yy][xx][bb][cc][aa] / __RT__);          
#endif
}

double EInternal_2x1(int xx, int aa, int cc, int uu, int vv, int bb, int yy) {
#ifdef NUSSINOV_JACOBSON
  return 0.0;
#else
  return asint1x2[vv][uu][yy][xx][bb][cc][aa];          
#endif
}

/********/

double boltzmannInternal_2x2(int xx, int aa, int cc, int uu, int vv, int dd, int bb, int yy) {
#ifdef NUSSINOV_JACOBSON
  return 1.0;
#else
  return exp(- sint4[xx][yy][uu][vv][aa][cc][dd][bb] / __RT__);          
#endif
}

double EInternal_2x2(int xx, int aa, int cc, int uu, int vv, int dd, int bb, int yy) {
#ifdef NUSSINOV_JACOBSON
  return 0.0;
#else
  return sint4[xx][yy][uu][vv][aa][cc][dd][bb];          
#endif
}

/********/

double boltzmannAddStemInMultiLoop(int nb) {
#ifdef NUSSINOV_JACOBSON
  return 1.0;
#else
  static int first_call_boltzmannAddStemInMultiLoop = 1,ii;
  static double addStemTab[20];

  if (first_call_boltzmannAddStemInMultiLoop) {
    first_call_boltzmannAddStemInMultiLoop = 0;

    addStemTab[0]=1.0;
    for (ii=1;ii<10;ii++) {
      addStemTab[ii] = exp(- ((ii) * helix_penalty_multi_loop) / __RT__);
    }
  }

  if (nb < 20) {
    return  addStemTab[nb];
  }
  else {
    return exp(- (nb * helix_penalty_multi_loop) / __RT__);
  }
#endif
}

double EAddStemInMultiLoop(int nb) {
#ifdef NUSSINOV_JACOBSON
  return 0.0;
#else
  static int first_call_boltzmannAddStemInMultiLoop_mfe = 1,ii;
  static double addStemTab_mfe[20];

  if (first_call_boltzmannAddStemInMultiLoop_mfe) {
    first_call_boltzmannAddStemInMultiLoop_mfe = 0;

    addStemTab_mfe[0]=1.0;
    for (ii=1;ii<10;ii++) {
      addStemTab_mfe[ii] = ii * helix_penalty_multi_loop;
    }
  }

  if (nb < 20) {
    return  addStemTab_mfe[nb];
  }
  else {
    return nb * helix_penalty_multi_loop;
  }
#endif
}

/********/

double boltzmannExtendMultiLoop(int len) {
#ifdef NUSSINOV_JACOBSON
  return 1.0;
#else
  static int first_call_boltzmannExtendMultiLoop = 1,ii;
  static double extendTab[20];

  if (first_call_boltzmannExtendMultiLoop) {
    first_call_boltzmannExtendMultiLoop = 0;

    extendTab[0] = 1.0;
    for (ii=1;ii<20;ii++) {
      extendTab[ii] = exp(- ((ii) * free_base_penalty_multi_loop) / __RT__);
    }
  }

  if (len < 20) {
    return  extendTab[len];
  }
  else {
    return exp(- (len * free_base_penalty_multi_loop) / __RT__);
  }
#endif
}

double EExtendMultiLoop(int len) {
#ifdef NUSSINOV_JACOBSON
  return 0.0;
#else
  static int first_call_boltzmannExtendMultiLoop_mfe = 1,ii;
  static double extendTab_mfe[20];

  if (first_call_boltzmannExtendMultiLoop_mfe) {
    first_call_boltzmannExtendMultiLoop_mfe = 0;

    extendTab_mfe[0] = 1.0;
    for (ii=1;ii<20;ii++) {
      extendTab_mfe[ii] = ii * free_base_penalty_multi_loop;
    }
  }

  if (len < 20) {
    return  extendTab_mfe[len];
  }
  else {
    return len * free_base_penalty_multi_loop;
  }
#endif
}

/********/

double boltzmannEndMultiLoop() {
#ifdef NUSSINOV_JACOBSON
  return 1.0;
#else
  static int first_call_boltzmannEndMultiLoop = 1;
  static double endValue;

  if (first_call_boltzmannEndMultiLoop) {
    first_call_boltzmannEndMultiLoop = 0;
    endValue = exp(- (offset_multi_loop + helix_penalty_multi_loop) / __RT__);
  }

  return endValue;
#endif
}

double EEndMultiLoop() {
#ifdef NUSSINOV_JACOBSON
  return 0.0;
#else
  static int first_call_boltzmannEndMultiLoop_mfe = 1;
  static double endValue_mfe;

  if (first_call_boltzmannEndMultiLoop_mfe) {
    first_call_boltzmannEndMultiLoop_mfe = 0;
    endValue_mfe = offset_multi_loop + helix_penalty_multi_loop;
  }

  return endValue_mfe;
#endif
}

