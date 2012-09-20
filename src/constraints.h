/*
 *  constraints.h
 *  
 *
 *  Created by Jerome Waldispuhl on 06/21/2010.
 *  Copyright 2010 McGill. All rights reserved.
 *
 */

/* predicats */

int ss_constraint_is_non_empty();
int is_hairpin(int i, int j);
int is_basepair(int i, int j);
int is_stack(int i, int j);
int is_internal_loop(int i, int m, int n, int j);
int is_triloop(int ii, int jj, char **ss_sample);
int is_tetraloop(int ii, int jj, char **ss_sample);
int is_ggg_hairpin(int ii, int jj, char **ss_sample);
int single_strand_length(int i);

/* fill table functions */

void init_structure_cst_tables(int length_seq);
void clean_structure_cst_tables(int length_seq);
void fill_structure_cst_tables(char *buffer);
int *init_mutation_cst_tables(char *buffer);
int **init_mutation_cst_segment(int *lintab);
int cmpt_open_mutations();
