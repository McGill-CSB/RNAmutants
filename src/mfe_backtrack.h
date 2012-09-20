int *mfe_uniform_random_config(int len, int nb_mutations);
void mfe_fill_random_mutations(int nb_mutations, int i, int j, char **ss_sample);
void backtrackHelix(int nb_mut_remaining,int ii,int jj, int lnt, int rnt, char **ss_sample, double mfe);
void backtrackMultiLoop(int nb_mut_remaining,int ii,int jj, int lnt, int rnt, int last_helix, char **ss_sample, double mfe);
void backtrackExteriorLoop(int nb_mut_remaining,int ii,int jj, int lnt, int rnt, char **ss_sample, double mfe);
void startBacktrackKSuperOptimal(int k, double mfe);
