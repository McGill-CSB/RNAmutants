int *uniform_random_config(int len, int nb_mutations);
void fill_random_mutations(int nb_mutations, int i, int j, char **ss_sample);
void fill_weighted_random_mutations(int nb_mutations, int i, int j, char **ss_sample);
double samplingHelix(int nb_mut_remaining,int ii,int jj, int lnt, int rnt, char **ss_sample);
double samplingMultiLoop(int nb_mut_remaining,int ii,int jj, int lnt, int rnt, int last_helix, char **ss_sample);
double samplingExteriorLoop(int nb_mut_remaining,int ii,int jj, int lnt, int rnt, char **ss_sample);
int startBasicSampling(int i, int j, char **ss_sample, double *Esample);
void startSamplingKmutant(int k, int i, int j, char **ss_sample, double *Esample);
int basicSamplingEngine(FILE *sample_file_flux,int nos, int stat_flag, int warning_flag, int compatible_neighbors);
int sampleFromFile(FILE *sample_file_flux,const char *commandfile, int stat_flag, int warning_flag, int *compatible_neighbors);
