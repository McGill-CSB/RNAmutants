/* generic memory allocation*/

void *xmalloc(int size);

/* mathematical functions */

inline int minimum(int a, int b);
inline double minimum_double(double a, double b);
inline int maximum(int a, int b);
double random_va();

/* RNA functions */

char *emptyRNAss();
void cleanRNAss(char *RNAss);
int canBasePair(int ind1, int ind2);
int validBasePair(int ind1, int ind2);
int convert2index(char carac);
inline int char2index(int pos);
char index2char(int pos, int index);
char writeNt(int pos, int lcase);
char writeNt_with_bias(int pos, int lcase);

/* check that double are equal under precision range */

int checkDouble(double a, double b);

/* math functions */

int kronecker(int,int);
double combinat(int,int);
double genereMutant(int,int);

/* miscellaneous */

void print_hcode(int, unsigned int);

/* global job ID tracking for error reporting */
const char *getJobID();
void setJobID(const char *new_job_id);

/* bias tabel for sequences w/o secondary structures */
/*** mutational bias ***/

inline double nucleotide_bias(int position, int index);
double sequence_bias(int lefti, int righti, int nmut);
int partition_mutations(int i_ext_left, int i_int_left, int i_int_right, int i_ext_right, int nb_mut_unk_loop);
int partition_mutations_uniform(int i_ext_left, int i_int_left, int i_int_right, int i_ext_right, int nb_mut_unk_loop);
