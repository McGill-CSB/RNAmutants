/* global variables */

int rna_len = 0; /* length of the arn*/
char *input_tape;
int max_mutations = 0; /* has to be set to number of mutations + 1 */

/* structural constraint */

const int __hairpin_min_length__ = 3;  /* the minimal length of a hairpin loop */
const int __hairpin_max_length__ = 60; /* the minimal length of a hairpin loop */
const int __vis_range__ = 4;           /* range of visibility */
const int __max_size_bulge__ = 30;     /* maximal size of bulge in loop */
const int max_sequence_size = 10000;

/* global variables */

struct Zcell {
  double pf;
  //double pfGC;
  double mfe;
};

struct Zcell *****Zms;
struct Zcell *****Zes;
struct Zcell *****Zs;
struct Zcell *****Zm;
struct Zcell *****Ze;
struct Zcell *****ZdangN;
struct Zcell *****ZdangB;
struct Zcell *****ZdangL;
struct Zcell *****ZdangR;

/** table size values **/

int nb_value_triloop_table;
double *triloop_cmpt_table;
double *triloop_weight_table;
int nb_value_tetraloop_table;
double *tetraloop_cmpt_table;
double *tetraloop_weight_table;

/**** functions to declare ****/

int build_special_loop_list_for_given_config(unsigned int *list_of_config, int tabindex, int ileft, int iright, int nt_loop_left, int nt_loop_right, int nb_mutations,int *config);

/*********************/

#define INDEX_A 0
#define INDEX_C 1
#define INDEX_G 2
#define INDEX_U 3
#define INDEX_GHOST 4
#define INDEX_STOP -1

double mutbias[4][4] = {{1.0, 1.0, 1.0, 1.0},{1.0, 1.0, 1.0, 1.0},{1.0, 1.0, 1.0, 1.0},{1.0, 1.0, 1.0, 1.0}};

int uncst_mutation[2][4][4] = {
{{INDEX_A,INDEX_STOP},{INDEX_C,INDEX_STOP},{INDEX_G,INDEX_STOP},{INDEX_U,INDEX_STOP}},
{{INDEX_C,INDEX_G,INDEX_U,INDEX_STOP},{INDEX_A,INDEX_G,INDEX_U,INDEX_STOP},
{INDEX_A,INDEX_C,INDEX_U,INDEX_STOP},{INDEX_A,INDEX_C,INDEX_G,INDEX_STOP}}
};

/* 0: mutations, 1: right mutation, 2: left mutation, 3: mutation on both sides */
int bp_mutation[4][4][4][5][2] =
  /* A-A  */
  {{{{{INDEX_STOP,INDEX_STOP}},
     /* A-C  */
     {{INDEX_STOP,INDEX_STOP}},
     /* A-G  */
     {{INDEX_STOP,INDEX_STOP}},
     /* A-U  */
     {{INDEX_A,INDEX_U},{INDEX_STOP,INDEX_STOP}}},
    /* C-A  */
    {{{INDEX_STOP,INDEX_STOP}},
     /* C-C  */
     {{INDEX_STOP,INDEX_STOP}},
     /* C-G  */
     {{INDEX_C,INDEX_G},{INDEX_STOP,INDEX_STOP}},
     /* C-U  */
     {{INDEX_STOP,INDEX_STOP}}},
    /* G-A  */
    {{{INDEX_STOP,INDEX_STOP}},
     /* G-C  */
     {{INDEX_G,INDEX_C},{INDEX_STOP,INDEX_STOP}},
     /* G-G  */
     {{INDEX_STOP,INDEX_STOP}},
     /* G-U  */
     {{INDEX_G,INDEX_U},{INDEX_STOP,INDEX_STOP}}},
    /* U-A  */
    {{{INDEX_U,INDEX_A},{INDEX_STOP,INDEX_STOP}},
     /* U-C  */
     {{INDEX_STOP,INDEX_STOP}},
     /* U-G  */
     {{INDEX_U,INDEX_G},{INDEX_STOP,INDEX_STOP}},
     /* U-U  */
     {{INDEX_STOP,INDEX_STOP}}}},
   /*********** 1 mutation on the right *************/
  /* A-A with 1 mutation */
   {{{{INDEX_A,INDEX_U},{INDEX_STOP,INDEX_STOP}},
     /* A-C with 1 mutation */
     {{INDEX_A,INDEX_U},{INDEX_STOP,INDEX_STOP}},
     /* A-G with 1 mutation */
     {{INDEX_A,INDEX_U},{INDEX_STOP,INDEX_STOP}},
     /* A-U with 1 mutation */
     {{INDEX_STOP,INDEX_STOP}}},
    /* C-A with 1 mutation */
    {{{INDEX_C,INDEX_G},{INDEX_STOP,INDEX_STOP}},
     /* C-C with 1 mutation */
     {{INDEX_C,INDEX_G},{INDEX_STOP,INDEX_STOP}},
     /* C-G with 1 mutation */
     {{INDEX_STOP,INDEX_STOP}},
     /* C-U with 1 mutation */
     {{INDEX_C,INDEX_G},{INDEX_STOP,INDEX_STOP}}},
    /* G-A with 1 mutation */
    {{{INDEX_G,INDEX_C},{INDEX_G,INDEX_U},{INDEX_STOP,INDEX_STOP}},
     /* G-C with 1 mutation */
     {{INDEX_G,INDEX_U},{INDEX_STOP,INDEX_STOP}},
     /* G-G with 1 mutation */
     {{INDEX_G,INDEX_C},{INDEX_G,INDEX_U},{INDEX_STOP,INDEX_STOP}},
     /* G-U with 1 mutation */
     {{INDEX_G,INDEX_C},{INDEX_STOP,INDEX_STOP}}},
    /* U-A with 1 mutation */
    {{{INDEX_U,INDEX_G},{INDEX_STOP,INDEX_STOP}},
     /* U-C with 1 mutation */
     {{INDEX_U,INDEX_A},{INDEX_U,INDEX_G},{INDEX_STOP,INDEX_STOP}},
     /* U-G with 1 mutation */
     {{INDEX_U,INDEX_A},{INDEX_STOP,INDEX_STOP}},
     /* U-U with 1 mutation */
     {{INDEX_U,INDEX_A},{INDEX_U,INDEX_G},{INDEX_STOP,INDEX_STOP}}}},
   /*********** 1 mutation on left *************/
  /* A-A with 1 mutation */
   {{{{INDEX_U,INDEX_A},{INDEX_STOP,INDEX_STOP}},
     /* A-C with 1 mutation */
     {{INDEX_G,INDEX_C},{INDEX_STOP,INDEX_STOP}},
     /* A-G with 1 mutation */
     {{INDEX_C,INDEX_G},{INDEX_U,INDEX_G},{INDEX_STOP,INDEX_STOP}},
     /* A-U with 1 mutation */
     {{INDEX_G,INDEX_U},{INDEX_STOP,INDEX_STOP}}},
    /* C-A with 1 mutation */
    {{{INDEX_U,INDEX_A},{INDEX_STOP,INDEX_STOP}},
     /* C-C with 1 mutation */
     {{INDEX_G,INDEX_C},{INDEX_STOP,INDEX_STOP}},
     /* C-G with 1 mutation */
     {{INDEX_U,INDEX_G},{INDEX_STOP,INDEX_STOP}},
     /* C-U with 1 mutation */
     {{INDEX_A,INDEX_U},{INDEX_G,INDEX_U},{INDEX_STOP,INDEX_STOP}}},
    /* G-A with 1 mutation */
    {{{INDEX_U,INDEX_A},{INDEX_STOP,INDEX_STOP}},
     /* G-C with 1 mutation */
     {{INDEX_STOP,INDEX_STOP}},
     /* G-G with 1 mutation */
     {{INDEX_C,INDEX_G},{INDEX_U,INDEX_G},{INDEX_STOP,INDEX_STOP}},
     /* G-U with 1 mutation */
     {{INDEX_A,INDEX_U},{INDEX_STOP,INDEX_STOP}}},
    /* U-A with 1 mutation */
    {{{INDEX_STOP,INDEX_STOP}},
     /* U-C with 1 mutation */
     {{INDEX_G,INDEX_C},{INDEX_STOP,INDEX_STOP}},
     /* U-G with 1 mutation */
     {{INDEX_C,INDEX_G},{INDEX_STOP,INDEX_STOP}},
     /* U-U with 1 mutation */
     {{INDEX_A,INDEX_U},{INDEX_G,INDEX_U},{INDEX_STOP,INDEX_STOP}}}},
   /*********** 2 mutations *************/
   /* A-A with 2 mutations */
   {{{{INDEX_C,INDEX_G},{INDEX_G,INDEX_C},{INDEX_G,INDEX_U},{INDEX_U,INDEX_G},{INDEX_STOP,INDEX_STOP}},
     /* A-C with 2 mutations */
     {{INDEX_C,INDEX_G},{INDEX_G,INDEX_U},{INDEX_U,INDEX_A},{INDEX_U,INDEX_G},{INDEX_STOP,INDEX_STOP}},
     /* A-G with 2 mutations */
     {{INDEX_G,INDEX_C},{INDEX_G,INDEX_U},{INDEX_U,INDEX_A},{INDEX_STOP,INDEX_STOP}},
     /* A-U with 2 mutations */
     {{INDEX_C,INDEX_G},{INDEX_G,INDEX_C},{INDEX_U,INDEX_A},{INDEX_U,INDEX_G},{INDEX_STOP,INDEX_STOP}}},
    /* C-A with 2 mutations */
    {{{INDEX_A,INDEX_U},{INDEX_G,INDEX_C},{INDEX_G,INDEX_U},{INDEX_U,INDEX_G},{INDEX_STOP,INDEX_STOP}},
     /* C-C with 2 mutations */
     {{INDEX_A,INDEX_U},{INDEX_G,INDEX_U},{INDEX_U,INDEX_A},{INDEX_U,INDEX_G},{INDEX_STOP,INDEX_STOP}},
     /* C-G with 2 mutations */
     {{INDEX_A,INDEX_U},{INDEX_G,INDEX_C},{INDEX_G,INDEX_U},{INDEX_U,INDEX_A},{INDEX_STOP,INDEX_STOP}},
     /* C-U with 2 mutations */
     {{INDEX_G,INDEX_C},{INDEX_U,INDEX_A},{INDEX_U,INDEX_G},{INDEX_STOP,INDEX_STOP}}},
    /* G-A with 2 mutations */
    {{{INDEX_A,INDEX_U},{INDEX_C,INDEX_G},{INDEX_U,INDEX_G},{INDEX_STOP,INDEX_STOP}},
     /* G-C with 2 mutations */
     {{INDEX_A,INDEX_U},{INDEX_C,INDEX_G},{INDEX_U,INDEX_A},{INDEX_U,INDEX_G},{INDEX_STOP,INDEX_STOP}},
     /* G-G with 2 mutations */
     {{INDEX_A,INDEX_U},{INDEX_U,INDEX_A},{INDEX_STOP,INDEX_STOP}},
     /* G-U with 2 mutations */
     {{INDEX_C,INDEX_G},{INDEX_U,INDEX_A},{INDEX_U,INDEX_G},{INDEX_STOP,INDEX_STOP}}},
    /* U-A with 2 mutations */
    {{{INDEX_A,INDEX_U},{INDEX_C,INDEX_G},{INDEX_G,INDEX_U},{INDEX_G,INDEX_C},{INDEX_STOP,INDEX_STOP}},
     /* U-C with 2 mutations */
     {{INDEX_A,INDEX_U},{INDEX_C,INDEX_G},{INDEX_G,INDEX_U},{INDEX_STOP,INDEX_STOP}},
     /* U-G with 2 mutations */
     {{INDEX_A,INDEX_U},{INDEX_G,INDEX_C},{INDEX_G,INDEX_U},{INDEX_STOP,INDEX_STOP}},
     /* U-U with 2 mutations */
     {{INDEX_C,INDEX_G},{INDEX_G,INDEX_C},{INDEX_STOP,INDEX_STOP}}}}};
    
