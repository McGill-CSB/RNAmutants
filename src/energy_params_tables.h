
/* dangle array store the data of dangle file of mfold */
/* index 1: position of the nt 0:up, 1:down            */
/* index 2: upper nt index  (X)                        */
/* index 3: lower nt index  (Y)                        */
/* index 3:  nt index (Z)                        */
/* example :                   XZ                      */
/*                             Y                       */

 double dangle[2][4][4][4];

/* sint2  array store the data of sint2 file of mfold  */
/* sint2[i][j][k][l][x][y]                             */
/*                                    X                */
/*                                   I K               */
/*                                   J L               */
/*                                    Y                */

 double sint2[4][4][4][4][4][4];

/* asint1x2 array store the data of asint1x2 file of mfold  */
/* asint1x2[i][j][k][l][x][y][z]                            */
/*                                5' ===> 3'                */
/*                                    X                     */
/*                                   I  K                   */
/*                                   J  L                   */
/*                                    YZ                    */
/*                                5' ===> 3'                */

 double asint1x2[4][4][4][4][4][4][4];

/* sint4 array store the data of sint4 file of mfold   */
/* sint4[i][j][k][l][u][v][x][y]                       */
/*                                    UV               */
/*                                   I  K              */
/*                                   J  L              */
/*                                    XY               */

 double sint4[4][4][4][4][4][4][4][4];

/* dangle array store the data of tstacki file of mfold  */
/* tstacki[i][j][k][l]                                   */
/*                                5' ===> 3'             */
/*                                   I K                 */
/*                                   J L                 */
/*                               3' <=== 5'              */

 double tstacki[4][4][4][4];

/* dangle array store the data of tstackh file of mfold  */
/* tstackh[i][j][k][l]                                   */
/*                                5' ===> 3'             */
/*                                   I K                 */
/*                                   J L                 */
/*                               3' <=== 5'              */

 double tstackh[4][4][4][4];

/* dangle array store the data of stack file of mfold  */
/* stack[i][j][k][l]                                   */
/*                                5' ===> 3'           */
/*                                   I K               */
/*                                   J L               */
/*                               3' <=== 5'            */

 double stack[4][4][4][4];

/* arrays store the data of loop file of mfold  */
/* index are for the length of the loop         */
/* size is not converted. Thus 0 index is empty */

double bulge[31];
double internal[31];
double hairpin[31];

/* listTetraLoop store the data of tloop file of mfold  */

typedef struct tetraLoop{
  char motif[6];
  int hcode; /* int representation */
  double bonus;
  struct tetraLoop *next;
} *TetraLoop;

TetraLoop listTetraLoop;
double refTableTetraLoop[4096]; /* this array is used to have a direct acces to the values of listTretraLoop */

/* listTriLoop store the data of tloop file of mfold  */

typedef struct triLoop{
  char motif[5];
  int hcode; /* int representation */
  double bonus;
  struct triLoop *next;
} *TriLoop;

TriLoop listTriLoop;
double refTableTriLoop[1024]; /* this array is used to have a direct acces to the values of listTretraLoop */

/* differents parameters used for loops */

/* for all */
 double param_large_loop;

/* internal loop*/
 double max_correction_asym_internal_loop;
 double array_ninio_internal_loop[5];

/* multi-loop */
 double offset_multi_loop;
 double free_base_penalty_multi_loop;
 double helix_penalty_multi_loop;

/* efn2 multi-loop */
 double offset_efn2_multi_loop;
 double free_base_penalty_efn2_multi_loop;
 double helix_penalty_efn2_multi_loop;

/* other */
 double au_penalty;
 double ggg_hairpin_bonus;
 double c_hairpin_slope;
 double c_hairpin_intersect;
 double c_hairpin_3;
 double inter_molecular_init_energy;
 int gail_rule;

/* functions used */

void readDangle(const char *filename);
void readSint2(const char *filename);
void readSint4(const char *filename);
void readAsint1x2(const char *filename);
void readTstacki(const char *filename);
void readTstackh(const char *filename);
void readStack(const char *filename);
void readLoop(const char *filename);
int readTetraLoop(const char *filename);
int readTriLoop(const char *filename);
void readMiscLoop(const char *filename);
