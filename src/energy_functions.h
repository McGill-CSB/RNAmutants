
/* physical values */

#define DEFAULT_TEMPERATURE 310.0 /* 37 Celsius degres (273.15 + 37 = 310.15)*/

#define GAS_CONSTANT 1.9872
#define JOULE_GAS_CONSTANT 8.314472 

double penalty_close_helix(int,int);
double boltzmann_penalty_close_helix(int,int);
double ggg_counter(int, int, int, int, int, int, int);
int triloop_counter(int, int, int, int, int, int, int);
int tetraloop_counter(int, int, int, int, int, int, int);
double boltzmannStack(int, int, int, int);
double boltzmannBulge(int, int, int, int, int);
double boltzmannHairpin(int, int, int, int, int);
double boltzmannInternal_generic(int, int, int, int, int, int, int, int, int, int);
double boltzmannInternal_1x1(int, int, int, int, int, int);
double boltzmannInternal_1x2(int, int, int, int, int, int, int);
double boltzmannInternal_2x1(int, int, int, int, int, int, int);
double boltzmannInternal_2x2(int, int, int, int, int, int, int, int);
double boltzmannAddStemInMultiLoop(int);
double boltzmannExtendMultiLoop(int);
double boltzmannEndMultiLoop();
/*****/
double EStack(int, int, int, int);
double EBulge(int, int, int, int, int);
double EHairpin(int, int, int, int, int);
double EInternal_generic(int, int, int, int, int, int, int, int, int, int);
double EInternal_1x1(int, int, int, int, int, int);
double EInternal_1x2(int, int, int, int, int, int, int);
double EInternal_2x1(int, int, int, int, int, int, int);
double EInternal_2x2(int, int, int, int, int, int, int, int);
double EAddStemInMultiLoop(int);
double EExtendMultiLoop(int);
double EEndMultiLoop();
/*****/
double ggg_bonus(int, int, int, int, int, int, int);
double triloop_bonus(int, int, int, int, int, int, int);
double tetraloop_bonus(int, int, int, int, int, int, int);
