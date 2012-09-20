#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<float.h> /* for DBL_MAX */
#include "reader_energy_params.h"
#include "util.h"

#define INDEX_A 0
#define INDEX_C 1
#define INDEX_G 2
#define INDEX_U 3

const int __size_buffer__=200;

/* dangle array store the data of dangle file of mfold */
/* index 1: position of the nt 0:up, 1:down            */
/* index 2: upper nt index  (X)                        */
/* index 3: lower nt index  (Y)                        */
/* index 3: extern nt index (Z)                        */
/*                              5' ==> 3'              */
/*                                 XZ                  */
/*                                 Y                   */
/*                             5' <== 3'               */

extern double dangle[2][4][4][4];

/* sint2  array store the data of sint2 file of mfold  */
/* sint2[i][j][k][l][x][y]                             */
/*                                5' ==> 3'            */
/*                                    X                */
/*                                   I K               */
/*                                   J L               */
/*                                    Y                */
/*                                5' <== 3'            */

extern double sint2[4][4][4][4][4][4];

/* asint1x2 array store the data of asint1x2 file of mfold  */
/* asint1x2[i][j][k][l][x][y][z]                            */
/*                                5' ===> 3'                */
/*                                    X                     */
/*                                   I  K                   */
/*                                   J  L                   */
/*                                    YZ                    */
/*                                5' <=== 3'                */

extern double asint1x2[4][4][4][4][4][4][4];

/* sint4 array store the data of sint4 file of mfold   */
/* sint4[i][j][k][l][u][v][x][y]                       */
/*                                5' ===> 3'           */
/*                                    UV               */
/*                                   I  K              */
/*                                   J  L              */
/*                                    YX               */
/*                                5' <=== 3'           */

extern double sint4[4][4][4][4][4][4][4][4];

/* dangle array store the data of tstacki file of mfold  */
/* tstacki[i][j][k][l]                                   */
/*                                5' ===> 3'             */
/*                                   I K                 */
/*                                   J L                 */
/*                               3' <=== 5'              */

extern double tstacki[4][4][4][4];

/* dangle array store the data of tstackh file of mfold  */
/* tstackh[i][j][k][l]                                   */
/*                                5' ===> 3'             */
/*                                   I K                 */
/*                                   J L                 */
/*                               3' <=== 5'              */

extern double tstackh[4][4][4][4];

/* dangle array store the data of stack file of mfold  */
/* stack[i][j][k][l]                                   */
/*                                5' ===> 3'           */
/*                                   I K               */
/*                                   J L               */
/*                               3' <=== 5'            */

extern double stack[4][4][4][4];

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

extern TetraLoop listTetraLoop;
extern double refTableTetraLoop[444445]; /* this array is used to have a direct acces to the values of listTretraLoop */

/* listTriLoop store the data of tloop file of mfold  */

typedef struct triLoop{
	char motif[5];
	int hcode; /* int representation */
	double bonus;
	struct triLoop *next;
} *TriLoop;

extern TriLoop listTriLoop;
extern double refTableTriLoop[44445]; /* this array is used to have a direct acces to the values of listTetraLoop */

/* differents parameters used for loops */

/* for all */
extern double param_large_loop;

/* internal loop*/
extern double max_correction_asym_internal_loop;
extern double array_ninio_internal_loop[5];

/* multi-loop */
extern double offset_multi_loop;
extern double free_base_penalty_multi_loop;
extern double helix_penalty_multi_loop;

/* efn2 multi-loop */
extern double offset_efn2_multi_loop;
extern double free_base_penalty_efn2_multi_loop;
extern double helix_penalty_efn2_multi_loop;

/* other */
extern double au_penalty;
extern double ggg_hairpin_bonus;
extern double c_hairpin_slope;
extern double c_hairpin_intersect;
extern double c_hairpin_3;
extern double inter_molecular_init_energy;
extern int gail_rule;


/*********************************************************************************************/
/*                                                                                           */
/* FUNCTIONS                                                                                 */
/*                                                                                           */
/*********************************************************************************************/


/*********************************************************************************************/
/*  read dangle.bt file                                                                       */
/*********************************************************************************************/

void readDangle(const char *filename) {
	
	FILE *f_data;
	char vs[4][50]; /* before evaluation, the fields must be checked (case of ".") */
	int nn,xx,yy,zz;
	
	if (!(f_data=fopen(filename,"r"))) {
		fprintf(stderr,"Cannot open file %s%s\n",filename,getJobID());
		exit(EXIT_FAILURE);
	}
	
	for (nn=0;nn<2;nn++) {
		for (xx=0;xx<4;xx++) {
			for (yy=0;yy<4;yy++) {
				for (zz=0;zz<4;zz++) {
					dangle[nn][xx][yy][zz] = -DBL_MAX;
				}
			}
		}
	}
	
	for (nn=0;nn<2;nn++) {
		for (yy=0;yy<4;yy++) {
			for (xx=0;xx<4;xx++) {
				if (fscanf(f_data,"%s%s%s%s",
						   &vs[0][0],&vs[1][0],&vs[2][0],&vs[3][0]) != 4)
				{ /* read failed */
					fprintf(stderr,"dangle: Corrupted file at [n=%d,x=%d,y=%d]%s\n",nn,xx,yy,getJobID());
					exit(EXIT_FAILURE);
				}
				for (zz=0;zz<4;zz++) {
					if (vs[zz][0] == '.') {
						dangle[nn][xx][yy][zz] = DBL_MAX;
					}
					else {
						dangle[nn][xx][yy][zz] = atof(&vs[zz][0]);	
					}
				}
			}   
		}
	}
	
	{
		char buff[20];
		if (fscanf(f_data,"%s",&buff[0])==1) {
			fprintf(stderr,"dangle: Corrupted file (read %s at the expected EOF)%s\n",buff,getJobID());
			exit(EXIT_FAILURE);
		}
	}
	
	fclose(f_data);
}


/*********************************************************************************************/
/*  read sint2.bt file                                                                       */
/*********************************************************************************************/

void readSint2(const char *filename) {
	
	FILE *f_data;
	char vs[4][20]; /* before evaluation, the fields must be checked (case of ".") */
	int ii,jj,kk,ll,xx,yy,indI,indJ,indK,indL,indX,indY;
	
	if (!(f_data=fopen(filename,"r"))) {
		fprintf(stderr,"Cannot open file %s%s\n",filename,getJobID());
		exit(EXIT_FAILURE);
	}
	
	
	/* init table */
	
	for(ii=0;ii<4;ii++) {
		for(jj=0;jj<4;jj++) {
			for(kk=0;kk<4;kk++) {
				for(ll=0;ll<4;ll++) {
					for(xx=0;xx<4;xx++) {
						for(yy=0;yy<4;yy++) {
							sint2[ii][jj][kk][ll][xx][yy]=-DBL_MAX; /* signal for an error in the main prog: should never be reach */
						}
					}
				}
			}
		}
	}
	
	/* read and store data */
	
	for (ii=0;ii<6;ii++) {
		
		switch(ii) { /* define I, J index */
			case 0:
				indI = INDEX_A;
				indJ = INDEX_U;
				break;
			case 1:
				indI = INDEX_C;
				indJ = INDEX_G;
				break;
			case 2:
				indI = INDEX_G;
				indJ = INDEX_C;
				break;
			case 3:
				indI = INDEX_U;
				indJ = INDEX_A;
				break;
			case 4:
				indI = INDEX_G;
				indJ = INDEX_U;
				break;
			case 5:
				indI = INDEX_U;
				indJ = INDEX_G;
				break;
			default:
				fprintf(stderr,"sint2: Corrupted file%s\n",getJobID());
				exit(EXIT_FAILURE);
		}
		
		for (indX=0;indX<4;indX++) { /* enumerate each line */
			
			for (jj=0;jj<6;jj++) { /* cut each line in 6 fields (closing K,L basepairs) */
				
				switch(jj) { /* define K, L index */
					case 0:
						indK = INDEX_A;
						indL = INDEX_U;
						break;
					case 1:
						indK = INDEX_C;
						indL = INDEX_G;
						break;
					case 2:
						indK = INDEX_G;
						indL = INDEX_C;
						break;
					case 3:
						indK = INDEX_U;
						indL = INDEX_A;
						break;
					case 4:
						indK = INDEX_G;
						indL = INDEX_U;
						break;
					case 5:
						indK = INDEX_U;
						indL = INDEX_G;
						break;
					default:
						fprintf(stderr,"sint2: Corrupted file%s\n",getJobID());
						exit(EXIT_FAILURE);
				}
				
				if (fscanf(f_data,"%s%s%s%s",&vs[0][0],&vs[1][0],&vs[2][0],&vs[3][0])!=4) {
					fprintf(stderr,"sint2: Corrupted file%s\n",getJobID());
					exit(EXIT_FAILURE); 
				}
				
				for(indY=0;indY<4;indY++) {
					if (vs[indY][0] == '.') {
						sint2[indI][indJ][indK][indL][indX][indY] = DBL_MAX;
					}
					else {
						sint2[indI][indJ][indK][indL][indX][indY] = atof(&vs[indY][0]);	
					}
				}
			}
		}
	}   
	
	
	{
		char buff[20];
		if (fscanf(f_data,"%s",&buff[0])==1) {
			fprintf(stderr,"sint2: Corrupted file (read %s at the expected EOF)%s\n",buff,getJobID());
			exit(EXIT_FAILURE);
		}
	}
	
	fclose(f_data);
}

/*********************************************************************************************/
/*  read asint1x2.bt file                                                                       */
/*********************************************************************************************/

void readAsint1x2(const char *filename) {
	
	FILE *f_data;
	char vs[4][20]; /* before evaluation, the fields must be checked (case of ".") */
	int ii,jj,kk,ll,xx,yy,zz,indI,indJ,indK,indL,indX,indY,indZ;
	
	
	if (!(f_data=fopen(filename,"r"))) {
		fprintf(stderr,"Cannot open file %s%s\n",filename,getJobID());
		exit(EXIT_FAILURE);
	}
	
	
	/* init table */
	
	for(ii=0;ii<4;ii++) {
		for(jj=0;jj<4;jj++) {
			for(kk=0;kk<4;kk++) {
				for(ll=0;ll<4;ll++) {
					for(xx=0;xx<4;xx++) {
						for(yy=0;yy<4;yy++) {
							for(zz=0;zz<4;zz++) {
								asint1x2[ii][jj][kk][ll][xx][yy][zz]=-DBL_MAX; /* signal for an error in the main prog: should never be reach */
							}
						}
					}
				}
			}
		}
	}
	
	/* read and store data */
	
	for (ii=0;ii<6;ii++) {
		
		switch(ii) { /* define I, J index */
			case 0:
				indI = INDEX_A;
				indJ = INDEX_U;
				break;
			case 1:
				indI = INDEX_C;
				indJ = INDEX_G;
				break;
			case 2:
				indI = INDEX_G;
				indJ = INDEX_C;
				break;
			case 3:
				indI = INDEX_U;
				indJ = INDEX_A;
				break;
			case 4:
				indI = INDEX_G;
				indJ = INDEX_U;
				break;
			case 5:
				indI = INDEX_U;
				indJ = INDEX_G;
				break;
			default:
				fprintf(stderr,"asint1x2: Corrupted file%s\n",getJobID());
				exit(EXIT_FAILURE);
		}
		
		for (indZ=0;indZ<4;indZ++) { /* enumerate each line */
			
			for (indX=0;indX<4;indX++) { /* enumerate each line */
				
				for (jj=0;jj<6;jj++) { /* cut each line in 6 fields (closing K,L basepairs) */
					
					switch(jj) { /* define K, L index */
						case 0:
							indK = INDEX_A;
							indL = INDEX_U;
							break;
						case 1:
							indK = INDEX_C;
							indL = INDEX_G;
							break;
						case 2:
							indK = INDEX_G;
							indL = INDEX_C;
							break;
						case 3:
							indK = INDEX_U;
							indL = INDEX_A;
							break;
						case 4:
							indK = INDEX_G;
							indL = INDEX_U;
							break;
						case 5:
							indK = INDEX_U;
							indL = INDEX_G;
							break;
						default:
							fprintf(stderr,"asint1x2: Corrupted file%s\n",getJobID());
							exit(EXIT_FAILURE);
					}
					
					if (fscanf(f_data,"%s%s%s%s",&vs[0][0],&vs[1][0],&vs[2][0],&vs[3][0])!=4) {
						fprintf(stderr,"asint1x2: Corrupted file%s\n",getJobID());
						exit(EXIT_FAILURE); 
					}
					
					for(indY=0;indY<4;indY++) {
						if (vs[indY][0] == '.') {
							asint1x2[indI][indJ][indK][indL][indX][indY][indZ] = DBL_MAX;
						}
						else {
							asint1x2[indI][indJ][indK][indL][indX][indY][indZ] = atof(&vs[indY][0]);	
						}
					}
				}
			}
		}
	}
	
	
	{
		char buff[20];
		if (fscanf(f_data,"%s",&buff[0])==1) {
			fprintf(stderr,"asint1x2: Corrupted file (read %s at the expected EOF)%s\n",buff,getJobID());
			exit(EXIT_FAILURE);
		}
	}
	
	
	fclose(f_data);
}

/*********************************************************************************************/
/*  read sint4.bt file                                                                       */
/*********************************************************************************************/


void readSint4(const char *filename) {
	
	FILE *f_data;
	char vs[100]; /* before evaluation, the fields must be checked (case of ".") */
	int ii,jj,kk,ll,uu,vv,xx,yy,indI,indJ,indK,indL;
	
	if (!(f_data=fopen(filename,"r"))) {
		fprintf(stderr,"Cannot open file %s%s\n",filename,getJobID());
		exit(EXIT_FAILURE);
	}
	
	
	/* init table */
	
	for(ii=0;ii<4;ii++) {
		for(jj=0;jj<4;jj++) {
			for(kk=0;kk<4;kk++) {
				for(ll=0;ll<4;ll++) {
					for(uu=0;uu<4;uu++) {
						for(vv=0;vv<4;vv++) {
							for(xx=0;xx<4;xx++) {
								for(yy=0;yy<4;yy++) {
									sint4[ii][jj][kk][ll][uu][vv][xx][yy]=-DBL_MAX; /* flag for an error in the main prog: -1 should never be used */
								}
							}
						}
					}
				}
			}
		}
	}
	
	/* read and store data */
	
	for (ii=0;ii<6;ii++) {
		
		switch(ii) { /* enumerate the left basepair */
			case 0:
				indI = INDEX_A;
				indJ = INDEX_U;
				break;
			case 1:
				indI = INDEX_C;
				indJ = INDEX_G;
				break;
			case 2:
				indI = INDEX_G;
				indJ = INDEX_C;
				break;
			case 3:
				indI = INDEX_U;
				indJ = INDEX_A;
				break;
			case 4:
				indI = INDEX_G;
				indJ = INDEX_U;
				break;
			case 5:
				indI = INDEX_U;
				indJ = INDEX_G;
				break;
			default:
				fprintf(stderr,"sint4: Corrupted file%s\n",getJobID());
				exit(EXIT_FAILURE);
		}
		
		for (jj=0;jj<6;jj++) { /* enumerate the right basepair */
			
			switch(jj) { /* define K, L index */
				case 0:
					indK = INDEX_A;
					indL = INDEX_U;
					break;
				case 1:
					indK = INDEX_C;
					indL = INDEX_G;
					break;
				case 2:
					indK = INDEX_G;
					indL = INDEX_C;
					break;
				case 3:
					indK = INDEX_U;
					indL = INDEX_A;
					break;
				case 4:
					indK = INDEX_G;
					indL = INDEX_U;
					break;
				case 5:
					indK = INDEX_U;
					indL = INDEX_G;
					break;
				default:
					fprintf(stderr,"sint4: Corrupted file%s\n",getJobID());
					exit(EXIT_FAILURE);
			}
			
			/* now, we can store the block */
			
			for(uu=0;uu<4;uu++){
				for(vv=0;vv<4;vv++){ /* enumerate lines of the block */
					for(xx=0;xx<4;xx++) { /* cut each line in 4 blocks */
						for(yy=0;yy<4;yy++) {
							
							if (fscanf(f_data,"%s",&vs[0])!=1) {
								fprintf(stderr,"sint4: Corrupted file%s\n",getJobID());
								exit(EXIT_FAILURE); 
							}
							
							if (vs[0] == '.') {
								sint4[indI][indJ][indK][indL][uu][xx][yy][vv] = DBL_MAX;
							}
							else {
								sint4[indI][indJ][indK][indL][uu][xx][yy][vv] = atof(&vs[0]);	
							}
						}
					}
				}
			}
		}
	}
	
	
	{
		char buff[20];
		if (fscanf(f_data,"%s",&buff[0])==1) {
			fprintf(stderr,"sint4: Corrupted file (read %s at the expected EOF)%s\n",buff,getJobID());
			exit(EXIT_FAILURE);
		}
	}
	
	fclose(f_data);
	
}

/*********************************************************************************************/
/*  read tstacki.bt file                                                                     */
/*********************************************************************************************/

void readTstacki(const char *filename) {
	
	FILE *f_data;
	char vs[4][20]; /* before evaluation, the fields must be checked (case of ".") */
	int ii,jj,kk,ll;
	
	if (!(f_data=fopen(filename,"r"))) {
		fprintf(stderr,"Cannot open file %s%s\n",filename,getJobID());
		exit(EXIT_FAILURE);
	}
	
	
	/* init table */
	
	for(ii=0;ii<4;ii++) {
		for(jj=0;jj<4;jj++) {
			for(kk=0;kk<4;kk++) {
				for(ll=0;ll<4;ll++) {
					tstacki[ii][jj][kk][ll]=-DBL_MAX; /* signal for an error in the main prog: -1 should never be used */
				}
			}
		}
	}
	
	/* read and store data */
	
	for (ii=0;ii<4;ii++) { /* enumerate each block lines */
		for (kk=0;kk<4;kk++) { /* enumerate each line of the blocks */
			for (jj=0;jj<4;jj++) { /* cut each line in 4 fields */
				if (fscanf(f_data,"%s%s%s%s",&vs[0][0],&vs[1][0],&vs[2][0],&vs[3][0])!=4) {
					fprintf(stderr,"tstacki: Corrupted file%s\n",getJobID());
					exit(EXIT_FAILURE); 
				}
				
				for(ll=0;ll<4;ll++) {
					if (vs[ll][0] == '.') {
						tstacki[ii][jj][kk][ll] = DBL_MAX;
					}
					else {
						tstacki[ii][jj][kk][ll] = atof(&vs[ll][0]);	
					}
				}
			}
		}
	}   
	
	
	{
		char buff[20];
		if (fscanf(f_data,"%s",&buff[0])==1) {
			fprintf(stderr,"tstacki: Corrupted file (read %s at the expected EOF)%s\n",buff,getJobID());
			exit(EXIT_FAILURE);
		}
	}
	
	fclose(f_data);
}

/*********************************************************************************************/
/*  read tstackh.bt file                                                                     */
/*********************************************************************************************/

void readTstackh(const char *filename) {
	
	FILE *f_data;
	char vs[4][20]; /* before evaluation, the fields must be checked (case of ".") */
	int ii,jj,kk,ll;
	
	if (!(f_data=fopen(filename,"r"))) {
		fprintf(stderr,"Cannot open file %s%s\n",filename,getJobID());
		exit(EXIT_FAILURE);
	}
	
	
	/* init table */
	
	for(ii=0;ii<4;ii++) {
		for(jj=0;jj<4;jj++) {
			for(kk=0;kk<4;kk++) {
				for(ll=0;ll<4;ll++) {
					tstackh[ii][jj][kk][ll]=-DBL_MAX; /* signal for an error in the main prog: -1 should never be used */
				}
			}
		}
	}
	
	/* read and store data */
	
	for (ii=0;ii<4;ii++) { /* enumerate each block lines */
		for (kk=0;kk<4;kk++) { /* enumerate each line of the blocks */
			for (jj=0;jj<4;jj++) { /* cut each line in 4 fields */
				if (fscanf(f_data,"%s%s%s%s",&vs[0][0],&vs[1][0],&vs[2][0],&vs[3][0])!=4) {
					fprintf(stderr,"tstackh: Corrupted file%s\n",getJobID());
					exit(EXIT_FAILURE); 
				}
				
				for(ll=0;ll<4;ll++) {
					if (vs[ll][0] == '.') {
						tstackh[ii][jj][kk][ll] = DBL_MAX;
					}
					else {
						tstackh[ii][jj][kk][ll] = atof(&vs[ll][0]);	
					}
				}
			}
		}
	}   
	
	
	{
		char buff[20];
		if (fscanf(f_data,"%s",&buff[0])==1) {
			fprintf(stderr,"tstackh: Corrupted file (read %s at the expected EOF)%s\n",buff,getJobID());
			exit(EXIT_FAILURE);
		}
	}
	
	fclose(f_data);
}


/*********************************************************************************************/
/*  read stack.bt file                                                                       */
/*********************************************************************************************/

void readStack(const char *filename) {
	
	FILE *f_data;
	char vs[4][20]; /* before evaluation, the fields must be checked (case of ".") */
	int ii,jj,kk,ll;
	
	if (!(f_data=fopen(filename,"r"))) {
		fprintf(stderr,"Cannot open file %s%s\n",filename,getJobID());
		exit(EXIT_FAILURE);
	}
	
	
	/* init table */
	
	for(ii=0;ii<4;ii++) {
		for(jj=0;jj<4;jj++) {
			for(kk=0;kk<4;kk++) {
				for(ll=0;ll<4;ll++) {
					stack[ii][jj][kk][ll]=-DBL_MAX; /* signal for an error in the main prog: -1 should never be used */
				}
			}
		}
	}
	
	/* read and store data */
	
	for (ii=0;ii<4;ii++) { /* enumerate each block lines */
		for (kk=0;kk<4;kk++) { /* enumerate each line of the blocks */
			for (jj=0;jj<4;jj++) { /* cut each line in 4 fields */
				if (fscanf(f_data,"%s%s%s%s",&vs[0][0],&vs[1][0],&vs[2][0],&vs[3][0])!=4) {
					fprintf(stderr,"stack: Corrupted file%s\n",getJobID());
					exit(EXIT_FAILURE); 
				}
				
				for(ll=0;ll<4;ll++) {
					if (vs[ll][0] == '.') {
						stack[ii][jj][kk][ll] = DBL_MAX;
					}
					else {
						stack[ii][jj][kk][ll] = atof(&vs[ll][0]);	
					}
				}
			}
		}
	}   
	
	
	{
		char buff[20];
		if (fscanf(f_data,"%s",&buff[0])==1) {
			fprintf(stderr,"stack: Corrupted file (read %s at the expected EOF)%s\n",buff,getJobID());
			exit(EXIT_FAILURE);
		}
	}
	
	fclose(f_data);
}

/*********************************************************************************************/
/*  read loop.bt file                                                                        */
/*********************************************************************************************/

void readLoop(const char *filename) {
	
	FILE *f_data;
	char *buffer;
	char vs[3][20]; /* before evaluation, the fields must be checked (case of ".") */
	int ii,cmptLine,size;
	
	if (!(buffer=(char *)malloc(__size_buffer__ * sizeof(char)))) {
		fprintf(stderr,"Cannot allocate memory to buffer%s\n",getJobID());
		exit(EXIT_FAILURE);
	}
	
	
	if (!(f_data=fopen(filename,"r"))) {
		fprintf(stderr,"Cannot open file %s%s\n",filename,getJobID());
		exit(EXIT_FAILURE);
	}
	
	/* init table in order to avoid double affectation */
	
	bulge[0] = DBL_MAX;
	internal[0] = DBL_MAX;
	hairpin[0] = DBL_MAX;
	
	for (ii=1;ii<31;ii++) {
		bulge[ii] = -1;
		internal[ii] = -1;
		hairpin[ii] = -1;
	}
	
	cmptLine=0;
	
	while(fgets(buffer,__size_buffer__,f_data)) {
		if (sscanf(buffer,"%d%s%s%s",&size,&vs[0][0],&vs[1][0],&vs[2][0]) != 4)
		{ /* read failed */
			fprintf(stderr,"loop: Corrupted file%s\n",getJobID());
			exit(EXIT_FAILURE);
		}
		
		/* check the first field (length) */
		
		if ((size<1)||(size>30)) { /* read fails */
			fprintf(stderr,"loop: Corrupted file (size)%s\n",getJobID());
			exit(EXIT_FAILURE);
		}
		
		/* second field is for internal loop */
		
		if (internal[size]!=-1) { /* read fails */
			fprintf(stderr,"loop: Corrupted file (internal loop)%s\n",getJobID());
			exit(EXIT_FAILURE);
		}
		
		if (vs[0][0] == '.') {
			internal[size] = DBL_MAX;
		}
		else {
			internal[size] = atof(&vs[0][0]);	
		}
		
		/* third field is for bulge */
		
		if (bulge[size]!=-1) { /* read fails */
			fprintf(stderr,"loop: Corrupted file (bulge)%s\n",getJobID());
			exit(EXIT_FAILURE);
		}
		
		if (vs[1][0] == '.') {
			bulge[size] = DBL_MAX;
		}
		else {
			bulge[size] = atof(&vs[1][0]);	
		}
		
		/* fourth field is for hairpin */
		
		if (hairpin[size]!=-1) { /* read fails */
			fprintf(stderr,"loop: Corrupted file (hairpin)%s\n",getJobID());
			exit(EXIT_FAILURE);
		}
		
		if (vs[2][0] == '.') {
			hairpin[size] = DBL_MAX;
		}
		else {
			hairpin[size] = atof(&vs[2][0]);	
		}
		
		cmptLine++;
	}
	
	if (cmptLine != 30) {
		fprintf(stderr,"loop: Corrupted file (data)%s\n",getJobID());
		exit(EXIT_FAILURE);
	}
	
	fclose(f_data);
	free(buffer);
}

/*********************************************************************************************/
/*  read tloop.bt file                                                                       */
/*********************************************************************************************/

int readTetraLoop(const char *filename) {
	
	FILE *f_data;
	char *buffer;
	char motif[20];
	int ii,cmptLine,tmp_hcode,nt_code;
	float fv;
	TetraLoop tmpCell = NULL;
	
	if (!(buffer=(char *)malloc(__size_buffer__ * sizeof(char)))) {
		fprintf(stderr,"Cannot allocate memory to buffer%s\n",getJobID());
		exit(EXIT_FAILURE);
	}
	
	if (!(f_data=fopen(filename,"r"))) {
		fprintf(stderr,"Cannot open file %s%s\n",filename,getJobID());
		exit(EXIT_FAILURE);
	}
	
	/* init table in order to avoid double affectation */
	
	for (ii=0;ii<4096;ii++) {
		refTableTetraLoop[ii] = 0;
	}
	
	cmptLine=0;
	
	while(fgets(buffer,__size_buffer__,f_data)) {
		if (sscanf(buffer,"%s%f",&motif[0],&fv) != 2)
		{ /* read failed */
			fprintf(stderr,"tloop: Corrupted file%s\n",getJobID());
			exit(EXIT_FAILURE);
		}
		
		/* check the first field */
		
		if (strlen(&motif[0])!=6) { /* read fails */
			fprintf(stderr,"tloop: Corrupted file (motif %s)%s\n",&motif[0],getJobID());
			exit(EXIT_FAILURE);
		}
		
		/* int tmpCell or create a new cell */
		
		if (tmpCell) {
			tmpCell->next = (TetraLoop)malloc(sizeof(struct tetraLoop));
			tmpCell = tmpCell->next;
			tmpCell->next = NULL;
		}
		else {
			listTetraLoop = tmpCell = (TetraLoop)malloc(sizeof(struct tetraLoop));
			tmpCell->next = NULL;
		}
		
		/* compute hcode */
		
		tmp_hcode=0;
		for (ii=0;ii<6;ii++) {
			switch (motif[ii]) {
				case 'A':
					nt_code=0;
					break;
				case 'C':
					nt_code=1;
					break;
				case 'G':
					nt_code=2;
					break;
				case 'U':
					nt_code=3;
					break;
				default:
					fprintf(stderr,"tloop: Unexpected character %c%s\n",motif[ii],getJobID());
					exit(EXIT_FAILURE);
			}
			tmp_hcode = (tmp_hcode<<2) | nt_code;
		}
		
		if (tmp_hcode>4095) {
			fprintf(stderr,"readTetraLoop: corrupted data (bad hcode)%s\n",getJobID());
			exit(EXIT_FAILURE);
		}
		
		/* fill cell */
		
		memcpy(tmpCell->motif,&motif[0],6*sizeof(char));
		tmpCell->hcode = tmp_hcode;
		tmpCell->bonus = (double)(fv);
		
		/* update quick acces table */
		
		if(refTableTetraLoop[tmp_hcode]) {
			fprintf(stderr,"readTetraLoop: corrupted data (entry conflict: %d)%s\n",tmp_hcode,getJobID());
			exit(EXIT_FAILURE);
		}
		
		refTableTetraLoop[tmp_hcode] = tmpCell->bonus;
		
		cmptLine++;
	}
	
	fclose(f_data);
	free(buffer);
	
	return cmptLine;
	
}

/*********************************************************************************************/
/*  read triloop.bt file                                                                     */
/*********************************************************************************************/

int readTriLoop(const char *filename) {
	
	FILE *f_data;
	char *buffer;
	char motif[20];
	int ii,cmptLine,tmp_hcode,nt_code;
	float fv;
	TriLoop tmpCell = NULL;
	
	if (!(buffer=(char *)malloc(__size_buffer__ * sizeof(char)))) {
		fprintf(stderr,"Cannot allocate memory to buffer%s\n",getJobID());
		exit(EXIT_FAILURE);
	}
	
	if (!(f_data=fopen(filename,"r"))) {
		fprintf(stderr,"Cannot open file %s%s\n",filename,getJobID());
		exit(EXIT_FAILURE);
	}
	
	/* init table in order to avoid double affectation */
	
	for (ii=0;ii<1024;ii++) {
		refTableTriLoop[ii] = 0;
	}
	
	cmptLine=0;
	
	while(fgets(buffer,__size_buffer__,f_data)) {
		if (sscanf(buffer,"%s%f",&motif[0],&fv) != 2)
		{ /* read failed */
			fprintf(stderr,"triloop: Corrupted file%s\n",getJobID());
			exit(EXIT_FAILURE);
		}
		
		/* check the first field */
		
		if (strlen(&motif[0])!=5) { /* read fails */
			fprintf(stderr,"triloop: Corrupted file (motif %s)%s\n",&motif[0],getJobID());
			exit(EXIT_FAILURE);
		}
		
		/* int tmpCell or create a new cell */
		
		if (tmpCell) {
			tmpCell->next = (TriLoop)malloc(sizeof(struct triLoop));
			tmpCell = tmpCell->next;
			tmpCell->next = NULL;
		}
		else {
			listTriLoop = tmpCell = (TriLoop)malloc(sizeof(struct triLoop));
			tmpCell->next = NULL;
		}
		
		/* compute hcode */
		
		tmp_hcode=0;
		for (ii=0;ii<5;ii++) {
			switch (motif[ii]) {
				case 'A':
					nt_code=0;
					break;
				case 'C':
					nt_code=1;
					break;
				case 'G':
					nt_code=2;
					break;
				case 'U':
					nt_code=3;
					break;
				default:
					fprintf(stderr,"triloop: Unexpected character %c%s\n",motif[ii],getJobID());
					exit(EXIT_FAILURE);
			}
			tmp_hcode = (tmp_hcode<<2) & nt_code;
		}
		
		if (tmp_hcode>1023) {
			fprintf(stderr,"readTriLoop: corrupted data (bad hcode)%s\n",getJobID());
			exit(EXIT_FAILURE);
		}
		
		/* fill cell */
		
		memcpy(tmpCell->motif,&motif[0],5*sizeof(char));
		tmpCell->hcode = tmp_hcode;
		tmpCell->bonus = (double)(fv);
		
		/* update quick acces table */
		
		if(refTableTriLoop[tmp_hcode]) {
			fprintf(stderr,"readTriLoop: corrupted data (entry conflict)%s\n",getJobID());
			exit(EXIT_FAILURE);
		}
		
		refTableTriLoop[tmp_hcode] = tmpCell->bonus;
		
		cmptLine++;
	}
	
	fclose(f_data);
	free(buffer);
	
	return cmptLine;
	
}

/*********************************************************************************************/
/*  read miscloop.bt file                                                                    */
/*********************************************************************************************/

void readMiscLoop(const char *filename) {
	
	FILE *f_data;
	float fv1,fv2,fv3,fv4;
	
	if (!(f_data=fopen(filename,"r"))) {
		fprintf(stderr,"Cannot open file %s%s\n",filename,getJobID());
		exit(EXIT_FAILURE);
	}
	
	/* parameter for general formula (zuker use 1.9872 for gas constant) */
	
	if (fscanf(f_data,"%f",&fv1) != 1) { /* read failed */
		fprintf(stderr,"miscloop: corrupted file (fgets 1)%s\n",getJobID());
		exit(EXIT_FAILURE);
	}
	param_large_loop = (double)(fv1);
	
	/* internal loop*/
	if (fscanf(f_data,"%f",&fv1) != 1) { /* read failed */
		fprintf(stderr,"miscloop: corrupted file (fgets 2)%s\n",getJobID());
		exit(EXIT_FAILURE);
	}
	max_correction_asym_internal_loop = (double)(fv1);
	
	
	if (fscanf(f_data,"%f%f%f%f",&fv1,&fv2,&fv3,&fv4) != 4) { /* read failed */
		fprintf(stderr,"miscloop: corrupted file (fgets 3)%s\n",getJobID());
		exit(EXIT_FAILURE);
	}
	array_ninio_internal_loop[0] = 0.0;
	array_ninio_internal_loop[1] = (double)(fv1);
	array_ninio_internal_loop[2] = (double)(fv2);
	array_ninio_internal_loop[3] = (double)(fv3);
	array_ninio_internal_loop[4] = (double)(fv4);
	
	/* multi-loop */
	if (fscanf(f_data,"%f%f%f",&fv1,&fv2,&fv3) != 3) { /* read failed */
		fprintf(stderr,"miscloop: corrupted file (fgets 4)%s\n",getJobID());
		exit(EXIT_FAILURE);
	}
	offset_multi_loop = (double)(fv1);
	free_base_penalty_multi_loop = (double)(fv2);
	helix_penalty_multi_loop = (double)(fv3);
	
	/* efn2 multi-loop */
	if (fscanf(f_data,"%f%f%f",&fv1,&fv2,&fv3) != 3) { /* read failed */
		fprintf(stderr,"miscloop: corrupted file (fgets 5)%s\n",getJobID());
		exit(EXIT_FAILURE);
	}
	offset_efn2_multi_loop = (double)(fv1);
	free_base_penalty_efn2_multi_loop = (double)(fv2);
	helix_penalty_efn2_multi_loop = (double)(fv3);
	
	/* other */
	if (fscanf(f_data,"%f",&fv1) != 1) { /* read failed */
		fprintf(stderr,"miscloop: corrupted file (fgets 6)%s\n",getJobID());
		exit(EXIT_FAILURE);
	}
	au_penalty = (double)(fv1);
	
	
	if (fscanf(f_data,"%f",&fv1) != 1) { /* read failed */
		fprintf(stderr,"miscloop: corrupted file (fgets 7)%s\n",getJobID());
		exit(EXIT_FAILURE);
	}
	ggg_hairpin_bonus = (double)(fv1);
	
	if (fscanf(f_data,"%f",&fv1) != 1) { /* read failed */
		fprintf(stderr,"miscloop: corrupted file (fgets 8)%s\n",getJobID());
		exit(EXIT_FAILURE);
	}
	c_hairpin_slope = (double)(fv1);
	
	if (fscanf(f_data,"%f",&fv1) != 1) { /* read failed */
		fprintf(stderr,"miscloop: corrupted file (fgets 9)%s\n",getJobID());
		exit(EXIT_FAILURE);
	}
	c_hairpin_intersect = (double)(fv1);
	
	if (fscanf(f_data,"%f",&fv1) != 1) { /* read failed */
		fprintf(stderr,"miscloop: corrupted file (fgets 10)%s\n",getJobID());
		exit(EXIT_FAILURE);
	}
	c_hairpin_3 = (double)(fv1);
	
	if (fscanf(f_data,"%f",&fv1) != 1) { /* read failed */
		fprintf(stderr,"miscloop: corrupted file (fgets 11)%s\n",getJobID());
		exit(EXIT_FAILURE);
	}
	inter_molecular_init_energy = (double)(fv1);

	
	if (fscanf(f_data,"%d",&gail_rule) != 1) { /* read failed */
		fprintf(stderr,"miscloop: corrupted file (fgets 12)%s\n",getJobID());
		exit(EXIT_FAILURE);
	}
	
	{ /* final checking */
		char buff[20];
		if (fscanf(f_data,"%s",&buff[0])==1) {
			fprintf(stderr,"miscloop: Corrupted file (read %s at the expected EOF)%s\n",buff,getJobID());
			exit(EXIT_FAILURE);
		}
	}
	
	fclose(f_data);
}
