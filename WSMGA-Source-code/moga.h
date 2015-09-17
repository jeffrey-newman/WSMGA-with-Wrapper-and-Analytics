#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "epanet2.h"

#ifndef var
#define var
#define square(x)  ((x)*(x))
#define SIGMA 1e-14
#define INF   1e14
#define PI    3.1415926
#define maxpop 2000  /*max population*/
#define maxchrom1 200  /*max chromosome length for integer chromosomes, which
						is also the max number of integer decision variables*/
#define maxno_option 100  /*max no. of options for each decision variable*/
#define maxchrom2  100  /*max chromosome length for real chromosomes, which is 
                        also the max number of real decision variables*/
#define maxobj  10 /*max number of objective functions*/
#define maxcons 10 /*max number of constraints functions*/
#define maxpro 55 /*max number of properties recorded*/

/*define global variables*/

typedef struct  /*individual properties*/
{
	int chrom1[maxchrom1];  /*integer chromosomes*/
	double xchrom1[maxchrom1],
	      ychrom1[maxchrom1],
		  zchrom1[maxchrom1],
		  uchrom1[maxchrom1],
		  vchrom1[maxchrom1];
	double chrom2[maxchrom2]; /*real-coded chromosomes*/
	double fit[maxobj],  /*fitness values of the individual*/
		cons[maxcons],  /*list of constraint violations of the individual*/
		property[maxpro],  /*list of properties (computed values) of the individual*/
	    cub_len,  /*crowding distance of the individual*/
	    overallcons;  /*overall constraint violations for the individual*/
	int rank, /*rank of the individual*/
		flag; /*flag for ranking*/
	int tag, //tag to track if the ind is modified
        eval, //no. of evaluation that the ind is generated
        gen; //no.of generation that the ind is generated	
} individual;  /*structure defining individual*/

typedef struct  /*population properties*/
{
	int maxrank;  /*max no. of ranks present in the population*/
	int rankno[maxpop];  /*individuals at different ranks*/
	individual ind[maxpop],  /*different individuals*/
		*ind_ptr;
} population;  /*population structure*/

int noeval; //no. of evaluations

#endif
	
